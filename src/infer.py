#!/usr/bin/env python

"""Generate TSV with inferred crossover window for each recombinant read

TODO: 
    - parallelize                      # done via samtools -@ <threads>
    - mask/exclude near chrom ends     # done via edge_mask_bp and reference .fai
"""

from __future__ import annotations
from typing import Dict, List, Tuple
from pathlib import Path
import gzip
import subprocess as sp
import pandas as pd
from loguru import logger


# CIGAR ops that consume reference/query
_REF_CONSUME = set("MDN=X")
_QUERY_CONSUME = set("MIS=X")


def load_phased_biallelic_snps(vcf_gz: Path) -> Dict[str, List[Tuple[int, str, str]]]:
    """Parse a VCF to collect phased, heterozygous, biallelic SNPs for the first sample.

    Returns
    -------
    dict[str, list[tuple[int, str, str]]]
        {chrom: [(pos1, allele0_base, allele1_base), ...]} sorted by pos.
    """
    phased: Dict[str, List[Tuple[int, str, str]]] = {}

    with gzip.open(vcf_gz, 'rt') as fh:
        gt_idx = None  # index of GT within FORMAT
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                # header columns: #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [SAMPLE...]
                header = line.rstrip("\n").split("\t")
                # require at least one sample
                if len(header) < 10:
                    return phased
                continue

            # Parse core VCF columns
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10:
                continue

            chrom, pos_s, _id, ref, alt, _qual, _filter, _info, fmt = cols[:9]
            sample0 = cols[9]

            # biallelic SNP only
            alt_list = alt.split(",")
            if len(ref) != 1 or len(alt_list) != 1 or len(alt_list[0]) != 1:
                continue

            # Find GT field index in FORMAT
            fmt_keys = fmt.split(":")
            if gt_idx is None:
                try:
                    gt_idx = fmt_keys.index("GT")
                except ValueError:
                    # No GT field at all
                    continue

            samp_fields = sample0.split(":")
            if gt_idx >= len(samp_fields):
                continue

            gt = samp_fields[gt_idx]
            # Require phased (uses '|'), non-missing, and simple diploid like "0|1" or "1|0"
            if "|" not in gt or "." in gt:
                continue
            a_str0, a_str1 = gt.split("|", 1)
            # Ensure both alleles are 0/1 and heterozygous
            if a_str0 not in ("0", "1") or a_str1 not in ("0", "1"):
                continue
            a0, a1 = int(a_str0), int(a_str1)
            if a0 == a1:
                continue

            pos = int(pos_s)
            alleles = [ref, alt_list[0]]
            # Map genotype indexes to actual bases (already ref-oriented)
            phased.setdefault(chrom, []).append((pos, alleles[a0], alleles[a1]))

    # Sort positions within each chromosome
    for c in list(phased.keys()):
        phased[c].sort(key=lambda x: x[0])

    return phased


def _parse_cigar(cigar: str) -> List[Tuple[str, int]]:
    """Return [('M',len), ('I',len), ...] from a CIGAR string."""
    ops: List[Tuple[str, int]] = []
    n = 0
    for ch in cigar:
        if ch.isdigit():
            n = n * 10 + ord(ch) - 48
        else:
            # ch is an op
            ops.append((ch, n))
            n = 0
    return ops


def _ref_end_from_cigar(pos1: int, cigar_ops: List[Tuple[str, int]]) -> int:
    """Compute 1-based inclusive reference end (matches your original)."""
    ref_len = sum(length for op, length in cigar_ops if op in _REF_CONSUME)
    # pos1 is 1-based leftmost; end inclusive is pos1 + ref_len - 1
    return pos1 + ref_len - 1


def _get_chrom_lengths(reference: Path) -> Dict[str, int]:
    """Get chromosome lengths from reference.fai, creating it with samtools faidx if needed."""
    fai = Path(str(reference) + ".fai")
    if not fai.exists():
        logger.info(f"Indexing reference with samtools faidx: {reference}")
        rc = sp.call(["samtools", "faidx", str(reference)])
        if rc != 0:
            raise RuntimeError("samtools faidx failed")

    lengths: Dict[str, int] = {}
    with fai.open() as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            chrom, length_s = fields[0], fields[1]
            lengths[chrom] = int(length_s)
    return lengths


def read_crossovers(
    bam_path: Path,
    phased_dict: Dict[str, List[Tuple[int,str,str]]],
    min_snps: int,
    chrom_lengths: Dict[str, int] | None = None,
    edge_mask_bp: int = 0,
    threads: int = 1,
) -> pd.DataFrame:
    """Read BAM without pysam; uses `samtools view -h -F 0x904` and parses CIGAR.

    Parameters
    ----------
    bam_path : Path
        Input BAM/CRAM.
    phased_dict : Dict[str, List[Tuple[int, str, str]]]
        {chrom: [(pos1, allele0, allele1), ...]} with pos1 1-based.
    min_snps : int
        Minimum phased SNPs observed on a read.
    chrom_lengths : dict, optional
        {chrom: length} from reference .fai (for edge masking).
    edge_mask_bp : int, optional
        Mask crossovers within this many bp of chrom ends. 0 = no masking.
    threads : int, optional
        Number of threads to give to samtools view (-@). Default 1.
    """
    # Stream SAM records; -h keeps headers so we can skip them; -F 0x904 filters unmapped (0x4),
    # secondary (0x100), supplementary (0x800) early for speed.
    cmd = ["samtools", "view", "-h", "-F", "0x904"]
    if threads and threads > 1:
        cmd.extend(["-@", str(threads)])
    cmd.append(str(bam_path))

    proc = sp.Popen(cmd, stdout=sp.PIPE, text=True, bufsize=1)

    rows = []
    try:
        assert proc.stdout is not None
        for line in proc.stdout:
            # skip headers
            if line.startswith("@"):
                continue

            # Parse SAM fields (minimum 11)
            # QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL [TAG...]
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 11:
                continue  # malformed
            qname = parts[0]
            rname = parts[2]
            pos1 = int(parts[3])          # 1-based leftmost
            cigar = parts[5]
            seq = parts[9]

            # Extra check (already filtered via -F 0x904)
            if rname == "*" or cigar == "*" or rname not in phased_dict:
                continue

            # parse cigar string to List[Tuple[str, int]]
            cigar_ops = _parse_cigar(cigar)
            ref_end1 = _ref_end_from_cigar(pos1, cigar_ops)  # 1-based inclusive

            # Build ref->base mapping for *matches only* (M,=,X) like get_aligned_pairs(matches_only=True)
            ref2base: Dict[int, str] = {}
            rpos = pos1                   # 1-based reference cursor
            qpos = 0                      # 0-based read cursor (in SEQ as stored in SAM)

            for op, length in cigar_ops:
                if op in ("M", "=", "X"):
                    # These consume both ref and query; record per-base mapping
                    for _ in range(length):
                        ref2base[rpos] = seq[qpos]
                        rpos += 1
                        qpos += 1
                elif op == "I":          # insertion relative to ref (consume query only)
                    qpos += length
                elif op in ("D", "N"):   # deletion or ref skip (consume ref only)
                    rpos += length
                elif op == "S":          # soft clip (consume query only)
                    qpos += length
                elif op in ("H", "P"):   # hard clip / pad (consume neither)
                    pass
                else:
                    # Unknown op; be conservative and skip
                    pass

            # Scan phased SNPs in read span [pos1, ref_end1]
            snp_pos: List[int] = []
            bits: List[int] = []
            for (snp_pos1, a0, a1) in phased_dict[rname]:
                if snp_pos1 < pos1 or snp_pos1 > ref_end1:
                    continue
                b = ref2base.get(snp_pos1)
                if b == a0:
                    snp_pos.append(snp_pos1)
                    bits.append(0)
                elif b == a1:
                    snp_pos.append(snp_pos1)
                    bits.append(1)
                # else: not informative/mismatch -> ignore

            if len(bits) < min_snps:
                continue

            # Detect crossover(s) as transitions in phased bits
            changes = [i for i in range(1, len(bits)) if bits[i] != bits[i - 1]]
            phased_str = "".join(map(str, bits))

            if len(changes) == 0:
                xl, xr = "NA", "NA"
            elif len(changes) == 1:
                i = changes[0]
                xl, xr = snp_pos[i - 1], snp_pos[i]
            else:
                # >1 transition: skip (consistent with your original)
                continue

            # Optional: mask crossovers near chromosome ends
            if (
                xl != "NA"
                and xr != "NA"
                and chrom_lengths is not None
                and edge_mask_bp > 0
            ):
                chrom_len = chrom_lengths.get(rname)
                if chrom_len is not None:
                    left_min = edge_mask_bp
                    right_max = chrom_len - edge_mask_bp
                    # If crossover window is outside safe interior region, skip
                    if not (left_min <= xl <= right_max and left_min <= xr <= right_max):
                        continue

            rows.append(
                {
                    "scaff": rname,
                    "start": pos1,
                    "end": ref_end1,
                    "nsnps": len(bits),
                    "phased_snps": phased_str,
                    "crossover_left": xl,
                    "crossover_right": xr,
                    "read": qname,
                }
            )
    finally:
        if proc.stdout:
            proc.stdout.close()
        rc = proc.wait()
        if rc != 0:
            raise RuntimeError(f"samtools view exited with code {rc}")

    return pd.DataFrame(
        rows,
        columns=[
            "scaff",
            "start",
            "end",
            "nsnps",
            "phased_snps",
            "crossover_left",
            "crossover_right",
            "read",
        ],
    )


def run_infer(
    bam_path: Path,
    vcf_gz: Path,
    reference: Path,
    outdir: Path,
    prefix: str,
    min_snps: int,
    threads: int = 1,
    edge_mask_bp: int = 0,
):
    """Infer crossover windows and write TSV.

    Parameters
    ----------
    bam_path : Path
        Input BAM with long-read alignments.
    vcf_gz : Path
        Phased, biallelic SNP VCF (bgzipped).
    reference : Path
        Reference FASTA used for alignment (for chrom lengths).
    outdir : Path
        Output directory (created if needed).
    prefix : str
        Output prefix, e.g. 'Chr3' or 'test1'.
    min_snps : int
        Minimum informative phased SNPs per read.
    threads : int, optional
        Number of threads to pass to samtools view (-@). Default 1.
    edge_mask_bp : int, optional
        Number of bp to mask at each chromosome end when considering crossovers.
        If 0, no masking is applied.
    """
    logger.info("inferring crossover positions")

    outdir = outdir.expanduser().absolute()
    outdir.mkdir(parents=True, exist_ok=True)

    bam_path = bam_path.expanduser().absolute()
    vcf_gz = vcf_gz.expanduser().absolute()
    reference = reference.expanduser().absolute()

    phased_dict = load_phased_biallelic_snps(vcf_gz)
    chrom_lengths = _get_chrom_lengths(reference)

    data = read_crossovers(
        bam_path=bam_path,
        phased_dict=phased_dict,
        min_snps=min_snps,
        chrom_lengths=chrom_lengths,
        edge_mask_bp=edge_mask_bp,
        threads=threads,
    )

    outpath = outdir / f"{prefix}.tsv"
    data.to_csv(outpath, sep="\t", index=False)
    print(data)
    logger.info(f"Wrote crossover table: {outpath}")

