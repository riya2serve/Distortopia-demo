#!/usr/bin/env python
"""
Generate TSV with inferred crossover window for each recombinant read.

Parallel per-contig version (mimics filter.py):
  - Determine contigs (prefer BAM idxstats; fallback to VCF header)
  - For each contig, run samtools view restricted to that contig
  - Infer crossover windows per read and write tmp/<contig>.tsv
  - Concatenate tmp TSVs in contig order into final <prefix>.tsv

Notes:
  - Uses conda-env samtools via Path(sys.prefix)/bin/samtools
  - threads parameter is used as number of worker processes (per contig)
"""

from __future__ import annotations

from typing import Dict, List, Tuple, Optional
from pathlib import Path
import os
import sys
import gzip
import subprocess as sp
from multiprocessing import Pool
from loguru import logger


# conda env binaries (like filter.py)
BIN = Path(sys.prefix) / "bin"
BIN_SAM = str(BIN / "samtools")

# CIGAR ops that consume reference/query
_REF_CONSUME = set("MDN=X")
_QUERY_CONSUME = set("MIS=X")  # kept for completeness; not used below

OUT_COLUMNS = [
    "scaff",
    "start",
    "end",
    "nsnps",
    "phased_snps",
    "crossover_left",
    "crossover_right",
    "read",
]


def contigs_from_bam(bam: Path) -> list[str]:
    """Use samtools idxstats to get contig list from BAM."""
    out = sp.check_output([BIN_SAM, "idxstats", str(bam)], text=True)
    contigs: list[str] = []
    for line in out.splitlines():
        name = line.split("\t", 1)[0]
        if name and name != "*":
            contigs.append(name)
    return contigs


def contigs_from_vcf(vcf_gz: Path) -> list[str]:
    """Fallback: parse ##contig=<ID=...> from VCF header (using gzip; no bcftools required)."""
    contigs: list[str] = []
    with gzip.open(vcf_gz, "rt") as fh:
        for line in fh:
            if line.startswith("##contig=") and "ID=" in line:
                try:
                    inside = line.split("<", 1)[1].rsplit(">", 1)[0]
                    parts = {}
                    for p in inside.split(","):
                        if "=" in p:
                            k, v = p.split("=", 1)
                            parts[k] = v
                    cid = parts.get("ID")
                    if cid:
                        contigs.append(cid)
                except Exception:
                    continue
            elif line.startswith("#CHROM"):
                break
    return contigs


def load_phased_biallelic_snps(vcf_gz: Path) -> Dict[str, List[Tuple[int, str, str]]]:
    phased: Dict[str, List[Tuple[int, str, str]]] = {}

    with gzip.open(vcf_gz, "rt") as fh:
        gt_idx = None
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10:
                continue

            chrom, pos_s, _id, ref, alt, *_rest, fmt = cols[:9]
            sample0 = cols[9]

            alt_list = alt.split(",")
            if len(ref) != 1 or len(alt_list) != 1 or len(alt_list[0]) != 1:
                continue

            fmt_keys = fmt.split(":")
            if gt_idx is None:
                try:
                    gt_idx = fmt_keys.index("GT")
                except ValueError:
                    continue

            samp_fields = sample0.split(":")
            if gt_idx >= len(samp_fields):
                continue

            gt = samp_fields[gt_idx]
            if "|" not in gt or "." in gt:
                continue

            a0, a1 = gt.split("|")
            if a0 == a1:
                continue

            pos = int(pos_s)
            alleles = [ref, alt_list[0]]
            phased.setdefault(chrom, []).append((pos, alleles[int(a0)], alleles[int(a1)]))

    for c in phased:
        phased[c].sort(key=lambda x: x[0])

    return phased


def _parse_cigar(cigar: str) -> List[Tuple[str, int]]:
    ops: List[Tuple[str, int]] = []
    n = 0
    for ch in cigar:
        if ch.isdigit():
            n = n * 10 + (ord(ch) - 48)
        else:
            ops.append((ch, n))
            n = 0
    return ops


def _ref_end_from_cigar(pos1: int, cigar_ops: List[Tuple[str, int]]) -> int:
    ref_len = sum(length for op, length in cigar_ops if op in _REF_CONSUME)
    return pos1 + ref_len - 1


def _get_chrom_lengths(reference: Path) -> Dict[str, int]:
    fai = Path(str(reference) + ".fai")
    if not fai.exists():
        sp.check_call([BIN_SAM, "faidx", str(reference)])

    lengths: Dict[str, int] = {}
    with fai.open() as fh:
        for line in fh:
            chrom, length_s, *_ = line.rstrip("\n").split("\t")
            lengths[chrom] = int(length_s)
    return lengths


def _infer_contig_worker(args: tuple[str, str, str, list[tuple[int, str, str]], int, int, int, str]) -> str:
    """
    Worker: infer crossovers on one contig and write tmp/<contig>.tsv.

    args:
      contig, bam_path, contig_name, phased_sites, chrom_len, min_snps, edge_mask_bp, tmpdir
    """
    contig, bam_path, contig_name, phased_sites, chrom_len, min_snps, edge_mask_bp, tmpdir = args

    tmpdir_p = Path(tmpdir)
    tmpdir_p.mkdir(parents=True, exist_ok=True)  # safety (especially under multiprocessing)
    out_tsv = tmpdir_p / f"{contig}.tsv"
    out_tmp = Path(str(out_tsv) + f".tmp.{os.getpid()}")

    # IMPORTANT: always create the temp file so os.replace cannot fail later
    # even if we end up writing zero rows.
    out_tmp.write_text("")

    # If no phased SNPs on this contig, publish empty file
    if not phased_sites:
        os.replace(out_tmp, out_tsv)
        return str(out_tsv)

    # samtools view for just this contig; keep primary alignments only (-F 0x904)
    cmd = [BIN_SAM, "view", "-h", "-F", "0x904", str(bam_path), contig_name]
    proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, text=True, bufsize=1)

    # phased_sites is sorted by position
    snp_i = 0  # advances monotonically because BAM is coordinate-sorted per contig

    # write buffered lines (TSV without header; main process adds header once)
    lines: List[str] = []

    try:
        assert proc.stdout is not None
        for line in proc.stdout:
            if line.startswith("@"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 11:
                continue

            qname = parts[0]
            rname = parts[2]
            if rname == "*" or rname != contig_name:
                continue

            pos1 = int(parts[3])
            cigar = parts[5]
            seq = parts[9]
            if cigar == "*" or seq == "*":
                continue

            cigar_ops = _parse_cigar(cigar)
            ref_end1 = _ref_end_from_cigar(pos1, cigar_ops)
            if ref_end1 < pos1:
                continue

            # advance snp_i to first SNP >= pos1 (monotonic)
            while snp_i < len(phased_sites) and phased_sites[snp_i][0] < pos1:
                snp_i += 1
            if snp_i >= len(phased_sites) or phased_sites[snp_i][0] > ref_end1:
                continue

            # walk CIGAR and “jump” to SNP positions (avoid building a full ref->base map)
            snp_pos: List[int] = []
            bits: List[int] = []

            j = snp_i  # local SNP pointer for this read
            rpos = pos1
            qpos = 0

            for op, length in cigar_ops:
                if length <= 0:
                    continue

                if op in ("M", "=", "X"):
                    block_start = rpos
                    block_end = rpos + length  # half-open [start, end)

                    # consume SNPs that fall in this match block
                    while j < len(phased_sites):
                        p, a0, a1 = phased_sites[j]
                        if p < block_start:
                            j += 1
                            continue
                        if p >= block_end:
                            break

                        offset = p - block_start
                        qi = qpos + offset
                        if 0 <= qi < len(seq):
                            b = seq[qi]
                            if b == a0:
                                snp_pos.append(p)
                                bits.append(0)
                            elif b == a1:
                                snp_pos.append(p)
                                bits.append(1)
                        j += 1

                    rpos += length
                    qpos += length

                elif op == "I":
                    qpos += length
                elif op in ("D", "N"):
                    rpos += length
                elif op == "S":
                    qpos += length
                else:
                    continue

                # early exit if we passed the end of the read region with respect to SNPs
                if j < len(phased_sites) and phased_sites[j][0] > ref_end1:
                    break

            if len(bits) < min_snps:
                continue

            changes = [i for i in range(1, len(bits)) if bits[i] != bits[i - 1]]
            if len(changes) > 1:
                continue

            if not changes:
                xl = xr = "NA"
            else:
                i = changes[0]
                xl, xr = snp_pos[i - 1], snp_pos[i]

            if xl != "NA" and edge_mask_bp > 0:
                # filter crossovers too close to chrom ends
                if not (edge_mask_bp <= int(xl) <= chrom_len - edge_mask_bp):
                    continue

            lines.append(
                f"{contig_name}\t{pos1}\t{ref_end1}\t{len(bits)}\t{''.join(map(str, bits))}\t{xl}\t{xr}\t{qname}\n"
            )

            # flush occasionally to avoid huge memory per worker
            if len(lines) >= 200_000:
                with out_tmp.open("a") as oh:
                    oh.writelines(lines)
                lines.clear()

    finally:
        # close stdout so proc can terminate cleanly
        if proc.stdout is not None:
            proc.stdout.close()
        _, stderr = proc.communicate()
        if proc.returncode != 0:
            # publish stderr to log before raising
            logger.error(f"samtools view failed on contig {contig_name} (rc={proc.returncode}). stderr:\n{stderr}")
            raise RuntimeError(f"samtools view failed on contig {contig_name} (rc={proc.returncode})")

    # final flush and publish atomically
    if lines:
        with out_tmp.open("a") as oh:
            oh.writelines(lines)
        lines.clear()

    os.replace(out_tmp, out_tsv)
    return str(out_tsv)


def infer_per_contig_to_tsv(
    bam_path: Path,
    phased_dict: Dict[str, List[Tuple[int, str, str]]],
    contigs: list[str],
    chrom_lengths: Dict[str, int],
    min_snps: int,
    edge_mask_bp: int,
    threads: int,
    out_tsv: Path,
    tmpdir: Path,
) -> Path:
    """Parallel per-contig infer, then concatenate tmp TSVs into out_tsv."""
    threads = max(1, int(threads))
    tmpdir.mkdir(parents=True, exist_ok=True)

    work: list[tuple[str, str, str, list[tuple[int, str, str]], int, int, int, str]] = []
    for c in contigs:
        phased_sites = phased_dict.get(c, [])
        clen = chrom_lengths.get(c, 0)
        work.append((c, str(bam_path), c, phased_sites, clen, int(min_snps), int(edge_mask_bp), str(tmpdir)))

    logger.info(f"Inferring crossovers on {len(contigs)} contigs with {threads} workers")
    with Pool(processes=threads) as pool:
        pool.map(_infer_contig_worker, work)

    # concatenate in contig order (header once)
    out_tmp = Path(str(out_tsv) + ".tmp")
    with out_tmp.open("w") as oh:
        oh.write("\t".join(OUT_COLUMNS) + "\n")
        for c in contigs:
            part = tmpdir / f"{c}.tsv"
            if not part.exists():
                continue
            with part.open("r") as fh:
                for ln in fh:
                    oh.write(ln)

    os.replace(out_tmp, out_tsv)
    logger.info(f"Wrote crossover table: {out_tsv}")
    return out_tsv


def run_infer(
    bam_path: Path,
    vcf_gz: Path,
    reference: Path,
    outdir: Path,
    prefix: str,
    min_snps: int,
    threads: int = 1,
    edge_mask_bp: int = 0,
    include_scaffs: Optional[List[str]] = None,
    chunk_size: int = 50_000,  # kept for API compatibility; not used in per-contig mode
):
    """Infer crossover windows using parallel per-contig processing (multiprocessing)."""
    logger.info("inferring crossover positions (parallel per contig)")

    outdir = outdir.expanduser().absolute()
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / f"{prefix}.tsv"
    tmpdir = outdir / f".{prefix}.infer_tmp"

    bam_path = bam_path.expanduser().absolute()
    vcf_gz = vcf_gz.expanduser().absolute()
    reference = reference.expanduser().absolute()

    phased_dict = load_phased_biallelic_snps(vcf_gz)
    chrom_lengths = _get_chrom_lengths(reference)

    # Determine contigs (prefer BAM)
    try:
        contigs = contigs_from_bam(bam_path)
    except Exception:
        logger.warning("Could not derive contigs from BAM idxstats; falling back to VCF header contigs")
        contigs = contigs_from_vcf(vcf_gz)

    if include_scaffs:
        keep = set(include_scaffs)
        contigs = [c for c in contigs if c in keep]

    # Only keep contigs that appear in reference lengths
    contigs = [c for c in contigs if c in chrom_lengths]
    if not contigs:
        raise RuntimeError("No contigs to process after filtering. Check BAM/VCF/reference consistency.")

    infer_per_contig_to_tsv(
        bam_path=bam_path,
        phased_dict=phased_dict,
        contigs=contigs,
        chrom_lengths=chrom_lengths,
        min_snps=int(min_snps),
        edge_mask_bp=int(edge_mask_bp),
        threads=int(threads),
        out_tsv=outpath,
        tmpdir=tmpdir,
    )

    # cleanup tmpdir
    for p in tmpdir.glob("*"):
        try:
            p.unlink()
        except Exception:
            pass
    try:
        tmpdir.rmdir()
    except Exception:
        pass
