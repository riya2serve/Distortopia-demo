#!/usr/bin/env python

"""Generate TSV with inferred crossover window for each recombinant read

Buffered pandas version:
- Uses pandas DataFrame chunks
- Uses DataFrame.to_csv(..., mode="a", chunksize=...)
- No csv module, no new imports

TODO:
    - parallelize                      # done via samtools -@ <threads>
    - mask/exclude near chrom ends     # done via edge_mask_bp and reference .fai
"""

from __future__ import annotations
from typing import Dict, List, Tuple, Iterator, Optional, Set
from pathlib import Path
import gzip
import subprocess as sp
import pandas as pd
from loguru import logger


# CIGAR ops that consume reference/query
_REF_CONSUME = set("MDN=X")
_QUERY_CONSUME = set("MIS=X")


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
    ops = []
    n = 0
    for ch in cigar:
        if ch.isdigit():
            n = n * 10 + ord(ch) - 48
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
        sp.check_call(["samtools", "faidx", str(reference)])

    lengths = {}
    with fai.open() as fh:
        for line in fh:
            chrom, length_s, *_ = line.rstrip("\n").split("\t")
            lengths[chrom] = int(length_s)
    return lengths


def read_crossovers(
    bam_path: Path,
    phased_dict: Dict[str, List[Tuple[int, str, str]]],
    min_snps: int,
    chrom_lengths: Dict[str, int],
    edge_mask_bp: int,
    threads: int,
    include_scaffs: Optional[Set[str]] = None,
) -> Iterator[Dict[str, object]]:
    """Yield crossover rows one at a time (pandas-friendly streaming)."""

    cmd = ["samtools", "view", "-h", "-F", "0x904"]
    if threads > 1:
        cmd.extend(["-@", str(threads)])
    cmd.append(str(bam_path))

    proc = sp.Popen(cmd, stdout=sp.PIPE, text=True, bufsize=1)

    try:
        for line in proc.stdout:
            if line.startswith("@"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 11:
                continue

            qname, rname, pos1, cigar, seq = (
                parts[0],
                parts[2],
                int(parts[3]),
                parts[5],
                parts[9],
            )

            if rname == "*" or rname not in phased_dict:
                continue
            if include_scaffs and rname not in include_scaffs:
                continue

            cigar_ops = _parse_cigar(cigar)
            ref_end1 = _ref_end_from_cigar(pos1, cigar_ops)

            ref2base = {}
            rpos, qpos = pos1, 0
            for op, length in cigar_ops:
                if op in ("M", "=", "X"):
                    for _ in range(length):
                        ref2base[rpos] = seq[qpos]
                        rpos += 1
                        qpos += 1
                elif op == "I":
                    qpos += length
                elif op in ("D", "N"):
                    rpos += length
                elif op == "S":
                    qpos += length

            snp_pos, bits = [], []
            for p, a0, a1 in phased_dict[rname]:
                if p < pos1:
                    continue
                if p > ref_end1:
                    break
                b = ref2base.get(p)
                if b == a0:
                    snp_pos.append(p)
                    bits.append(0)
                elif b == a1:
                    snp_pos.append(p)
                    bits.append(1)

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
                clen = chrom_lengths[rname]
                if not (edge_mask_bp <= xl <= clen - edge_mask_bp):
                    continue

            yield {
                "scaff": rname,
                "start": pos1,
                "end": ref_end1,
                "nsnps": len(bits),
                "phased_snps": "".join(map(str, bits)),
                "crossover_left": xl,
                "crossover_right": xr,
                "read": qname,
            }

    finally:
        proc.stdout.close()
        proc.wait()


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
    chunk_size: int = 50_000,
):
    """Infer crossover windows using buffered pandas writes."""

    logger.info("inferring crossover positions (buffered pandas)")

    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / f"{prefix}.tsv"

    phased_dict = load_phased_biallelic_snps(vcf_gz)
    chrom_lengths = _get_chrom_lengths(reference)

    buffer: List[Dict[str, object]] = []
    wrote_header = False

    for row in read_crossovers(
        bam_path=bam_path,
        phased_dict=phased_dict,
        min_snps=min_snps,
        chrom_lengths=chrom_lengths,
        edge_mask_bp=edge_mask_bp,
        threads=threads,
        include_scaffs=set(include_scaffs) if include_scaffs else None,
    ):
        buffer.append(row)

        if len(buffer) >= chunk_size:
            df = pd.DataFrame(buffer, columns=OUT_COLUMNS)
            df.to_csv(
                outpath,
                sep="\t",
                index=False,
                mode="w" if not wrote_header else "a",
                header=not wrote_header,
                chunksize=chunk_size,
            )
            wrote_header = True
            buffer.clear()

    if buffer:
        df = pd.DataFrame(buffer, columns=OUT_COLUMNS)
        df.to_csv(
            outpath,
            sep="\t",
            index=False,
            mode="w" if not wrote_header else "a",
            header=not wrote_header,
        )

    logger.info(f"Wrote crossover table: {outpath}")

