#!/usr/bin/env python
"""Filter mapcall VCF for biallelic SNPs and heterozygous sites only.

Runs after:  disto mapcall
Runs before: disto infer

Strategy:
  - Determine contigs (if BAM provided use BAM idxstats)
  - For each contig, run:
      bcftools view -r contig -v snps -m2 -M2 -i 'GT=het' -Oz -o tmp/contig.vcf.gz
      bcftools index tmp/contig.vcf.gz
    in parallel (multiprocessing)
  - Concatenate contig VCFs in order and index final output
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
import subprocess as sp
from multiprocessing import Pool
from loguru import logger


BIN = Path(sys.prefix) / "bin"
BIN_SAM = str(BIN / "samtools")
BIN_BCF = str(BIN / "bcftools")


_HET_EXPR = 'GT="0/1" || GT="1/0" || GT="0|1" || GT="1|0"'


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
    """Fallback: parse ##contig=<ID=...> from VCF."""
    out = sp.check_output([BIN_BCF, "view", "-h", str(vcf_gz)], text=True)
    contigs: list[str] = []
    for line in out.splitlines():
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
    return contigs


def _filter_contig_worker(args: tuple[str, str, str]) -> str:
    """Worker to filter one contig at a time and write tmp/contig.vcf.gz (+ index)."""
    contig, vcf_in, tmpdir = args
    tmpdir_p = Path(tmpdir)

    out_vcf = tmpdir_p / f"{contig}.vcf.gz"
    out_tmp = Path(str(out_vcf) + ".tmp")

    # Filter: SNPs only, biallelic only, heterozygous only
    cmd = [
        BIN_BCF, "view",
        "-r", contig,
        "-v", "snps",
        "-m2", "-M2",
        "-i", _HET_EXPR,
        "-Oz",
        "-o", str(out_tmp),
        str(vcf_in),
    ]
    sp.run(cmd, check=True)

    # Validate compressed stream, then publish
    sp.run(["gzip", "-t", str(out_tmp)], check=True)
    os.replace(out_tmp, out_vcf)

    # Index per-contig
    sp.run([BIN_BCF, "index", "-f", str(out_vcf)], check=True)
    return str(out_vcf)


def filter_vcf_per_contig(
    vcf_gz: Path,
    out_vcf_gz: Path,
    contigs: list[str],
    threads: int,
    tmpdir: Path,
) -> Path:
    """Run per-contig filtering in parallel, then concat into out_vcf_gz."""
    threads = max(1, int(threads))
    tmpdir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Filtering {len(contigs)} contigs with {threads} workers")
    work = [(c, str(vcf_gz), str(tmpdir)) for c in contigs]

    with Pool(processes=threads) as pool:
        pool.map(_filter_contig_worker, work)

    # Ensure concat is in contig order
    part_files = [str(tmpdir / f"{c}.vcf.gz") for c in contigs if (tmpdir / f"{c}.vcf.gz").exists()]
    if not part_files:
        raise RuntimeError("No per-contig VCFs produced. Check bcftools output/errors.")

    out_tmp = Path(str(out_vcf_gz) + ".tmp")

    cmd_concat = [BIN_BCF, "concat", "-Oz", "-o", str(out_tmp), *part_files]
    sp.run(cmd_concat, check=True)

    sp.run(["gzip", "-t", str(out_tmp)], check=True)
    os.replace(out_tmp, out_vcf_gz)

    sp.run([BIN_BCF, "index", "-f", str(out_vcf_gz)], check=True)
    logger.info(f"Wrote filtered VCF: {out_vcf_gz}")
    return out_vcf_gz


def run_filter(
    vcf_gz: Path,
    outdir: Path,
    prefix: str,
    threads: int,
    bam: Path | None = None,
    keep_tmp: bool = False,
) -> Path:
    """Main entry: filter mapcall VCF -> biallelic SNPs + het only, parallel per contig."""
    outdir = outdir.expanduser().absolute()
    outdir.mkdir(parents=True, exist_ok=True)

    vcf_gz = vcf_gz.expanduser().absolute()
    if bam is not None:
        bam = bam.expanduser().absolute()

    base = outdir / (prefix if prefix else vcf_gz.stem.replace(".vcf", ""))

    out_vcf_gz = base.with_suffix(".biallelic_het.vcf.gz")
    tmpdir = outdir / f".{base.name}.filter_tmp"

    # Determine contigs
    if bam is not None:
        logger.info("deriving contigs from BAM idxstats")
        contigs = contigs_from_bam(bam)
    else:
        logger.info("deriving contigs from VCF")
        contigs = contigs_from_vcf(vcf_gz)

    if not contigs:
        raise RuntimeError("Could not determine contigs. Provide --bam or ensure VCF header has ##contig lines.")

    # Filter + concat
    filter_vcf_per_contig(vcf_gz, out_vcf_gz, contigs, threads, tmpdir)

    # Cleanup
    if not keep_tmp:
        for p in tmpdir.glob("*"):
            try:
                p.unlink()
            except Exception:
                pass
        try:
            tmpdir.rmdir()
        except Exception:
            pass

    return out_vcf_gz
