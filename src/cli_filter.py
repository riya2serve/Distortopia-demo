#!/usr/bin/env python
"""CLI for disto filter."""

from __future__ import annotations
import argparse
from pathlib import Path

from filter import run_filter


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--vcf", required=True, type=Path, help="Input VCF (.vcf.gz) from disto mapcall")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory")
    parser.add_argument("--prefix", default="", type=str, help="Output prefix (default: derived from input)")
    parser.add_argument("--threads", default=1, type=int, help="Workers for per-contig multiprocessing")
    parser.add_argument("--bam", default=None, type=Path, help="Optional BAM to derive contigs (recommended)")
    parser.add_argument("--keep-tmp", action="store_true", help="Keep per-contig temp VCFs")


def main(args: argparse.Namespace) -> Path:
    return run_filter(
        vcf_gz=args.vcf,
        outdir=args.outdir,
        prefix=args.prefix,
        threads=args.threads,
        bam=args.bam,
        keep_tmp=args.keep_tmp,
    )

