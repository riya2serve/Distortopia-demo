#!/usr/bin/env python
"""CLI for disto filter."""

from __future__ import annotations
import argparse
from pathlib import Path


def _setup_filter_subparser(subparser, description: str):
    p = subparser.add_parser(
        "filter",
        description=description,
        help="Filter mapcall VCF to biallelic SNPs + heterozygous sites (parallel per contig).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=True,
    )

    p.add_argument("--vcf", required=True, type=Path, help="Input VCF (.vcf.gz) from disto mapcall")
    p.add_argument("--bam", required=False, type=Path, help="Optional BAM (.bam) to derive contigs (recommended)")
    p.add_argument("-o", "--out", required=True, type=Path, help="Output directory")
    p.add_argument("-p", "--prefix", default="", type=str, help="Output prefix (default: derived from input)")
    p.add_argument("-t", "--threads", default=1, type=int, help="Workers for per-contig multiprocessing")
    p.add_argument("--keep-tmp", action="store_true", help="Keep per-contig temp files")

    return p

