#!/usr/bin/env python

import argparse
from pathlib import Path
from .make_wide import make_wide


def _setup_infer_subparser(subparsers: argparse._SubParsersAction, header: str = None) -> None:
    """Setup and return an ArgumentParser w/ subcommands."""
    parser = subparsers.add_parser(
        "infer",
        description=header,
        help="Infer crossover map",
        formatter_class=make_wide(argparse.RawDescriptionHelpFormatter),
    )
    parser.add_argument(
        "-r", "--reference", metavar="Path", type=Path, required=True,
        help="Path to the reference genome fasta.",
    )
    parser.add_argument(
        "-v", "--vcf", metavar="Path", type=Path, required=True,
        help="Path to the phased VCF from disto mapcall.",
    )
    parser.add_argument(
        "-b", "--bam", metavar="Path", type=Path, required=True,
        help="Path to the sorted bam from disto mapcall.",
    )    
    parser.add_argument(
        "-o", "--out", metavar="Path", type=Path, default=".",
        help="Output directory.",
    )
    parser.add_argument(
        "-p", "--prefix", metavar="str", type=str, default="test",
        help="Prefix for output files.",
    )
    parser.add_argument(
        "-m", "--min-snps", metavar="int", type=int, default=5,
        help="Minimum num bi-allelic SNPs per read. [default=5]",
    )
    parser.add_argument(
        "-l", "--log-level", metavar="str", type=str, default="INFO",
        help="Log level (DEBUG, INFO, WARN, ERROR) [default=INFO]",
    )
    parser.add_argument(
        "-L", "--log-file", metavar="Path", type=Path,
        help="Log file. Logging to stdout is also appended to this file. [default=None]."
    )
    return parser
