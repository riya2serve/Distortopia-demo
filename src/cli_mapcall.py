#!/usr/bin/env python

import argparse
from pathlib import Path


def _setup_mapcall_subparser(subparsers: argparse._SubParsersAction, header: str = None) -> None:
    """Setup and return an ArgumentParser w/ subcommands."""
    parser = subparsers.add_parser(
        "mapcall",
        description=header,
        help="Call variants from long read gamete pool",
        formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),        
    )
    parser.add_argument(
        "-r", "--reference", metavar="Path", type=Path, required=True,
        help="Path to the reference genome fasta.",
    )
    parser.add_argument(
        "-g", "--gametes", metavar="Path", type=Path, required=True,
        help="Path to the pooled long read gamete fastq sequences.",
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
        "-q", "--min-map-q", metavar="int", type=int, default=10,
        help="Minimum alignment score. [default=10]",
    )
    parser.add_argument(
        "-Q", "--min-base-q", metavar="int", type=int, default=20,
        help="Minimum base quality score [default=20]",
    )
    parser.add_argument(
        "-t", "--threads", metavar="int", type=int, default=3,
        help="Run c/t multi-threaded jobs concurrently. Larger -t reduces RAM and I/O. [default=3]",
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
