#!/usr/bin/env python

import argparse
from pathlib import Path
from .make_wide import make_wide


def _setup_plot_subparser(subparsers: argparse._SubParsersAction, header: str = None) -> None:
    """Setup and return an ArgumentParser w/ subcommands."""
    parser = subparsers.add_parser(
        "plot",
        description=header,
        help="Plot crossover position distribution",
        formatter_class=make_wide(argparse.RawDescriptionHelpFormatter),
    )
    parser.add_argument(
        "-t", "--tsv", metavar="Path", type=Path, required=True,
        help="Path to TSV file with crossover positions",
    )
    parser.add_argument(
        "-o", "--out", metavar="Path", type=Path, default=".",
        help="Output directory.",
    )
    parser.add_argument(
        "-c", "--chrom",
        action="append",
        default=None,
        help="Chromosome/contig to plot (can repeat: -c chr1 -c chr2). Default: all."
    )
    parser.add_argument(
        "--bins", 
        type=int, 
        default=25, 
        help="Histogram bins per chromosome (default: 25)"
    )
    parser.add_argument(
        "-p", "--prefix", metavar="str", type=str, default="test",
        help="Prefix for output files.",
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
