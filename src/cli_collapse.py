#!/usr/bin/env python

import argparse
from pathlib import Path
from .make_wide import make_wide


def _setup_collapse_subparser(subparsers: argparse._SubParsersAction, header: str = None) -> None:
    """Setup the `collapse` subcommand parser."""
    parser = subparsers.add_parser(
        "collapse",
        description=header,
        help="Collapse read-level crossovers into locus-level events",
        formatter_class=make_wide(argparse.RawDescriptionHelpFormatter),
    )

    parser.add_argument(
        "-i", "--infer-tsv", metavar="Path", type=Path, required=True,
        help="Input TSV from `disto infer` (read-level crossovers).",
    )
    parser.add_argument(
        "-o", "--out", metavar="Path", type=Path, required=True,
        help="Output locus TSV (collapsed crossover loci).",
    )
    parser.add_argument(
        "-g", "--merge-gap", metavar="int", type=int, default=0,
        help="Merge loci if within this many bp. [default=0]",
    )
    parser.add_argument(
        "--keep-na", action="store_true",
        help="Keep NA crossover rows (nonrecombinant reads). [default: drop NA rows]",
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

