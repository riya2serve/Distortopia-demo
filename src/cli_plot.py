#!/usr/bin/env python

import argparse
from pathlib import Path

from .make_wide import make_wide


def _setup_plot_subparser(subparsers: argparse._SubParsersAction, header: str = None) -> argparse.ArgumentParser:
    """Setup and return the `disto plot` subcommand parser."""
    parser = subparsers.add_parser(
        "plot",
        description=header,
        help="Plot crossover position distribution",
        formatter_class=make_wide(argparse.RawDescriptionHelpFormatter),
    )

    # --- Required / core ---
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
        help="Chromosome/contig to plot (can repeat: -c chr1 -c chr2). Default: all.",
    )
    parser.add_argument(
        "-p", "--prefix", metavar="str", type=str, default="test",
        help="Prefix for output files.",
    )

    # --- Binning controls ---
    parser.add_argument(
        "--bins",
        type=int,
        default=25,
        help="Histogram bin COUNT per chromosome (default: 25). Ignored if --bin-bp is set.",
    )
    parser.add_argument(
        "--bin-bp",
        type=int,
        default=None,
        help="Fixed genomic bin size in bp (e.g., 100000 for 100kb, 50000 for 50kb). "
             "If set, overrides --bins.",
    )

    # --- Error bars + rolling stats (read-level only) ---
    parser.add_argument(
        "--errorbars",
        action="store_true",
        help="Add SD error bars per bin (read-level infer TSV only). Requires --bin-bp.",
    )
    parser.add_argument(
        "--rolling-stats",
        action="store_true",
        help="Overlay rolling mean and rolling median across bins.",
    )
    parser.add_argument(
        "--rolling-window-bins",
        type=int,
        default=7,
        help="Window size (in bins) for rolling mean/median (centered). Default: 7.",
    )

    # --- Logging ---
    parser.add_argument(
        "-l", "--log-level", metavar="str", type=str, default="INFO",
        help="Log level (DEBUG, INFO, WARN, ERROR) [default=INFO]",
    )
    parser.add_argument(
        "-L", "--log-file", metavar="Path", type=Path,
        help="Log file. Logging to stdout is also appended to this file. [default=None].",
    )

    # NOTE: Validation like "--errorbars requires --bin-bp" is best done in the command handler
    # (where you call run_plot), because argparse subparser setup doesn't always have access to args.

    return parser
