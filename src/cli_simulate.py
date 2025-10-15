#!/usr/bin/env python

import argparse
from pathlib import Path
from .make_wide import make_wide

def _setup_simulate_subparser(subparsers: argparse._SubParsersAction, header: str = None) -> None:
    """Setup and return an ArgumentParser w/ subcommands."""
    parser = subparsers.add_parser(
        "simulate",
        description=header,
        help="Simulate recombinant gametes",
        formatter_class=make_wide(argparse.RawDescriptionHelpFormatter),
    )
    parser.add_argument(
        "-r", "--reference", metavar="Path", type=Path, required=True,
        help="Path to the reference genome fasta.",
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
        "-e", "--heterozygosity", metavar="float", type=float, default=1e-3,
        help="heterozygosity (per-site divergence between genome haplotypes).",
    )
    parser.add_argument(
        "-x", "--crossover-rate", metavar="float", type=float, default=1e-8,
        help="crossover rate per site per chrom per gamete.",
    )
    parser.add_argument(
        "-c", "--chromosomes", metavar="str", type=str, nargs="*",
        help="subselect one or more chromosomes from reference.",
    )    
    # parser.add_argument(
    #     "-u", "--uniform", action="store_true",
    #     help="model recomb rate as uniform, versus sinusoid.",
    # )
    parser.add_argument(
        "-n", "--nreads", metavar="int", type=int, default=10_000_000,
        help="Number of long reads to simulate. [default=1_000_000]",
    )
    parser.add_argument(
        "-g", "--read-length", metavar="int", type=int, default=100_000,
        help="Long read lengths [default=100_000]",
    )
    parser.add_argument(
        "-u", "--interference", metavar="float", type=float, default=5.0,
        help="Strength of crossover interference [default=5.0]",
    )    
    parser.add_argument(
        "-s", "--random-seed", metavar="int", type=int, default=None,
        help="Random seed [default=None]",
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
