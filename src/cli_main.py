#!/usr/bin/env python

import sys
import argparse
from loguru import logger

# CLI subparser setup functions
from .cli_simulate import _setup_simulate_subparser
from .cli_mapcall import _setup_mapcall_subparser
from .cli_filter import _setup_filter_subparser
from .cli_infer import _setup_infer_subparser
from .cli_collapse import _setup_collapse_subparser
from .cli_plot import _setup_plot_subparser

# Workhorses
from .simulate import run_simulate
from .mapcall import run_mapcall
from .filter import run_filter
from .infer import run_infer
from .collapse import run_collapse
from .plot import run_plot

"""
disto simulate -r REF.fa -p GAMETES -n 10_000_000 -l 100_000 -s 123
disto mapcall -r REF.fa -g GAMETES.fastq.gz -o . -q 10 -Q 20
disto filter --vcf GAMETES.vcf.gz --bam GAMETES.sorted.bam -o .
disto infer -r REF -v GAMETES.phased.vcf.gz -o . -p RATES
disto collapse -i RATES.tsv -o RATES.loci.tsv
disto plot -t RATES.tsv -o .
"""

VERSION = "0.0.1"
HEADER = f"""
-------------------------------------------------------------
disto [v.{VERSION}]
Infer crossover rate from sequenced gametes
-------------------------------------------------------------\
"""

DESCRIPTION = "disto command line tool. Select a positional subcommand:"


def setup_parsers() -> argparse.ArgumentParser:
    """Setup and return an ArgumentParser w/ subcommands."""
    parser = argparse.ArgumentParser(
        prog="disto",
        description=f"{HEADER}\n{DESCRIPTION}",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
    )
    parser.add_argument("-h", "--help", action="help", help=argparse.SUPPRESS)
    parser.add_argument("-v", "--version", action="version", version=f"disto {VERSION}")

    subparser = parser.add_subparsers(help="sub-commands", dest="subcommand")

    # add subcommands: these messages are subcommand headers
    _setup_simulate_subparser(
        subparser,
        f"{HEADER}\ndisto simulate: demultiplex pooled reads to sample files by index/barcode",
    )
    _setup_mapcall_subparser(
        subparser,
        f"{HEADER}\ndisto mapcall: map reads, call variants, find recombinants",
    )
    _setup_filter_subparser(
        subparser,
        f"{HEADER}\ndisto filter: filter mapcall VCF to biallelic SNPs + heterozygous sites only",
    )
    _setup_infer_subparser(
        subparser,
        f"{HEADER}\ndisto infer: infer crossover map from recombinants",
    )
    _setup_collapse_subparser(
        subparser,
        f"{HEADER}\ndisto collapse: collapse read-level COs into locus-level events",
    )
    _setup_plot_subparser(
        subparser,
        f"{HEADER}\ndisto plot: plot crossover distributions",
    )
    return parser


def main() -> None:
    try:
        command_line()
    except KeyboardInterrupt:
        logger.error("interrupted by user. Shutting down.")
        sys.exit(1)
    except Exception as exc:
        logger.exception("Unexpected error: see traceback below.")
        raise exc


def command_line() -> None:
    parser = setup_parsers()
    args = parser.parse_args()

    # If no subcommand was provided
    if args.subcommand is None:
        parser.print_help()
        sys.exit(0)

    allowed = {"simulate", "mapcall", "filter", "infer", "collapse", "plot"}
    if args.subcommand not in allowed:
        parser.print_help()
        sys.exit(0)

    run_subcommand(args)


def run_subcommand(args: argparse.Namespace) -> None:
    if args.subcommand == "simulate":
        logger.info("------------------------------------------------------------")
        logger.info("----- disto simulate: simulate recombinant gamete reads ----")
        logger.info("------------------------------------------------------------")
        run_simulate(
            reference=args.reference,
            outdir=args.out,
            prefix=args.prefix,
            heterozygosity=args.heterozygosity,
            crossover_rate=args.crossover_rate,
            nreads=args.nreads,
            read_length=args.read_length,
            random_seed=args.random_seed,
            chromosomes=args.chromosomes,
            interference=args.interference,
        )
        sys.exit(0)

    if args.subcommand == "mapcall":
        logger.info("----------------------------------------------------------------")
        logger.info("----- disto mapcall: map long reads, call and phase variants ---")
        logger.info("----------------------------------------------------------------")
        run_mapcall(
            reference=args.reference,
            reads=args.gametes,
            outdir=args.out,
            prefix=args.prefix,
            max_depth=None if args.max_depth == 0 else args.max_depth,
            min_map_q=args.min_map_q,
            min_base_q=args.min_base_q,
            threads=args.threads,
            whatshap_threads_per_worker=1,
        )
        sys.exit(0)

    if args.subcommand == "filter":
        logger.info("----------------------------------------------------------------")
        logger.info("----- disto filter: biallelic SNPs + heterozygous sites only ---")
        logger.info("----------------------------------------------------------------")
        run_filter(
            vcf_gz=args.vcf,
            outdir=args.out,
            prefix=args.prefix,
            threads=args.threads,
            bam=args.bam,
            keep_tmp=args.keep_tmp,
        )
        sys.exit(0)

    if args.subcommand == "infer":
        logger.info("----------------------------------------------------------------")
        logger.info("----- disto infer: infer crossover map from recombinants -------")
        logger.info("----------------------------------------------------------------")
        run_infer(
            reference=args.reference,
            vcf_gz=args.vcf,
            bam_path=args.bam,
            outdir=args.out,
            prefix=args.prefix,
            min_snps=args.min_snps,
            threads=args.threads,
            edge_mask_bp=args.edge_mask_bp,
        )
        sys.exit(0)

    if args.subcommand == "collapse":
        logger.info("----------------------------------------------------------------")
        logger.info("----- disto collapse: collapse read-level crossovers into loci -----")
        logger.info("----------------------------------------------------------------")
        run_collapse(
            tsv=args.infer_tsv,   # must match cli_collapse.py dest
            out_tsv=args.out,
            merge_gap=args.merge_gap,
            keep_na=args.keep_na,
        )
        sys.exit(0)

    if args.subcommand == "plot":
        logger.info("----------------------------------------------------------------")
        logger.info("----- disto plot: plot crossover distribution -----------------")
        logger.info("----------------------------------------------------------------")
        run_plot(
            tsv=args.tsv,
            outdir=args.out,
            prefix=args.prefix,
            chroms=args.chrom,
            bins=args.bins,
        )
        sys.exit(0)

    # Should never reach here
    raise SystemExit(f"Unknown subcommand: {args.subcommand}")


if __name__ == "__main__":
    main()

