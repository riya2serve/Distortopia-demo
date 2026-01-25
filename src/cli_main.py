#!/usr/bin/env python


import sys
import argparse
# from .make_wide import make_wide
from .cli_simulate import _setup_simulate_subparser
from .cli_mapcall import _setup_mapcall_subparser
from .cli_infer import _setup_infer_subparser
from .cli_plot import _setup_plot_subparser
from .simulate import run_simulate
from .mapcall import run_mapcall
from .infer import run_infer
from .plot import run_plot
# from ..utils.exceptions import DistoError
# from setup_logger import set_log_level
from loguru import logger

"""
disto simulate -r REF.fa -p GAMETES -n 10_000_000 -l 100_000 -s 123
disto call -r REF.fa -g GAMETES.fastq.gz -o . -q 10 -Q 20
disto infer -r REF -v GAMETES.vcf -o . -p RATES
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
        # epilog=EPILOG,
        formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
        add_help=False,
    )
    parser.add_argument('-h', '--help', action='help', help=argparse.SUPPRESS)
    parser.add_argument("-v", "--version", action='version', version=f"ipyrad {VERSION}")
    subparser = parser.add_subparsers(help="sub-commands", dest="subcommand")

    # add subcommands: these messages are subcommand headers
    _setup_simulate_subparser(subparser, f"{HEADER}\ndisto simulate: demultiplex pooled reads to sample files by index/barcode")
    _setup_mapcall_subparser(subparser, f"{HEADER}\ndisto mapcall: map reads, call variants, find recombinants")
    _setup_infer_subparser(subparser, f"{HEADER}\ndisto infer: infer crossover map from recombinants")
    _setup_plot_subparser(subparser, f"{HEADER}\ndisto plot: plot crossover distributions")    
    return parser


def main():
    try:
        command_line()
    except KeyboardInterrupt:
        logger.error("interrupted by user. Shutting down.")
        sys.exit(1)
    # expected error, only report message no traceback
    # except DistoError as exc:
    #     logger.error(f"Error: {exc}")
    #     logger.error("see error message above. Shutting down.")
    #     sys.exit(1)
    # raise with traceback
    except Exception as exc:
        logger.exception("Unexpected error: see traceback below.")
        raise exc


def command_line():
    parser = setup_parsers()
    args = parser.parse_args()

    # LOGGING: -----------------------------------------------------
    if hasattr(args, "log_level"):
        pass
        # set_log_level(args.log_level, args.log_file)
        # logger.info("HI")

    if args.subcommand not in ["simulate", "mapcall", "infer", "plot"]:
        # NO SUBCOMMAND: print help
        parser.print_help()
        sys.exit(0)
    else:
        run_subcommand(args)


def run_subcommand(args):
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
            # log_level=args.log_level,
        )
        sys.exit(0)

    if args.subcommand == "mapcall":
        logger.info("----------------------------------------------------------------")
        logger.info("----- disto mapcall: map long reads, call and phase variants ---")
        logger.info("----------------------------------------------------------------")
        run_mapcall(
            reference=args.reference,
            gametes=args.gametes,
            outdir=args.out,
            prefix=args.prefix,
            min_map_q=args.min_map_q,
            min_base_q=args.min_base_q,
            threads=args.threads,
            # log_level=args.log_level,
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
            # log_level=args.log_level,
            threads=args.threads,
            edge_mask_bp=args.edge_mask_bp,
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
            # min_snps=args.min_snps,
            # log_level=args.log_level,
        )
        sys.exit(0)     



if __name__ == "__main__":

    main()
