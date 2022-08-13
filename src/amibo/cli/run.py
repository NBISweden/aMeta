"""
Run ancient microbiome snakemake workflow

Information on how to setup workflow here

Begin by initializing and editing a configuration file:

  amibo config --init

Once configuration settings are to your liking, test the workflow

  amibo run -j 1
"""
import logging
import subprocess as sp

from . import SNAKEFILE


logger = logging.getLogger(__name__)


def run(args):
    cmd = f"snakemake -s {SNAKEFILE} {' '.join(args.extra_options)}"
    logger.info(f"running command '{cmd}'")
    try:
        sp.run(cmd, check=True, shell=True)
    except sp.CalledProcessError:
        logger.error(f"{cmd} failed")
        raise


def add_run_subcommand(subparsers):
    parser = subparsers.add_parser(
        "run",
        help=__doc__.split("\n", maxsplit=2)[1],
        description=__doc__,
    )
    parser.add_argument(
        "--configfile",
        action="store",
        default="config/config.yaml",
        help="configuration file",
    )
    parser.add_argument(
        "--test", action="store_true", default=False, help="run small test workflow"
    )
    parser.set_defaults(runner=run)
