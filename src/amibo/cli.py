"""Console script for datasources"""
import logging
import os
import subprocess as sp
import sys
from argparse import ArgumentParser

from . import __version__

__author__ = "Per Unneberg"


logger = logging.getLogger(__name__)

ROOT = os.path.abspath(os.path.dirname(__file__))

SNAKEFILE = os.path.join(ROOT, "workflow", "Snakefile")
TESTDIR = os.path.join(ROOT, ".test")

__RUNDOC__ = [
    "Run ancient microbiome snakemake workflow"
    " "
    "Information on how to setup workflow here"
    " "
    "Begin by initializing and editing a configuration file: "
    "  amibo config --init"
    " "
    "Once configuration settings are to your liking, test the workflow"
    "  amibo run -j 1"
]


def run(args):
    if args.test:
        args.extra_options += ["--directory", TESTDIR]
    cmd = f"snakemake -s {SNAKEFILE} {' '.join(args.extra_options)}"
    logger.info(f"running command '{cmd}'")
    try:
        sp.run(cmd, check=True, shell=True)
    except sp.CalledProcessError:
        logger.error(f"{cmd} failed")
        raise


def main(arg_list=None):
    if arg_list is None:
        arg_list = sys.argv[1:]
    logging.basicConfig(
        level=logging.INFO, format="%(levelname)s [%(name)s:%(funcName)s]: %(message)s"
    )
    top_parser = ArgumentParser(description=__doc__, prog="amibo")
    top_parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )
    top_parser.add_argument(
        "--debug", action="store_true", default=False, help="Print debug messages"
    )

    subparsers = top_parser.add_subparsers(dest="subcommand")
    subparsers.required = True

    parser = subparsers.add_parser(
        "run",
        help="Run ancient microbiome snakemake workflow",
        description="\n".join(__RUNDOC__),
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

    args, extra = top_parser.parse_known_args(arg_list)

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    del args.debug

    args.extra_options = extra

    args.runner(args)

    return 0
