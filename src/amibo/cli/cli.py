"""Console script for amibo"""
import logging
import sys
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

from amibo.cli.run import add_run_subcommand
from amibo.cli.test import add_test_subcommand
from amibo.templates import env

from .. import __version__

__author__ = "Per Unneberg"


logger = logging.getLogger(__name__)


class HelpfulArgumentParser(ArgumentParser):
    """An ArgumentParser that prints full help on errors."""

    def __init__(self, *args, **kwargs):
        if "formatter_class" not in kwargs:
            kwargs["formatter_class"] = RawDescriptionHelpFormatter
        super().__init__(*args, **kwargs)

    def error(self, message):
        self.print_help(sys.stderr)
        args = {"prog": self.prog, "message": message}
        self.exit(2, "%(prog)s: error: %(message)s\n" % args)


def sample_config(args):
    template = env.get_template("amibo.yaml.j2")
    print(template.render())


def add_sample_config_subcommand(subparsers):
    parser = subparsers.add_parser(
        "sample_config",
        help="Produce a sample amibo.yaml configuration file",
        description=None,
    )
    parser.set_defaults(runner=sample_config)


def main(arg_list=None):
    if arg_list is None:
        arg_list = sys.argv[1:]
    logging.basicConfig(
        level=logging.INFO, format="%(levelname)s [%(name)s:%(funcName)s]: %(message)s"
    )
    top_parser = HelpfulArgumentParser(description=__doc__, prog="amibo")
    top_parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )
    top_parser.add_argument(
        "--debug", action="store_true", default=False, help="Print debug messages"
    )

    subparsers = top_parser.add_subparsers(dest="subcommand")
    subparsers.required = True
    add_run_subcommand(subparsers)
    add_sample_config_subcommand(subparsers)
    add_test_subcommand(subparsers)

    args, extra = top_parser.parse_known_args(arg_list)

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    del args.debug

    args.extra_options = extra

    args.runner(args)

    return 0
