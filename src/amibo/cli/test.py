"""Run small test workflow

The workflow includes a small test data set. By default, the tests are
run in the package test directory, but this behaviour can be changed
with the --test-dir option. Prior to the first run, KrakenUniq and
krona databases are downloaded and installed. Also, the malt vmoptions
may have to be modified for systems with small memory (option
--adjust-malt-max-memory-usage).

The tests assume that all required programs and packages are
available. Use --use-conda to make the workflow install isolated
software stacks. Alternatively, if on a HPC where conda is unavailable
but there is a module system, you can apply the --use-envmodules
command. For this to work, you need to set the environment variable
ANCIENT_MICROBIOME_ENVMODULES to point to an environment modules file.

"""
import copy
import json
import logging
import os
import shutil
import subprocess as sp
import sys

import pkg_resources
import snakemake
import yaml
from amibo.utils import cd
from snakemake.deployment.conda import Env
from snakemake.deployment.env_modules import EnvModules
from snakemake.shell import shell
from snakemake.workflow import Workflow

from . import TESTDIR
from .run import run

version = pkg_resources.parse_version(snakemake.__version__)

logger = logging.getLogger(__name__)


def get_environment_path(envname, dname):
    # Deduce yaml file using snakemake conda
    snakefile = pkg_resources.resource_filename("amibo", "workflow/Snakefile")
    yaml = pkg_resources.resource_filename("amibo", f"workflow/envs/{envname}.yaml")
    if version >= pkg_resources.parse_version("7.8.0"):
        workflow = Workflow(
            snakefile=snakefile, use_conda=True, rerun_triggers=["mtime", "input"]
        )
    else:
        workflow = Workflow(snakefile=snakefile, use_conda=True)
    env_dir = os.path.join(os.path.realpath(dname), ".snakemake/conda")
    env = Env(workflow, env_file=yaml, env_dir=env_dir)
    return os.path.join(env_dir, env.hash)


def build_krakenuniq_database(args, envmodules=None):
    krakenuniq_build = "krakenuniq-build"
    jellyfish = "jellyfish"
    if args.use_conda:
        path = get_environment_path("krakenuniq", args.test_dir)
        krakenuniq_build = f"{path}/bin/krakenuniq-build"
        jellyfish = f"{path}/bin/jellyfish"
    cmd = (
        f"{krakenuniq_build} --db resources/KrakenUniq_DB "
        f"--kmer-len 21 --minimizer-len 11 --jellyfish-bin {jellyfish}"
    )
    if args.use_envmodules:
        env_modules = EnvModules(*envmodules.get("krakenuniq"))
    with cd(args.test_dir):
        shell(cmd, env_modules=env_modules)


def build_krona_taxonomy(args):
    if not args.use_conda:
        return
    updateTaxonomy = "updateTaxonomy.sh"
    path = get_environment_path("krona", args.test_dir)
    updateTaxonomy = os.path.join(path, "updateTaxonomy.sh")
    krona_home = os.path.join(path, "opt", "krona")
    cmd = f"{updateTaxonomy} taxonomy"
    with cd(krona_home):
        sp.run(cmd, check=True, shell=True)


def adjust_malt_max_memory_usage(args):
    if not args.use_conda:
        # Assume malt memory is set correctly in all but conda
        # environments
        return
    path = get_environment_path("malt", args.test_dir)
    cmd = f"conda list -p {path} malt --json"
    res = sp.check_output(cmd.split(" ")).decode()
    version = json.loads(res)[0]["version"]
    with cd(os.path.join(path, "opt", f"malt-{version}")):
        sp.run(
            'sed -i -e "s/-Xmx64GB/-Xmx3G/" malt-build.vmoptions',
            shell=True,
            check=True,
        )
        sp.run(
            'sed -i -e "s/-Xmx64GB/-Xmx3G/" malt-run.vmoptions', shell=True, check=True
        )


def snakemake_init_conda_envs(args):
    logger.debug("Initializing conda environments")
    args.extra_options = [
        "--show-failed-logs",
        "-j",
        "1",
        "--conda-cleanup-pkgs",
        "cache",
        "--conda-create-envs-only",
        "--directory",
        args.test_dir,
    ]
    run(args)


def copytree_testdir(dname):
    """Copy test directory"""
    try:
        shutil.copytree(TESTDIR, dname)
    except Exception as e:
        print(e)
        raise


def test(args):
    envmodules = dict()
    if args.info:
        logger.info("Test data files are available at:")
        logger.info(f"  {TESTDIR}")
        logger.info("TODO: more info coming")
        sys.exit(0)
    if args.test_dir is None:
        args.test_dir = TESTDIR
    if args.use_envmodules:
        envfile = (
            args.envmodules_file
            if args.envmodules_file is not None
            else os.environ.get("ANCIENT_MICROBIOME_ENVMODULES")
        )
        with open(envfile) as fh:
            try:
                envmodules = yaml.safe_load(fh)
            except AttributeError as exc:
                print(exc)
                raise
            except yaml.YAMLError as exc:
                print(exc)
                raise
    if not args.no_init:
        if os.path.exists(args.test_dir):
            logger.debug(f"test directory {args.test_dir} exists: skip testdir setup")
        else:
            copytree_testdir(args.test_dir)
        if "--use-conda" in args.extra_options:
            snakemake_init_conda_envs(copy.deepcopy(args))

            # Setup databases
        build_krakenuniq_database(args, envmodules.get("envmodules"))
        build_krona_taxonomy(args)
        adjust_malt_max_memory_usage(args)

    args.extra_options += ["--directory", args.test_dir]
    if args.use_conda:
        args.extra_options += ["--use-conda"]
    if args.use_envmodules:
        args.extra_options += ["--use-envmodules"]
        if (
            os.environ.get("ANCIENT_MICROBIOME_ENVMODULES") is None
            and args.envmodules_file is not None
        ):
            os.environ["ANCIENT_MICROBIOME_ENVMODULES"] = args.envmodules_file

    run(args)


def add_test_subcommand(subparsers):
    parser = subparsers.add_parser(
        "test",
        help=__doc__.split("\n", maxsplit=1)[0],
        description=__doc__,
    )

    parser.add_argument(
        "--test-dir",
        action="store",
        default=None,
        help="Run tests in specified directory",
    )
    parser.add_argument(
        "--no-init",
        action="store_true",
        default=False,
        help="Skip test dir init in case of rerun",
    )
    parser.add_argument(
        "--info",
        action="store_true",
        default=False,
        help="Show test data directory location and information",
    )
    parser.add_argument(
        "-a",
        "--adjust-malt-max-memory-usage",
        action="store_true",
        default=False,
        help=(
            "Set max memory usage in malt-run.vmoptions"
            " and malt-build.vmoptions to -Xmx3G"
        ),
    )
    parser.add_argument(
        "--use-conda",
        action="store_true",
        default=False,
        help=(
            "If defined in the rule, run job in a conda environment."
            " If this flag is not set, the conda directive"
            " is ignored."
        ),
    )
    parser.add_argument(
        "--use-envmodules",
        action="store_true",
        default=False,
        help=(
            "If defined in the rule, run job within the given"
            " environment modules, loaded in the given"
            " order. This can be combined with --use-conda,"
            " which will then be only"
            " used as a fallback for rules which don't"
            " define environment modules."
        ),
    )
    parser.add_argument(
        "--envmodules-file",
        action="store",
        default=None,
        help=(
            "Load environment module definitions from this file."
            " Alternatively, set the ANCIENT_MICROBIOMES_MODULES."
        ),
    )
    parser.set_defaults(runner=test)
