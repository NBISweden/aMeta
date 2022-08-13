"""
Run small test workflow

The workflow includes a small test data set. By default, run tests in
the test directory. Note that the test automatically adds the
--use-conda flag to generate isolated software stack environments.

"""
import copy
import json
import logging
import os
import shutil
import subprocess as sp

import pkg_resources
from amibo.utils import cd
from snakemake.deployment.conda import Env
from snakemake.workflow import Workflow

from . import TESTDIR
from .run import run

logger = logging.getLogger(__name__)


def get_environment_path(envname, dname):
    # Deduce yaml file using snakemake conda
    snakefile = pkg_resources.resource_filename("amibo", "workflow/Snakefile")
    yaml = pkg_resources.resource_filename("amibo", f"workflow/envs/{envname}.yaml")
    workflow = Workflow(
        snakefile=snakefile, use_conda=True, rerun_triggers=["mtime", "input"]
    )
    env_dir = os.path.join(os.path.realpath(dname), ".snakemake/conda")
    env = Env(workflow, env_file=yaml, env_dir=env_dir)
    return os.path.join(env_dir, env.hash)


def build_krakenuniq_database(dname):
    path = get_environment_path("krakenuniq", dname)
    cmd = (
        f"{path}/bin/krakenuniq-build --db resources/KrakenUniq_DB "
        f"--kmer-len 21 --minimizer-len 11 --jellyfish-bin {path}/bin/jellyfish"
    )
    with cd(dname):
        sp.run(cmd, check=True, shell=True)


def build_krona_taxonomy(dname):
    path = get_environment_path("krona", dname)
    cmd = "./updateTaxonomy.sh taxonomy"
    with cd(os.path.join(path, "opt", "krona")):
        sp.run(cmd, check=True, shell=True)


def adjust_malt_max_memory_usage(dname):
    path = get_environment_path("malt", dname)
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
        "--use-conda",
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
    if args.test_dir is None:
        args.test_dir = TESTDIR
    if not args.no_init:
        if os.path.exists(args.test_dir):
            logger.debug(f"test directory {args.test_dir} exists: skip testdir setup")
        else:
            copytree_testdir(args.test_dir)
        snakemake_init_conda_envs(copy.deepcopy(args))

        # Setup databases
        build_krakenuniq_database(args.test_dir)
        build_krona_taxonomy(args.test_dir)
        adjust_malt_max_memory_usage(args.test_dir)

    args.extra_options += ["--directory", args.test_dir, "--use-conda"]
    run(args)


def add_test_subcommand(subparsers):
    parser = subparsers.add_parser(
        "test",
        help=__doc__.split("\n", maxsplit=2)[1],
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
    parser.set_defaults(runner=test)
