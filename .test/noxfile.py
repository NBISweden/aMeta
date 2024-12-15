"""Nox sessions."""

import os
import re
import shlex
import shutil
import sys
from pathlib import Path
from textwrap import dedent

import nox


TESTDIR = Path(os.path.dirname(os.path.abspath(__file__)))
ROOTDIR = TESTDIR.parent
SNAKEFILE = "../workflow/Snakefile"


try:
    from nox import session
except ImportError:
    message = f"""\
    Nox failed to import the 'nox' package.

    Please install it using the following command:

    {sys.executable} -m pip install nox"""
    raise SystemExit(dedent(message)) from None

python_snakemake_versions = [
    ("3.12", "8.25.5"),
    ("3.11", "7.32.4"),
    ("3.8", "6.3.0"),
]
nox.needs_version = ">= 2021.6.6"


def pip_install_pulp(session, snakemake_version):
    """Pip install correct pulp version.

    cf https://github.com/snakemake/snakemake/issues/2607.
    """
    major, minor, patch = map(int, snakemake_version.split("."))
    if (major >= 8) and (minor >= 1) and (patch >= 2):
        session.install("pulp>=2.8")
    else:
        session.install("pulp<2.8")


def parse_sem(version):
    """Parse semantic version."""
    return tuple(map(int, version.split(".")))


def get_conda_frontend(session):
    """Determine the conda frontend."""
    major, minor, patch = map(int, session.run("conda", "--version", silent=True).split()[1].split("."))
    if (major >= 24):
        return "conda"
    if (major >= 23) and (minor >= 10) and (patch >= 0):
        return "conda"
    return "mamba"


def conda_env_dir(snakemake_version, python_version):
    """Get conda environment directory."""
    pver = "-".join([str(x) for x in parse_sem(python_version)[0:2]])
    sver = "-".join([str(x) for x in parse_sem(snakemake_version)])
    envdir = TESTDIR / "envs" / f".nox/snakemake-python-{pver}-snakemake-{sver}"
    return envdir


def snakemake_install_deps(session, conda_env_dir, conda_frontend):
    """Install snakemake dependencies."""
    cmd = [
        "snakemake",
        "-s", SNAKEFILE,
        "--use-conda", "--conda-create-envs-only",
        "--conda-frontend", conda_frontend,
        "--show-failed-logs",
        "--conda-prefix", str(conda_env_dir),
        "--cores", "1",
    ]
    session.run(*cmd)


def setup_krona(session, conda_env_dir, conda_frontend):
    """Setup krona."""
    envfile = ROOTDIR / "workflow" / "envs" / "krona.yaml"
    result = session.run(
        "snakemake", "-s", SNAKEFILE,
        "--use-conda", "--conda-frontend", conda_frontend,
        "--list-conda-envs", "--conda-prefix", str(conda_env_dir),
        "--cores", "1", silent=True
    )
    kronaenv = [x.split("\t") for x in result.split("\n") if re.search("krona.yaml", x)]
    if len(kronaenv) == 0:
        raise ValueError("No krona environment found")
    env = ROOTDIR / ".test" / kronaenv[0][2]
    updateTaxonomy = env / "opt" / "krona" / "updateTaxonomy.sh"
    taxonomy = env / "opt" / "krona" / "taxonomy"
    session.run("bash", str(updateTaxonomy), str(taxonomy), external=True)


def adjust_malt_memory_usage(session, conda_env_dir, conda_frontend):
    """Adjust memory usage for malt."""
    envfile = ROOTDIR / "workflow" / "envs" / "malt.yaml"
    result = session.run(
        "snakemake", "-s", SNAKEFILE,
        "--use-conda", "--conda-frontend", conda_frontend,
        "--list-conda-envs", "--conda-prefix", str(conda_env_dir),
        "--cores", "1", silent=True
    )
    maltenv = [x.split("\t") for x in result.split("\n") if re.search("malt.yaml", x)]
    if len(maltenv) == 0:
        raise ValueError("No malt environment found")
    env = ROOTDIR / ".test" / maltenv[0][2]
    m = re.search(
        "MALT \(version ([0-9.]+),",
        str(session.run(str(env / "bin" / "malt-run"), "-h", silent=True))
    )
    malt_version = m.group(1)
    major, minor, patch = parse_sem(malt_version)
    malt_opt = env / "opt" / f"malt-{major}.{minor}{patch}"
    malt_build = malt_opt / "malt-build.vmoptions"
    txt = malt_build.read_text()
    txt = re.sub(r"-Xmx[0-9]+[gG]", "-Xmx3G", txt)
    malt_build.write_text(txt)
    malt_run = malt_opt / "malt-run.vmoptions"
    txt = malt_run.read_text()
    txt = re.sub(r"-Xmx[0-9]+[gG]", "-Xmx3G", txt)
    malt_run.write_text(txt)


@session(venv_backend="conda")
@nox.parametrize(
    "python,snakemake",
    [
        (python, snakemake)
        for python, snakemake in python_snakemake_versions
    ],
)
def snakemake(session, snakemake):
    """Run snakemake."""
    conda_frontend = get_conda_frontend(session)

    session.conda_install("pytest")
    session.conda_install(
        f"snakemake={snakemake}",
        channel=["conda-forge", "bioconda"],
    )
    if parse_sem(snakemake) >= (8, 0, 0):
        session.conda_install(
            "snakemake-storage-plugin-http",
            channel=["conda-forge", "bioconda"],
        )
    envdir = conda_env_dir(snakemake, session.python)
    snakemake_install_deps(session, envdir, conda_frontend)
    setup_krona(session, envdir, conda_frontend)
    adjust_malt_memory_usage(session, envdir, conda_frontend)

    session.run("pytest", "tests.py")
