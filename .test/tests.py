import os
import re
import pytest
import shutil
import subprocess as sp
from pathlib import Path


TESTDIR = Path(os.path.dirname(os.path.abspath(__file__)))
ROOTDIR = TESTDIR.parent
SNAKEFILE = "../workflow/Snakefile"


def get_conda_frontend():
    """Determine the conda frontend."""
    major, minor, patch = map(
        int, sp.run(["conda", "--version"], capture_output=True, text=True)
        .stdout.split()[1].split("."))
    if (major >= 24):
        return "conda"
    if (major == 23) and (minor >= 10) and (patch >= 0):
        return "conda"
    return "mamba"


def parse_sem(version):
    """Parse semantic version."""
    return tuple(map(int, version.split(".")))


testfiles = [
    "config/config.yaml",
    "config/envmodules.yaml",
    "data/bar.fq.gz",
    "data/barfoo.fq.gz",
    "data/foo.fq.gz",
    "data/foobar.fq.gz",
    "resources/KrakenUniq_DB/library/ref.fa",
    "resources/KrakenUniq_DB/seqid2taxid.map",
    "resources/KrakenUniq_DB/taxonomy/names.dmp",
    "resources/KrakenUniq_DB/taxonomy/nodes.dmp",
    "resources/accession2taxid.map",
    "resources/pathogenomesFound.tab",
    "resources/ref.fa",
    "resources/samples.tsv",
    "resources/seqid2taxid.pathogen.map",
]


@pytest.fixture
def conda_frontend():
    """Return conda frontend."""
    return get_conda_frontend()


@pytest.fixture
def snakemake_version():
    """Return Snakemake version."""
    version = sp.run(["snakemake", "--version"], capture_output=True, text=True)
    return version.stdout.strip()


@pytest.fixture
def python_version():
    """Return Python version."""
    version = sp.run(["python", "--version"], capture_output=True, text=True)
    return version.stdout.strip().split()[1]


@pytest.fixture
def conda_env_dir(snakemake_version, python_version):
    """Return conda environment directory."""
    pver = "-".join([str(x) for x in parse_sem(python_version)[0:2]])
    sver = "-".join([str(x) for x in parse_sem(snakemake_version)])
    envdir = TESTDIR / "envs" / f".nox/snakemake-python-{pver}-snakemake-{sver}"
    return envdir


@pytest.fixture(name="deps")
def install_deps(conda_env_dir, conda_frontend, num_cores):
    """Install conda environment dependencies."""
    options = [
        "-s", SNAKEFILE,
        "--conda-create-envs-only", "--conda-frontend",
        conda_frontend, "--use-conda", "--show-failed-logs",
        "--conda-prefix", str(conda_env_dir),
        "--cores", str(num_cores)
    ]
    cmd = ["snakemake"] + options

    sp.run(cmd, check=True, text=True)


@pytest.fixture
def testdir(tmpdir_factory):
    """Copy test files."""
    dname = tmpdir_factory.mktemp("test")
    for testfile in testfiles:
        outfile = dname.join(testfile)
        outfile.dirpath().ensure_dir()
        shutil.copy(testfile, outfile)
    return dname


@pytest.fixture(name="krakenuniq")
def setup_krakenuniq_db(conda_env_dir, deps, conda_frontend, testdir):
    """Set up KrakenUniq database."""
    envfile = ROOTDIR / "workflow" / "envs" / "krakenuniq.yaml"
    result = sp.run(["snakemake", "--use-conda", "-s", SNAKEFILE,
            "--conda-frontend", conda_frontend,
            "--list-conda-envs",
            "--conda-prefix", str(conda_env_dir), "-j", "1"], check=True,
                    text=True, capture_output=True)
    krakenenv = [x.split("\t") for x in result.stdout.strip().split("\n") if re.search("krakenuniq.yaml", x)]
    if len(krakenenv) == 0:
        raise ValueError("No krakenuniq environment found")
    env = ROOTDIR / ".test" / krakenenv[0][2]
    krakenuniq_build = env / "bin" / "krakenuniq-build"
    jellyfish_bin = env / "bin" / "jellyfish"
    krakenuniq_db = testdir / "resources/KrakenUniq_DB"
    result = sp.run([str(krakenuniq_build), "--db", str(krakenuniq_db),
            "--kmer-len", "21", "--minimizer-len", "11",
            "--jellyfish-bin", str(jellyfish_bin)], check=True,
           text=True)


@pytest.fixture(name="krona")
def setup_krona():
    """Set up Krona."""
    kronadir = ROOTDIR / "workflow" / "envs" / "krona"
    sp.run(["./updateTaxonomy.sh"], cwd=kronadir, check=True, text=True)


def test_workflow(testdir, conda_env_dir, conda_frontend, num_cores, krakenuniq):
    """Test workflow."""
    options = [
        "--conda-frontend", conda_frontend,
        "--conda-prefix", str(conda_env_dir),
        "--use-conda", "--show-failed-logs", "--cores", str(num_cores),
        "-s", "../workflow/Snakefile",
        "--directory", str(testdir)
    ]
    cmd = ["snakemake"] + options

    sp.run(cmd, check=True, text=True)
