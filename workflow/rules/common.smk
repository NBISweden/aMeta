import os
import sys
import subprocess as sp
from pathlib import Path
from snakemake.utils import validate, logger
import pandas as pd
import contextlib
from config import WORKFLOW_DIR

# context manager for cd
@contextlib.contextmanager
def cd(path, logger):
    CWD = os.getcwd()
    logger.info("Changing directory from {} to {}".format(CWD, path))

    os.chdir(path)
    try:
        yield
    except Exception as e:
        logger.warning(e)
        logger.warning("Exception caught: ".format(sys.exc_info()[0]))
    finally:
        logger.info("Changing directory back to {}".format(CWD))
        os.chdir(CWD)


##### load config and sample sheets #####
configfile: "config/config.yaml"


validate(config, schema="../schemas/config.schema.yaml")

kw = {"sep": "\t" if config["samplesheet"].endswith(".tsv") else ","}
samples = pd.read_csv(config["samplesheet"], **kw).set_index("sample", drop=False)

samples.index.names = ["sample_id"]
validate(samples, schema="../schemas/samples.schema.yaml")


##############################
## Store some workflow metadata
##############################
config["__workflow_basedir__"] = workflow.basedir
config["__workflow_workdir__"] = os.getcwd()
config["__worfklow_commit__"] = None
config["__workflow_commit_link__"] = None

try:
    with cd(workflow.basedir, logger):
        commit = sp.check_output(["git", "rev-parse", "HEAD"]).decode().strip()
        commit_short = (
            sp.check_output(["git", "rev-parse", "--short", "HEAD"]).decode().strip()
        )
        config["__workflow_commit__"] = commit_short
        config[
            "__workflow_commit_link__"
        ] = f"https://github.com/NBISweden/manticore-smk/commit/{commit}"
except Exception as e:
    print(e)
    raise


##############################
# Global variables
##############################
#
# Store some config values in all-caps global vars
#
SAMPLES = samples["sample"].tolist()


##############################
# Wildcard constraints
##############################
#
# Restrict some globally used wildcards for enhanced performance
#
wildcard_constraints:
    sample=f"({'|'.join(samples['sample'].tolist())})",


##############################
# Input collection functions
##############################
def all_input(wildcards):
    d = {
        "multiqc.before": rules.MultiQC_BeforeTrimming.output,
        "multiqc.after": rules.MultiQC_AfterTrimming.output,
        "mapdamage": mapdamage_input(wildcards),
        "krakenuniq.krona": expand(
            "results/KRAKENUNIQ/{sample}/taxonomy.krona.html", sample=SAMPLES
        ),
        "malt.abundance": "results/MALT_ABUNDANCE_MATRIX/malt_abundance_matrix.txt",
        "auth": authentication_input(wildcards),
    }
    return d


def mapdamage_input(wildcards):
    if not config["analyses"]["mapdamage"]:
        return []
    return expand("results/MAPDAMAGE/{sample}", sample=SAMPLES)


def authentication_input(wildcards):
    if not config["analyses"]["authentication"]:
        return []
    return expand("results/AUTHENTICATION/{sample}", sample=SAMPLES)
