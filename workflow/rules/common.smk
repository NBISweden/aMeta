import os
import sys
import subprocess as sp
from pathlib import Path
from snakemake.utils import validate, logger
import pandas as pd
import contextlib
from config import WORKFLOW_DIR
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

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


if workflow.use_env_modules:
    envmodules = os.getenv("ANCIENT_MICROBIOME_ENVMODULES", "config/envmodules.yaml")

    configfile: envmodules


validate(config, schema="../schemas/config.schema.yaml")


kw = {"sep": "\t" if config["samplesheet"].endswith(".tsv") else ","}
samples = pd.read_csv(config["samplesheet"], **kw).set_index("sample", drop=False)
samples.index.names = ["sample_id"]
if "exclude" in config["samples"]:
    logger.info(
        "Excluding samples {exclude} from analysis".format(
            exclude=",".join(f"'{x}'" for x in config["samples"]["exclude"])
        )
    )
    samples = samples[~samples.index.isin(config["samples"]["exclude"])]
if "include" in config["samples"]:
    logger.info(
        "Restricting analysis to samples {incl}".format(
            incl=",".join(f"'{x}'" for x in config["samples"]["include"])
        )
    )
    samples = samples[samples.index.isin(config["samples"]["include"])]

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
        "multiqc.after": rules.MultiQC.output,
        "mapdamage": mapdamage_input(wildcards),
        "krakenuniq.krona": krona_input(wildcards),
        "malt.abundance": malt_input(wildcards),
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
    return expand("results/AUTHENTICATION/.{sample}_done", sample=SAMPLES)


def malt_input(wildcards):
    if not config["analyses"]["malt"]:
        return []
    return (
        "results/MALT_ABUNDANCE_MATRIX_SAM/malt_abundance_matrix_sam.txt",
        "results/MALT_ABUNDANCE_MATRIX_RMA6/malt_abundance_matrix_rma6.txt",
    )


def krona_input(wildcards):
    if not config["analyses"]["krona"]:
        return []
    return expand("results/KRAKENUNIQ/{sample}/taxonomy.krona.html", sample=SAMPLES)


def multiqc_input(wildcards):
    """Collect all inputs to multiqc"""
    d = {
        "fastqc_before_trimming": expand(
            "results/FASTQC_BEFORE_TRIMMING/{sample}_fastqc.zip", sample=SAMPLES
        ),
        "fastqc_after_trimming": expand(
            "results/FASTQC_AFTER_TRIMMING/{sample}.trimmed_fastqc.zip",
            sample=SAMPLES,
        ),
        "cutadapt": expand(
            "logs/CUTADAPT_ADAPTER_TRIMMING/{sample}.log", sample=SAMPLES
        ),
        "bowtie2": expand("logs/BOWTIE2/{sample}.log", sample=SAMPLES),
    }
    return d


def aggregate_PMD(wildcards):
    checkpoint_output = checkpoints.Extract_TaxIDs.get(sample=wildcards.sample).output[
        0
    ]
    return expand(
        "results/AUTHENTICATION/{sample}/{taxid}/PMD_plot.frag.pdf",
        sample=wildcards.sample,
        taxid=glob_wildcards(os.path.join(checkpoint_output, "{taxid,[0-9]+}")).taxid,
    )


def aggregate_plots(wildcards):
    checkpoint_output = checkpoints.Extract_TaxIDs.get(sample=wildcards.sample).output[
        0
    ]
    return expand(
        "results/AUTHENTICATION/{sample}/{taxid}/authentic_Sample_{sample}.trimmed.rma6_TaxID_{taxid}.pdf",
        sample=wildcards.sample,
        taxid=glob_wildcards(os.path.join(checkpoint_output, "{taxid,[0-9]+}")).taxid,
    )


def aggregate_post(wildcards):
    checkpoint_output = checkpoints.Extract_TaxIDs.get(sample=wildcards.sample).output[
        0
    ]
    return expand(
        "results/AUTHENTICATION/{sample}/{taxid}/{sample}.trimmed.rma6_MaltExtract_output/analysis.RData",
        sample=wildcards.sample,
        taxid=glob_wildcards(os.path.join(checkpoint_output, "{taxid,[0-9]+}")).taxid,
    )


def get_ref_id(wildcards):
    ref_id = {wildcards.taxid}
    with open(
        f"results/AUTHENTICATION/{wildcards.sample}/{wildcards.taxid}/{wildcards.sample}.trimmed.rma6_MaltExtract_output/default/readDist/{wildcards.sample}.trimmed.rma6_additionalNodeEntries.txt"
    ) as f:
        contents = f.readlines()
        try:
            ref_id = contents[-1].split(";")[1][1:]
        except:
            pass
    return ref_id
