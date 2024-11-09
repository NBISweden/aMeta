import os
import sys
import subprocess as sp
from pathlib import Path
from snakemake.utils import validate, logger
import pandas as pd
import contextlib
from config import WORKFLOW_DIR
try:
    from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
except:
    pass
from snakemake.io import Wildcards

try:
    HTTP = HTTPRemoteProvider()
except:
    pass

def storage_wrapper(url):
    try:
        fn = HTTP.remote(url, keep_local=True)
    except:
        fn = storage.http(f"https://{url}")
    return fn


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


use_env_modules = False
try:
    if workflow.use_env_modules:
        use_env_modules = True
except:
    from snakemake_interface_executor_plugins.settings import DeploymentMethod
    if DeploymentMethod.ENV_MODULES in workflow.deployment_settings.deployment_method:
        use_env_modules = True

if use_env_modules:
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
        ] = f"https://github.com/NBISweden/aMeta/commit/{commit}"
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
# Set the adapter list
ADAPTERS = config["adapters"].get("custom", list())
if config["adapters"]["illumina"]:
    ADAPTERS.append("AGATCGGAAGAG")
if config["adapters"]["nextera"]:
    ADAPTERS.append("CTGTCTCTTATA")


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
        "summary": summary_input(wildcards),
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

def summary_input(wildcards):
    return "results/overview_heatmap_scores.pdf"

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


def aggregate_maltextract(wildcards):
    """Collect maltextract output directories"""
    checkpoint_output = checkpoints.Create_Sample_TaxID_Directories.get(
        sample=wildcards.sample
    ).output[0]
    taxid = glob_wildcards(
        os.path.join(os.path.dirname(checkpoint_output), "{taxid,[0-9]+}")
    ).taxid
    return expand(
        "results/AUTHENTICATION/{sample}/{taxid}/MaltExtract_output/log.txt",
        sample=wildcards.sample,
        taxid=taxid,
    )


def _aggregate_utils(fmt, wildcards):
    """Collect common output for all aggregate functions. Returns a tuple
    of lists sample, and taxid"""
    logger.debug(
        f"Running _aggregate_utils for format '{fmt}', wildcards '{dict(wildcards)}'"
    )
    res = []
    checkpoint_output = checkpoints.Create_Sample_TaxID_Directories.get(
        sample=wildcards.sample
    ).output[0]
    taxid = glob_wildcards(
        os.path.join(os.path.dirname(checkpoint_output), "{taxid,[0-9]+}")
    ).taxid
    sample = []
    refid = []
    taxid_out = []
    for tid in taxid:
        wc = Wildcards(fromdict={"sample": wildcards.sample, "taxid": tid})
        _refid = get_ref_id(wc)
        if _refid is not None and _refid != tid:
            refid.append(_refid)
            taxid_out.append(tid)
            sample.append(wildcards.sample)
    if len(refid) > 0:
        res = expand(fmt, zip, sample=sample, taxid=taxid_out)
    return res


def aggregate_PMD(wildcards):
    fmt = "results/AUTHENTICATION/{sample}/{taxid}/PMD_plot.frag.pdf"
    return _aggregate_utils(fmt, wildcards)

def aggregate_plots(wildcards):
    fmt = "results/AUTHENTICATION/{sample}/{taxid}/authentic_Sample_{sample}.trimmed.rma6_TaxID_{taxid}.pdf"
    return _aggregate_utils(fmt, wildcards)

def aggregate_scores(wildcards):
    fmt = "results/AUTHENTICATION/{sample}/{taxid}/authentication_scores.txt"
    return _aggregate_utils(fmt, wildcards)

def aggregate_post(wildcards):
    fmt = "results/AUTHENTICATION/{sample}/{taxid}/MaltExtract_output/analysis.RData"
    return _aggregate_utils(fmt, wildcards)

def get_ref_id(wildcards):
    """Return reference id for a given taxonomy id"""
    ref_id = wildcards.taxid
    infile = f"results/AUTHENTICATION/{wildcards.sample}/{wildcards.taxid}/MaltExtract_output/default/readDist/{wildcards.sample}.trimmed.rma6_additionalNodeEntries.txt"
    if not os.path.exists(infile):
        logger.debug(f"No such file {infile}; cannot extract refid")
        return None
    with open(infile) as f:
        contents = f.readlines()
        try:
            ref_id = contents[1].split(";")[1][1:]
        except:
            logger.warning(
                f"Failed to extract ref_id from {infile}; returning taxid {wildcards.taxid}"
            )
            pass
    return ref_id


def format_maltextract_output_directory(wildcards):
    """Format MaltExtract output directory name"""
    return f"results/AUTHENTICATION/{wildcards.sample}/{wildcards.taxid}/MaltExtract_output/"
