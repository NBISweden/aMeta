![Logo](aMeta.png)

# aMeta: an accurate and memory-efficient ancient Metagenomic profiling workflow

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.10.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Tests](https://github.com/NBISweden/ancient-microbiome-smk/actions/workflows/main.yaml/badge.svg)](https://github.com/NBISweden/ancient-microbiome-smk/actions/workflows/main.yaml)

## About

aMeta is a Snakemake workflow for identifying microbial sequences in ancient DNA shotgun metagenomics samples. The workflow performs:

- trimming adapter sequences and removing reads shorter than 30 bp with Cutadapt
- quaity control before and after trimming with FastQC and MultiQC
- taxonomic sequence kmer-based classification with KrakenUniq
- sequence alignment with Bowtie2 and screening for common microbial pathogens
- deamination pattern analysis with MapDamage2
- Lowest Common Ancestor (LCA) sequence alignment with Malt
- authentication and validation of identified microbial species with MaltExtract

When using aMeta and / or pre-built databases provided together with the wokflow for your research projects, please cite our preprint: https://www.biorxiv.org/content/10.1101/2022.10.03.510579v1

## Authors

* Nikolay Oskolkov (@LeandroRitter) nikolay.oskolkov@scilifelab.se
* Claudio Mirabello (@clami66) claudio.mirabello@scilifelab.se
* Per Unneberg (@percyfal) per.unneberg@scilifelab.se

## Installation

Clone the repository, then create and activate aMeta conda environment:

    git clone https://github.com/NBISweden/aMeta
    cd aMeta
    conda env create -f workflow/envs/environment.yaml
    # alternatively: mamba env create -f workflow/envs/environment.yaml
    conda activate aMeta

Run a test to make sure that the workflow was installed correctly:

    cd .test
    ./runtest.sh -j 20

Here, and below, by `-j` you can specify the number of threads that the workflow can use.

## Quick start

To run the worflow you need to prepare a sample-file `config/samples.tsv` and a configuration file `config/config.yaml`, below we provide examples for both files. 

Here is an example of `samples.tsv`, this implies that the fastq-files files are located in `aMeta/data` folder:

    sample	fastq
    foo	data/foo.fq.gz
    bar	data/bar.fq.gz

Below is an example of `config.yaml`, here you will need to download a few databases that we made public (or build databases yourself).

    samplesheet: "config/samples.tsv"

    # KrakenUniq Microbial NCBI NT database (if you are interested in prokaryotes only)
    # can be downloaded from https://doi.org/10.17044/scilifelab.20518251
    krakenuniq_db: resources/DBDIR_KrakenUniq_MicrobialNT

    # KrakenUniq full NCBI NT database (if you are interested in prokaryotes and eukaryotes)
    # can be downloaded from https://doi.org/10.17044/scilifelab.20205504
    #krakenuniq_db: resources/DBDIR_KrakenUniq_Full_NT

    # Bowtie2 index and helping files for following up microbial pathogens 
    # can be downloaded from https://doi.org/10.17044/scilifelab.21185887
    bowtie2_patho_db: resources/library.pathogen.fna
    pathogenomesFound: resources/pathogensFound.very_inclusive.tab
    pathogenome_seqid2taxid_db: resources/seqid2taxid.pathogen.map

    # Bowtie2 index for full NCBI NT (for quick followup of prokaryotes and eukaryotes)
    # can be downloaded from https://doi.org/10.17044/scilifelab.21070063
    #bowtie2_patho_db: resources/library.fna

    # Helping files for building Malt database 
    # can be downloaded from https://doi.org/10.17044/scilifelab.21070063
    malt_nt_fasta: resources/library.fna
    malt_seqid2taxid_db: resources/seqid2taxid.map.orig
    malt_accession2taxid: resources/nucl_gb.accession2taxid

    # A path for downloading NCBI taxonomy files (performed automatically)
    # you do not need to change this line
    ncbi_db: resources/ncbi

    # Breadth and depth of coverage filters 
    # default thresholds are very conservative, can be tuned by users
    n_unique_kmers: 1000
    n_tax_reads: 200


After you have prepared the sample- and configration-file, the workflow can be run using the following command line:

    cd aMeta
    snakemake --snakefile workflow/Snakefile --use-conda -j 20


In the next sections we will give more information about the parameters in the configuration file as well as instructions on how to run the workflow in a computer cluster enviroment.


## More configuration options

Within `config.yaml` one can specify what samples to analyse in the `samples` section through the `include` and `exclude` keys, so that a global samplesheet can be reused multiple times.

Analyses `mapdamage`, `authentication`, `malt`, and `krona` can be individually turned on and off in the `analyses` section.

Adapter sequences can be defined in the `adapters` section of `config.yaml`. 
The keys `config['adapters']['illumina']` (default `true`) and `config['adapters']['nextera']` (default `false`) are switches
that turn on/off adapter trimming of illumina (`AGATCGGAAGAG`) and nextera (`AGATCGGAAGAG`) adapter sequences. Addional custom adapter sequences can be set in the configuration key
`config['adapters']['custom']` which must be an array of strings.

An example snippet that can optionally be added to the configuration file `config.yaml` is shown below:

    # you can include or exclude samples
    samples:
      include:
        - foo
      exclude:
        - bar

    # you can include or exclude certain types of analysis
    analyses:
      mapdamage: true
      authentication: true
      malt: true
      krona: true

    # you can specify type of adapters to trim
    adapters:
      illumina: true
      nextera: false
      # custom is a list of adapter sequences
      custom: []


### Environment module configuration

To run the workflow in a computer cluster environemnt you should specify environmental modules and runtimes via `--profile` as follows:

    snakemake --snakefile workflow/Snakefile -j 100 --profile .profile --use-envmodules

If the workflow is run on a HPC with the `--use-envmodules` option
(see
[using-environment-modules](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#using-environment-modules)),
the workflow will check for an additional file that configures environment modules. By default, the file is `config/envmodules.yaml`, but a custom location can be set with the
environment variable `ANCIENT_MICROBIOME_ENVMODULES`.

Environmental modules configurations are placed in a configuration section
`envmodules` with key-value pairs that map a dependency set to a list
of environment modules. The dependency sets are named after the rule's
corresponding conda environment file, such that a dependency set may
affect multiple rules. For instance, the following example shows how
to define modules for rules depending on fastqc, as it would be
implemented on the [uppmax](https://uppmax.uu.se/) compute cluster:

    envmodules:
      fastqc:
        - bioinfo-tools
        - FastQC

See the configuration schema file
(`workflows/schema/config.schema.yaml`) for more information.

### Runtime configuration

Most individual rules define the number of threads to run. Although
the number of threads for a given rule can be tweaked on the command
line via the option `--set-threads`, it is advisable to put all
runtime configurations in a
[profile](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html).
At its simplest, a profile is a directory (e.g. `.profile` in the
working directory) containing a file `config.yaml` which consists of
command line option settings. In addition to customizing threads, it
enables the customization of resources, such as runtime and memory. An
example is shown here:

    # Rerun incomplete jobs
    rerun-incomplete: true
    # Restart jobs once on failure
    restart-times: 1
    # Set threads for mapping and fastqc
    set-threads:
      - Bowtie2_Pathogenome_Alignment=10
      - FastQC_BeforeTrimming=5
    # Set resources (runtime in minutes, memory in mb) for malt
    set-resources:
      - Malt:runtime=7200
      - Malt:mem_mb=512000
    # Set defalt resources that apply to all rules
    default-resources:
      - runtime=120
      - mem_mb=16000
      - disk_mb=1000000

For more advanced profiles for different hpc systems, see [Snakemake-Profiles github page](https://github.com/snakemake-profiles).

