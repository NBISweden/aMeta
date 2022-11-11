![Logo](aMeta.png =400x400)

# Snakemake workflow: aMeta

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.10.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Tests](https://github.com/NBISweden/ancient-microbiome-smk/actions/workflows/main.yaml/badge.svg)](https://github.com/NBISweden/ancient-microbiome-smk/actions/workflows/main.yaml)

## About

Snakemake workflow for identifying microbe sequences in ancient DNA
samples. The workflow does:

- adapter trimming of sequences
- FastQC before and after trimming
- taxonomic sequence classification with KrakenUniq
- sequence alignment with Malt
- sequence damage analysis with Mapdamage2
- authentication of identified sequences

## Quickstart

Clone the repo, create and edit the configuration files (see below)
and run

    cd /path/to/workdir
    snakemake -s /path/to/repo/workflow/Snakefile -j 100 --profile .profile --use-envmodules

## Authors

* Nikolay Oskolkov (@LeandroRitter)
* Claudio Mirabello (@clami66)
* Per Unneberg (@percyfal)

## Configuration

The workflow requires a configuration file, by default residing in
`config/config.yaml` relative to the working directory, that defines
location of samplesheet, what samples and analyses to run, and
location of databases. The configuration file is validated against a
schema (`workflow/schemas/config.schema.yaml`) that can be consulted
for more detailed information regarding configuration properties.

The `samplesheet` key points to a samplesheet file that consists of at
minimum two columns, sample and fastq:

    sample	fastq
    foo     data/foo.fq.gz
    bar     /path/to/data/bar.fq.gz

What samples to analyse can be constrained in the `samples` section
through the `include` and `exclude` keys, so that a global samplesheet
can be reused multiple times.

Analyses `mapdamage`, `authentication`, `malt`, and `krona` can be
individually turned on and off in the `analyses` section.

Adapter sequence can be defined in the `adapters` configuration
section. The keys `config['adapters']['illumina']` (default `true`)
and `config['adapters']['nextera']` (default `false`) are switches
that turn on/off adapter trimming of illumina (`AGATCGGAAGAG`) and
nextera (`AGATCGGAAGAG`) adapter sequences. Addional custom adapter
sequences can be set in the configuration key
`config['adapters']['custom']` which must be an array of strings.

Database locations are defined by the following keys:

`krakenuniq_db`: path to KrakenUniq database

`bowtie2_patho_db`: Full path to Bowtie2 pathogenome database

`pathogenome_path`: Path to Bowtie2 pathogenome database, excluding
the database name

`pathogenomesFound`: List of pathogens to keep when filtering
KrakenUniq output

`malt_seqid2taxid_db`: Sequence id to taxonomy mapping

`malt_nt_fasta`: Fasta library

`malt_accession2taxid`: Accession to taxonomy id mapping

A minimal configuration example is shown below:

    samplesheet: resources/samples.tsv
    samples:
      include:
        - foo
        - bar
      exclude:
        - foobar

    analyses:
      mapdamage: false
      authentication: false
      malt: false
	  
	adapters:
	  illumina: true
	  nextera: false
	  # custom is a list of adapter sequences
	  custom: []

    # Databases
    krakenuniq_db: resources/KrakenUniq_DB
    bowtie2_patho_db: resources/ref.fa
    pathogenome_path: resources
    pathogenomesFound: resources/pathogenomesFound.tab
    malt_seqid2taxid_db: resources/KrakenUniq_DB/seqid2taxid.map
    malt_nt_fasta: resources/ref.fa
    malt_accession2taxid: resources/accession2taxid.map

### Environment module configuration

If the workflow is run on a HPC with the `--use-envmodules` option
(see
[using-environment-modules](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#using-environment-modules)),
the workflow will check for an additional configuration file that
configures environment modules. By default, the file is
`config/envmodules.yaml`, but a custom location can be set with the
environment variable `ANCIENT_MICROBIOME_ENVMODULES`.

envmodules configurations are placed in a configuration section
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

For more advanced profiles for different hpc systems, see
[Snakemake-Profiles github
page](https://github.com/snakemake-profiles).

## Usage

If you use this workflow in a paper, don't forget to give credits to
the authors by citing the URL of this (original) repository and, if
available, its DOI (see above).

### Step 1: Obtain a copy of this workflow

1. Create a new github repository using this workflow [as a
   template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
2. [Clone](https://help.github.com/en/articles/cloning-a-repository)
   the newly created repository to your local system, into the place
   where you want to perform the data analysis.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files
in the `config/` folder. Adjust `config.yaml` to configure the
workflow execution, and `samples.tsv` to specify your sample setup.

### Step 3: Install Snakemake

Install Snakemake using
[conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda create -c bioconda -c conda-forge -n snakemake snakemake

For installation details, see the [instructions in the Snakemake
documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --cluster qsub --jobs 100

or

    snakemake --use-conda --drmaa --jobs 100

If you not only want to fix the software stack but also the underlying OS, use

    snakemake --use-conda --use-singularity

in combination with any of the modes above. See the [Snakemake
documentation](https://snakemake.readthedocs.io/en/stable/executable.html)
for further details.

### Step 5: Investigate results

After successful execution, you can create a self-contained
interactive HTML report with all results via:

    snakemake --report report.html

This report can, e.g., be forwarded to your collaborators. An example
(using some trivial test data) can be seen
[here](https://cdn.rawgit.com/snakemake-workflows/rna-seq-kallisto-sleuth/master/.test/report.html).

### Step 6: Commit changes

Whenever you change something, don't forget to commit the changes back
to your github copy of the repository:

    git commit -a
    git push

### Step 7: Obtain updates from upstream

Whenever you want to synchronize your workflow copy with new
developments from upstream, do the following.

1. Once, register the upstream repository in your local copy: `git
   remote add -f upstream
   git@github.com:snakemake-workflows/ancient-microbiome-smk.git` or
   `git remote add -f upstream
   https://github.com/snakemake-workflows/ancient-microbiome-smk.git`
   if you do not have setup ssh keys.
2. Update the upstream version: `git fetch upstream`.
3. Create a diff with the current version: `git diff HEAD
   upstream/master workflow > upstream-changes.diff`.
4. Investigate the changes: `vim upstream-changes.diff`.
5. Apply the modified diff via: `git apply upstream-changes.diff`.
6. Carefully check whether you need to update the config files: `git
   diff HEAD upstream/master config`. If so, do it manually, and only
   where necessary, since you would otherwise likely overwrite your
   settings and samples.


### Step 8: Contribute back

In case you have also changed or added steps, please consider
contributing them back to the original repository:

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the
   original repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository)
   the fork to your local system, to a different place than where you
   ran your analysis.
3. Copy the modified files from your analysis to the clone of your
   fork, e.g., `cp -r workflow path/to/fork`. Make sure to **not**
   accidentally copy config file contents or sample sheets. Instead,
   manually update the example config files if necessary.
4. Commit and push your changes to your fork.
5. Create a [pull
   request](https://help.github.com/en/articles/creating-a-pull-request)
   against the original repository.

## Testing

Test cases are in the subfolder `.test`. They are automatically
executed via continuous integration with [Github
Actions](https://github.com/features/actions).
