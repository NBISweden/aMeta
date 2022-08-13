ancient-microbiome-smk results
==============================

ancient-microbiome-smk_ is a Snakemake workflow for detection of
microbiomes in ancient DNA samples. This analysis is based on
commit version {{ snakemake.config["__workflow_commit__"] }}_.


The analysis can be rerun with the following command:

.. code-block:: shell

   cd {{ snakemake.config["__workflow_workdir__"] }}
{% if snakemake.config["__workflow_basedir__"] == snakemake.config["__workflow_workdir__"] %}
   snakemake -j 1
{% else %}
   snakemake -j 1 -s {{ snakemake.config["__workflow_basedir__"] }}/Snakefile
{% endif %}

and the report

.. code-block:: shell

   cd {{ snakemake.config["__workflow_workdir__"] }}
{% if snakemake.config["__workflow_basedir__"] == snakemake.config["__workflow_workdir__"] %}
   snakemake --report report.html
{% else %}
   snakemake -j 1 -s {{ snakemake.config["__workflow_basedir__"] }}/Snakefile --report report.html
{% endif %}

.. note::

   Since the workflow is still work in progress, make sure to first
   run the commands with the `--dry-run` (`-n`) flag to make sure you
   don't inadvertedly have to regenerate large parts of the results.
   Many workflow dependencies are complex and in particular when
   running smaller parts of the workflow, unexpected things may
   happen.

Workflow summary
================

- adapter trimming of sequences with cutadapt_
- fastqc_ before and after trimming
- taxonomic sequence classification with krakenuniq_
- sequence alignment with malt_
- sequence damage analysis with mapdamage2_
- authentication of identified sequences


MultiQC report
=================

The `MultiQC report`_ collects results from fastqc_, bowtie2_, and
cutadapt_.


Workflow graph
==============


.. _ancient-microbiome-smk: https://github.com/NBISweden/ancient-microbiome-smk
.. _{{ snakemake.config["__workflow_commit__"] }}: {{ snakemake.config["__workflow_commit_link__"] }}
.. _MultiQC report: ./results/MULTIQC/multiqc_report.html
.. _fastqc: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _bowtie2: https://github.com/BenLangmead/bowtie2
.. _cutadapt: https://cutadapt.readthedocs.io/en/stable/
.. _krakenuniq: https://github.com/fbreitwieser/krakenuniq
.. _malt: https://github.com/husonlab/malt
.. _mapdamage2: https://ginolhac.github.io/mapDamage/
