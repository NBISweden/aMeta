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

.. note::

   Since the workflow is still work in progress, make sure to first
   run the commands with the `--dry-run` (`-n`) flag to make sure you
   don't inadvertedly have to regenerate large parts of the results.
   Many workflow dependencies are complex and in particular when
   running smaller parts of the workflow, unexpected things may
   happen.

Workflow graph
==============


.. _ancient-microbiome-smk: https://github.com/NBISweden/ancient-microbiome-smk
.. _{{ snakemake.config["__workflow_commit__"] }}: {{ snakemake.config["__workflow_commit_link__"] }}
