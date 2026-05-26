#!/bin/bash
set -xeuo pipefail

SNAKEMAKE_VERSION=$(snakemake --version)
SNAKEMAKE_MAJOR_VERSION=$(echo $SNAKEMAKE_VERSION | cut -d. -f1)

if [ $SNAKEMAKE_MAJOR_VERSION -lt 8 ]; then
    ENVIRONMENTFILE=environment.yaml
else
    ENVIRONMENTFILE=environment.v8.yaml
fi

echo "generating conda-linux-64.lock from $ENVIRONMENTFILE"
conda-lock --kind explicit --platform linux-64 -f workflow/envs/$ENVIRONMENTFILE

if [ $SNAKEMAKE_MAJOR_VERSION -lt 8 ]; then
    mv conda-linux-64.lock conda-linux-64.v7.lock
else
    mv conda-linux-64.lock conda-linux-64.v8.lock
fi
