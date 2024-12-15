#!/bin/bash
set -xeuo pipefail
echo "generating conda-linux-64.lock from environment.yml"
conda-lock --kind explicit --platform linux-64 -f workflow/envs/environmentv8.yaml
mv conda-linux-64.lock conda-linux-64.v8.lock
