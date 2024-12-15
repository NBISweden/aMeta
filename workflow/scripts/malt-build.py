#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Per Unneberg"
__copyright__ = "Copyright 2022, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

import re
from snakemake.shell import shell
import subprocess as sp

out = sp.run(["malt-build", "--help"], capture_output=True)
regex = re.compile(r"version (?P<major>\d+)\.(?P<minor>\d+)")
m = regex.search(out.stderr.decode())
if m is None:
    # Assume minor version 4
    minor = 4
else:
    minor = int(m.groupdict()["minor"])

a2t_option = "-a2taxonomy" if minor <= 4 else "-a2t"

log = snakemake.log_fmt_shell(stdout=False, stderr=True, append=True)
options = snakemake.params.get("extra", "")
unique_taxids = snakemake.input.unique_taxids

seqid2taxid = snakemake.params.seqid2taxid
nt_fasta = snakemake.params.nt_fasta
accession2taxid = snakemake.params.accession2taxid

output = snakemake

shell(
    "grep -wFf {unique_taxids} {seqid2taxid} > {snakemake.output.seqid2taxid_project}; "
    "cut -f1 {snakemake.output.seqid2taxid_project} > {snakemake.output.seqids_project}; "
    "grep -Ff {snakemake.output.seqids_project} {nt_fasta} | sed 's/>//g' > {snakemake.output.project_headers}; "
    "seqtk subseq {nt_fasta} {snakemake.output.project_headers} > {snakemake.output.project_fasta} {log}; "
    "unset DISPLAY; "
    "malt-build -i {snakemake.output.project_fasta} {a2t_option} {accession2taxid} -s DNA -t {snakemake.threads} -d {snakemake.output.db} {log}"
)
