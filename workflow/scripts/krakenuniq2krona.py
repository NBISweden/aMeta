#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sys import argv
import pandas as pd

krakenuniq_output = argv[1]
krakenuniq_seqs = argv[2]

# Read KrakenUniq output, previously filtered
kraken_output_df = pd.read_csv(krakenuniq_output, delimiter="\t", comment="#")
kraken_output_df.taxID.to_csv(
    f"{krakenuniq_output}_taxIDs_kmers1000.txt", index=False, header=False
)
selected_taxids = kraken_output_df.taxID.values

count = 0
with open(krakenuniq_seqs, "r") as seqs, open(
    f"{krakenuniq_seqs}_kmers1000.txt", "w"
) as filtered_seqs:
    for line in seqs:
        if int(line.split("\t")[2]) in selected_taxids:
            filtered_seqs.write(line)
            count += 1

print(
    f"Sequence data set dimensions after selecting reads corresponding to filtered KrakenUniq output: {count}"
)
