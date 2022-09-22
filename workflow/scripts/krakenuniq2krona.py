#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sys import argv
import pandas as pd

krakenuniq_output = argv[1]
krakenuniq_seqs = argv[2]

n_unique_kmers = int(argv[3])
n_tax_reads = int(argv[4])

# Read and filter KrakenUniq output
kraken_output_df = pd.read_csv(krakenuniq_output, delimiter="\t", comment="#")
print(f"Original data set dimensions: {kraken_output_df.shape}")

kraken_output_df = kraken_output_df[kraken_output_df["kmers"] > n_unique_kmers]

print(f"Data set dimensions after breadth of coverage filter: {kraken_output_df.shape}")

kraken_output_df = kraken_output_df[kraken_output_df["taxReads"] > n_tax_reads]
kraken_output_df = kraken_output_df[kraken_output_df["rank"] == "species"]
kraken_output_df.rename(columns={"%": "Pers_Reads"}, inplace=True)
kraken_output_df = kraken_output_df.sort_values(by="Pers_Reads", ascending=False)

print(f"Data set dimensions after depth of coverage filter: {kraken_output_df.shape}")
kraken_output_df.to_csv(f"{krakenuniq_output}.filtered", sep="\t", index=False)
kraken_output_df.taxID.to_csv(
    f"{krakenuniq_output}_taxIDs_kmers{n_unique_kmers}.txt", index=False, header=False
)
selected_taxids = kraken_output_df.taxID.values

count = 0
with open(krakenuniq_seqs, "r") as seqs, open(
    f"{krakenuniq_seqs}_kmers{n_unique_kmers}.txt", "w"
) as filtered_seqs:
    for line in seqs:
        if int(line.split("\t")[2]) in selected_taxids:
            filtered_seqs.write(line)
            count += 1

print(
    f"Sequence data set dimensions after selecting reads corresponding to filtered KrakenUniq output: {count}"
)
