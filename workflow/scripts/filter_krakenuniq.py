#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sys import argv
import pandas as pd

krakenuniq_output = argv[1]

n_unique_kmers = int(argv[2])
n_tax_reads = int(argv[3])
organism_list = argv[4]

# Read and filter KrakenUniq output
kraken_output_df = pd.read_csv(krakenuniq_output, delimiter="\t", comment="#")
organism_list_df = pd.read_csv(organism_list, delimiter="\t", header=None)
print(organism_list_df)
print(f"Original data set dimensions: {kraken_output_df.shape}")

kraken_output_df = kraken_output_df[kraken_output_df["kmers"] > n_unique_kmers]

print(
    f"Data set dimensions after breadth of coverage filter: {kraken_output_df.shape}"
)
kraken_output_df = kraken_output_df[kraken_output_df["taxReads"] > n_tax_reads]
kraken_output_df = kraken_output_df[kraken_output_df["rank"] == "species"]
kraken_output_df = kraken_output_df.sort_values(by="%", ascending=False)

print(
    f"Data set dimensions after depth of coverage filter: {kraken_output_df.shape}"
)
kraken_output_df.to_csv(f"{krakenuniq_output}.filtered", sep="\t", index=False)

kraken_output_df = kraken_output_df[
    kraken_output_df.taxID.isin(organism_list_df.iloc[:, 0])
]
kraken_output_df.to_csv(
    f"{krakenuniq_output}.pathogens", sep="\t", index=False
)
