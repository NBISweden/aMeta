#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gzip
from sys import argv
import re

sam = argv[1]
ids = argv[2]

counts = {}
read_ids = set()
taxids = []

with open(ids) as id_file:
    for id in id_file:
        counts[id.strip()] = 0
        taxids.append(id.strip())

exp = re.compile("[|]tax[|]([0-9]+)[|]*")

with gzip.open(sam, "rt") as sam_file:
    # M_ST-E00266:253:HCK7GCCXY:1:1101:20577:16463       0       AK125553.1|tax|63363    ...
    for read in sam_file:
        if "@" not in read:
            read_id = read.split()[0]
            taxid = re.findall(exp, read)[0]
            if (
                taxid in counts and (taxid, read_id) not in read_ids
            ):  # want to avoid counting duplicate reads
                counts[taxid] += 1
                read_ids.add((taxid, read_id))

for taxid in taxids:
    print(counts[taxid])
