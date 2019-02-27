#!/usr/bin/env python3
"""
Written by Will Overholt
12 June 2018

This script converts a uniprot counts file to a KO counts file. The counts file can be created using the enveomics script BlastTab.seqdepth.pl on a standard tab-delimited blastout from reads against uniprot.

The script needs the mapping file that was created using the script "~/Projects/GATech/scripts/kegg_gene_2_KO_lookupAPI.py". This is EXTRMELY SLOW (4 days) due to API restrictions in biopython. The results file can be found at ~/Projects/GATech/Databases/kegg_gene_to_KO_map.txt

It is run:
Uniprot_KEGG_link.py counts_file.txt kegg_gene_to_KO_map.txt > output

NOTE!!
This will only grab those KO that map to your blast file. It will give different results for each counts file that will be resolved in the next step (making the table)
"""

import sys, os, re
import collections

count = sys.argv[1]
mapping = sys.argv[2]

kegg_dict = {}
with open(mapping, "r") as f:
    for line in f:
        line = line.rstrip()
        elems = line.split("\t")
        kegg_dict[elems[0]] = elems[2]
kegg_dict = collections.OrderedDict(sorted(kegg_dict.items()))

"""
for key, value in kegg_dict.items():
    print(key, value)
"""
results = {} 
with open(count, "r") as f:
    for line in f:
        line = line.rstrip()
        elems = line.split("\t")
        UID = elems[0]
        UID = UID.split("|")
        UID = UID[1]
        if UID in kegg_dict:
            if kegg_dict[UID] in results:
                results[kegg_dict[UID]] = float(results[kegg_dict[UID]]) + float(elems[3])
            else:
                results[kegg_dict[UID]] = float(elems[3])
       
results = collections.OrderedDict(sorted(results.items()))

for key, value in results.items():
    print(key, value)

