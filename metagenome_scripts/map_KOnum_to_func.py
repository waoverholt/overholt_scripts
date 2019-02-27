#!/usr/bin/env python3

import sys, os
import pandas as pd

"""
The annotations come from using the hmm collection that Calle's friend created. 

hmmsearch --cpu 8 -E 1e-5 --tblout kegg.txt ~/Projects/GATech/Databases/CalleFriendKO_HMM/genes.hmm.gz all-chambers-geneseqs.faa
grep -v "^#" kegg.txt > kegg_annotations.txt
awk -v OFS='\t' '{print $1, "KEGG", $3, $5}' kegg_annotations.txt > kegg_annotations_anvio.txt
sed -i -r 's/(K[0-9]+)_.*\t/\1\t/' kegg_annotations_anvio.txt

"""
kegg_annotations_from_hmms = sys.argv[1]
kegg_database_map = sys.argv[2]
OUT = sys.argv[3]

kegg_annotations = pd.read_csv(kegg_annotations_from_hmms, header=None, sep="\t", names = ["gene_callers_id", "source", "accession", "e_value"])
kegg_function = pd.read_csv(kegg_database_map, header=None, sep="\t", names = ["accession", "function", "Category1", "Category2", "Category3"])

new_df = kegg_annotations.set_index('accession').join(kegg_function.set_index('accession'))

new_df["cat_funct"] = new_df['function'].astype(str) + ";" + new_df['Category1'].astype(str) + ";" + new_df["Category2"].astype(str) + ";" + new_df["Category3"].astype(str)
new_df.reset_index(drop=False, inplace = True)
new_df = new_df.iloc[:, [1,2,0,8,3] ]
new_df.columns.values[3] = 'function'
new_df.to_csv(OUT, sep="\t", index=False)
