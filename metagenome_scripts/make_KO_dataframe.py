#!/usr/bin/env python3
"""
Written by Will Overholt
13 June 2018

Usage:
make_KO_dataframe.py input_directory output_file.tab 

Previous steps:
(1) reads to uniprot using diamondblast with tabular output
(2) convert to counts files using enveomics BlastTab.seqdepth.pl
(3) convert uniprot names to KO with Uniprot_KEGG_link.py

"""

import os, sys
import pandas as pd
indir = sys.argv[1]

narrow_results = pd.DataFrame(columns=["KO_Name", "Counts", "Sample"])
for file in os.listdir(indir):
    sample_id = file.split(".")[0]
    df = pd.read_csv(os.path.join(indir,file), header=None, sep=" ", names = ["KO_Name", "Counts"])
    df["Sample"] = sample_id
    print("Processing: {0}\tDimensions: {1}".format(file, df.shape))
    narrow_results = narrow_results.append(df)
    #print(GO_df.memory_usage(index=True).sum())

wide_results = narrow_results.pivot("KO_Name", "Sample", "Counts")
wide_results = wide_results.fillna(0)
print("\nFinalTable: {0}\tDimensions: {1}".format(sys.argv[2], wide_results.shape)) 
wide_results.to_csv(sys.argv[2], sep="\t")
