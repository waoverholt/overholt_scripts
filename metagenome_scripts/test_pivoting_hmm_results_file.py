#!/usr/bin/env python3

"""
I am trying to summarize the output from search_bins_for_hmms.py to better merge with a checkM/gtdbtk summary spreadsheet.

Ideally will give:
bin    hmm_mod1    hmm_mod2 ...
bin1      3            1    ....
"""

import pandas as pd
import sys
import os

df = pd.read_csv("/data/Projects/AquaDiva/Probst_MG/MiSeq_Metagenomes/transfers_from_lumos/metawrap_bins/hmm_model_search_results.txt",
                 sep="\t", index_col=None)
wide_df = df.pivot_table(index='Bin', columns='hmm_mod', aggfunc='size', fill_value=0)
print(wide_df)