#!/usr/bin/env python3
"""
Written by Will Overholt
25 June 2018

This script merges samples with identical names that were sequenced twice
"""

import sys, os, re
import shutil

#Inputs
RUN1 = sys.argv[1]
RUN2 = sys.argv[2]
OUTDIR = sys.argv[3]


for file1 in sorted(os.listdir(RUN1)):
    #run through sorted directory
    file1_name = file1.split("_")[0]
    for file2 in sorted(os.listdir(RUN2)):
        file2_name = file2.split("_")[0]
        #only merge if file names are the same
        if file1_name == file2_name:
            #byte based copy to OUTPUT file
            OUT = open(os.path.join(OUTDIR, file1), "wb")
            with open(os.path.join(os.path.abspath(RUN1), file1), 'rb') as f1:
                shutil.copyfileobj(f1, OUT, 1024*1024*10)
            with open(os.path.join(os.path.abspath(RUN2), file2), 'rb') as f2:
                shutil.copyfileobj(f2, OUT, 1024*1024*10)
            OUT.close()
            #print to terminal the files that are merged (can pipe to file)
            print(os.path.join(os.path.abspath(RUN1), file1), os.path.join(os.path.abspath(RUN2), file2), sep = "\t")
