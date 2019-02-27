#!/usr/bin/env python3
"""
Written by Will Overholt
10 December 2018

Usage:
script.py <directory_of_checkm_results_bins> <outdir>
"""

import os, sys

indir = sys.argv[1]
outdir = sys.argv[2]

for dir in os.listdir(indir):
    bin_id = dir.split(".")[0]
    cur_fp = os.path.join(os.path.abspath(indir), bin_id, "genes.faa")
    outdir_bin_path = os.path.join(os.path.abspath(outdir), bin_id)
    if not os.path.exists(outdir_bin_path):
        os.makedirs(outdir_bin_path)
    outfile_p = os.path.join(outdir_bin_path, "genes.faa")
    outfile = open(outfile_p, "w")
    with open(cur_fp) as seq_fp:
        for line in seq_fp:
            if line[0] == ">":
                new_name = line[1:]
                new_name = "".join((">", bin_id, " # ", new_name))
                outfile.write(new_name)
            else:
                outfile.write(line)
    outfile.close()

   

