#!/usr/bin/env python3
"""
Written by Will A. Overholt
17 March 2020

Usage:
filter_ko.py <modified_kofamscan_results_file.txt>

Info:
This script was written to bridge kofamscan and KEGGDecoder.
It essentially keeps only the good hits from kofamscan.
'good hits' are defined as either:
(1) the score is > threshold*0.8
(2) score is > 100, when there is no threshold (thrshld == '-')

Generating the input file:
(1) Concatenate all the gene calls from all the bins into 1 document:
for file in $(find ./ -name "*faa"); do bin=$(basename $file); bin=${bin%_gene_seqs.faa}; sed 's/^>/>'"$bin"'_/' $file >> ../all_bin_gene_calls.faa; done

(2) Run kofamscan on the renamed AA sequences
exec_annotation -o kofamscan/all_bins_detailed_output_prok.txt -p ~/databases/kofamscan/db/profiles/prokaryote.hal --cpu 10 -f detail all_bin_gene_calls.faa

(3) Delete the stupid * in the output file
sed -i 's/*/ /' all_bins_detailed_output_prok_fixed.txt
"""

import sys,os
import pandas as pd
import itertools

#in_file = "/home/li49pol/data/Projects/Probst_MG/Nextseq_MG/10_Probst_Bins_Anvio/overholt_concatenated_summary_file/kofamscan/all_bins_detailed_output_prok_fixed.txt"
in_file = sys.argv[1]

ko_hits = {}
with open(in_file, 'r') as f:
    for line in itertools.islice(f,3, None):
        #elems = line.split()
        #print(elems)
        gene_name, KO, thrshld, score, evalue = line.split()[0:5]
        try:
            ko_hits[gene_name].append([KO, thrshld, score, evalue])
        except KeyError:
            ko_hits[gene_name] = [[KO, thrshld, score, evalue]]

# I actually think this is redundant, the list already seems to be sorted from kofamscam, oops
best_hits = {}
for gene in ko_hits:
    best_hit = []
    high_score = -100 #some of the hits were -10....

    for hit in ko_hits[gene]:
        if float(hit[2]) > high_score:
            best_hit = hit
            high_score = float(hit[2])
    best_hits[gene] = best_hit

# Filter the hits
for hit in best_hits:
    try:
        KO, thrshld, score, evalue = best_hits[hit]
        #some KOs don't have a threshold
        if thrshld == '-':
            #the median thrshld score from kofamscan is 298.4
            #the mean is 441
            #Q1 is 140
            if float(score) > 100:
                print(hit, KO)
        else:
            #found the threshold overly conservative, using 80% of the value
            if float(score) > float(thrshld)*0.8:
                print(hit, KO)
    except ValueError:
        sys.exit("The gene: {0} did not have a best_hit stored correctly.\nThe value stored is: {1}. It should be [KO, thrshld, score, evalue].\nExiting now".format(hit, best_hits[hit]))
