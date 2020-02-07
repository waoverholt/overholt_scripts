#!/usr/bin/env python3
"""
Written by Will Overholt
20 June 2019

Purpose: Run a subset (or the entire KEGG KO collection) of single hmm models against a directory of predicted protein sequences from genome bins
& generate a consolidated table of hits.

It makes the most sense to generate your own .hal file that contains each of the KO.hmm models you are interested in (e.g. all annamox genes).

This is modified from my "search_bins_for_hmms.py" program.

This script expects the checkM style output, where each bin contains a file named "genes.faa".
Each file within the bins does need to have the same name, but this can be modified if needed.

It can easily be modified to work on any faa file with the -f flag, if they are still
amino acid sequences.

usage:
python search_bins_for_hmms.py -i <input directory> -p <list_hmms.hal> -t <input_file_name>
                               -o <optional: output txt file>
                               --sum_outfile
"""

####################
#IMPORTS
import sys
import os
import subprocess
import pandas as pd
import argparse
import shutil

####################
#Argument Parser
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_dir', required=True,
                    help='directory of bins containing a "genes.faa" file, typically generated from checkM. This can be modified with -t')
parser.add_argument('-p', '--profile', required=True,
                    help='a kofamscan style .hal profile')
parser.add_argument('-o', '--output_file', required=False,
                     help='a tab delimited file of the results, sorted by Bin name and KO; if omitted unsorted'
                          'values will be printed to the terminal (can be sorted by piping to sort -k1,1 -k2,2')
parser.add_argument('--sum_outfile', help="a tab delimited file that summarizes the results. Generates a wide format"
                    "data table of Bin x hmm_mods, with the number of hits for each model as the cell value. I use this"
                    "to add information to a summary sheet for each bin (eg. checkM values, gtdbtk tax, interesting"
                    "features. This is just a boolean flag, it will use your outfile basename", required=False,
                    action='store_true')
parser.add_argument('-f', '--file_name', required=False, nargs='?', default="genes.faa", type=str,
                    help='The name of the protein sequence file found within each bin')
parser.add_argument('-e', '--evalue', required=False, nargs='?', default="1e-6", type=float, 
                    help='Minimum evalue for a hit to be counted')
args = parser.parse_args()
#####################
#Functions
def kofamscan_call(filename, profile):
    cmd = ['exec_annotation', '-p', profile, '--tmp-dir=./kofamscan_tmp', filename]
    ret = subprocess.Popen(cmd, stdin=subprocess.PIPE, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = ret.communicate()
    if err:
        print(err)
        quit()
    return out.decode("utf-8"), err, ret.returncode

######################
#Main Script
##Some genes are missing a threshold score, meaning I miss true hits (see anammox genes...)
##So I'm going to parse hits that have an evalue < 1e-6, which is very likely too low!! (??)
if __name__ == '__main__':
    df = pd.DataFrame(columns=["Bin", "GeneSeqID", "KO", "thrshld", "score", "evalue", "KO Definition"])
    for dir in os.listdir(args.input_dir):
        gene_seqs = os.path.join(os.path.abspath(args.input_dir), dir, args.file_name)
        result = kofamscan_call(gene_seqs, args.profile)[0].split('\n')
        hits = [hit.strip() for hit in result if not hit.startswith('#')]
        for hit in hits:
            items = hit.split()
            if len(items) > 1:
                if args.output_file:
                    #if you want all hits, you need to remove the column that kofamscan prints with the "*" [hit above thrshold]
                    #here I assume you only want true hits, and I skip the astrick (items[0])
                    if items[0] == "*":
                        df = df.append({'Bin': dir, 'GeneSeqID': items[1], 'KO': items[2],
                                    'thrshld': items[3], 'score': items[4], 'evalue': items[5],
                                    'KO Definition': " ".join(items[6:])}, ignore_index=True)
                    else:
                        if float(items[4]) < args.evalue:
                            df = df.append({'Bin': dir, 'GeneSeqID': items[0], 'KO': items[1],
                                    'thrshld': items[2], 'score': items[3], 'evalue': items[4],
                                    'KO Definition': " ".join(items[5:])}, ignore_index=True)
                else:
                    if items[0] == "*":
                        print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}'.format(dir, items[1], items[2], items[3], items[4], items[5], " ".join(items[6:])))
                    else:
                        if float(items[4]) < args.evalue:
                            print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}'.format(dir, items[0], items[1], items[2], items[3], items[4], " ".join(items[5:])))
    if args.output_file:
        df = df.sort_values(by=['Bin', 'KO'])
        df.to_csv(args.output_file, sep="\t", index=False)
        if args.sum_outfile:
            wide_df = df.pivot_table(index='Bin', columns='KO', aggfunc='size', fill_value=0)
            sum_fh = os.path.join(os.path.dirname(args.output_file),
                        os.path.splitext(os.path.basename(args.output_file))[0] + "_summary.txt")
            wide_df.to_csv(sum_fh, sep="\t")

##Clean up temp directory produced by kofamscan
kofamscan_temp_dir = os.path.join(os.getcwd(), "kofamscan_tmp")
shutil.rmtree(kofamscan_temp_dir)