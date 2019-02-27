#!/usr/bin/env python3
"""
Written by Will Overholt
26 Feb 2019

Purpose: Run a collection of single hmm models against a directory of predicted protein sequences from genome bins
& generate a consolidated table of hits.

The hmms I was working with when I wrote this script were the collection prepared by Karthik Anantharaman
and available here: https://github.com/banfieldlab/metabolic-hmms

The collection of hmm models can have an optional tab delimited file that gives a user specified cut-off value.
Otherwise a default e-value < 1e-5 is used.

This script expects the checkM style output, where each bin contains a file named "genes.faa"
It can easily be modified to work on any faa file, or a different extension if they are still
amino acid sequences.

usage:
python search_bins_for_hmms.py -i <input directory> -m <hmm model or directory of hmm models>
                               -c <optional: tsv file of bitscore cutoffs for each hmm> -o <optional: output txt file>
                               --sum_outfile
"""

####################
#IMPORTS
import sys
import os
import subprocess
import tempfile
import pandas as pd
import argparse
####################
#Argument Parser
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_dir', required=True,
                    help='directory of bins containing a "genes.faa" file, typically generated from checkM')
parser.add_argument('-m', '--hmm_mods', required=True,
                    help='an hmm model or a directory containing hmm models')
parser.add_argument('-c', '--cut_off_values', required=False,
                     help='a tab delimited file of hmm_model_name\tcut-off score to use')
parser.add_argument('-o', '--output_file', required=False,
                     help='a tab delimited file of the results, sorted by Bin name and hmm_model; if omitted unsorted'
                          'values will be printed to the terminal (can be sorted by piping to sort -k1,1 -k2,2')
parser.add_argument('--sum_outfile', help="a tab delimited file that summarizes the results. Generates a wide format"
                    "data table of Bin x hmm_mods, with the number of hits for each model as the cell value. I use this"
                    "to add information to a summary sheet for each bin (eg. checkM values, gtdbtk tax, interesting"
                    "features. This is just a boolean flag, it will use your outfile basename", required=False,
                    action='store_true')
args = parser.parse_args()
#####################
#Functions
def list_hmm_mods(arg):
    """
    If user gives are directory of hmm models, return list for each model
    :param arg:
    :return:
    """
    if os.path.isdir(arg):
        l_mods = [hmm for hmm in os.listdir(arg) if hmm.endswith(".hmm")]
        l_mods = [os.path.join(os.path.abspath(arg), hmm) for hmm in l_mods]
        return l_mods
    else:
        return [os.path.abspath(arg)]

def parse_cutoff_file(in_file):
    mods = {}
    with open(in_file) as f:
        for line in f:
            (mod, cutoff) = line.split()
            mods[mod] = cutoff
        return mods

def match_mod_with_cutoff(mod_path, cutoff_dict):
    mod = os.path.basename(mod_path)
    return cutoff_dict[mod]

def hmmsearch_call(filename, model, temp_outfile, cutoff=None, force_evalue=False):
    #Does not check to see if hmmsearch worked without errors!
    if force_evalue:
        #print("using evalue")
        cmd = ['hmmsearch', '--noali', '--tblout', temp_outfile, '-E', '1e-5', model, filename]
    elif cutoff:
        #print("using user cutoff of: {}".format(cutoff))
        cmd=['hmmsearch', '--noali', '--tblout', temp_outfile, '-T', cutoff, model, filename]
    else:
        #print("using default evalue 1e-5")
        cmd = ['hmmsearch', '--noali', '--tblout', temp_outfile, '-E', '1e-5', model, filename]
    ret = subprocess.Popen(cmd, stdin=subprocess.PIPE, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = ret.communicate()
    return out, err, ret.returncode

def find_hmm_hits(tmp_file):
    with open(tmp_file) as f:
        lines = [line.strip() for line in f if not line.startswith("#")]
        return lines

######################
#Main Script

if __name__ == '__main__':

    df = pd.DataFrame(columns=["Bin", "GeneSeqID", "hmm_mod", "Evalue", "bitScore"])
    hmm_mods_to_run = list_hmm_mods(args.hmm_mods)
    force_evalue=False
    if args.cut_off_values:
        cutoff_dict = parse_cutoff_file(args.cut_off_values)
    else:
        force_evalue=True

    for hmm_mod in hmm_mods_to_run:
        if force_evalue:
            "user did not provide cutoff values so we'll just use the default evalue"
        else:
            mod_cutoff = match_mod_with_cutoff(hmm_mod, cutoff_dict)
            for dir in os.listdir(args.input_dir):
                gene_seqs = os.path.join(os.path.abspath(args.input_dir), dir, "genes.faa")
                tempout = tempfile.NamedTemporaryFile(suffix="_hmmresults.txt", dir="/tmp")
                hmmsearch_result = hmmsearch_call(gene_seqs, hmm_mod, tempout.name, mod_cutoff, force_evalue)
                if hmmsearch_result[2] == 0:
                    "If hmmer returns with errorcode=0 [no error]"
                    hmm_hits = find_hmm_hits(tempout.name)
                    if hmm_hits:
                        for item in hmm_hits:
                            items = item.split()
                            if args.output_file:
                                df = df.append({'Bin': dir, 'GeneSeqID': items[0], 'hmm_mod': items[2],
                                                'Evalue': items[4], 'bitScore': items[5]}, ignore_index=True)
                            else:
                                print('{0}\t{1}\t{2}\t{3}\t{4}'.format(dir, items[0], items[2], items[4], items[5]))

                    tempout.close()
                else:
                    "If there is an error, print it and exit the code"
                    print(hmmsearch_result[1])
                    exit()

    if args.output_file:
        df = df.sort_values(by=['Bin', 'hmm_mod'])
        df.to_csv(args.output_file, sep="\t", index=False)
        if args.sum_outfile:
            wide_df = df.pivot_table(index='Bin', columns='hmm_mod', aggfunc='size', fill_value=0)
            sum_fh = os.path.join(os.path.dirname(args.output_file),
                      os.path.splitext(os.path.basename(args.output_file))[0] + "_summary.txt")
            wide_df.to_csv(sum_fh, sep="\t", index=False)