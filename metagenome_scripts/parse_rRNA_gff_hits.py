#!/usr/bin/env python3
"""
Written by Will A. Overholt
2019.04.11

This script has 2 primary functions.
(1) Summarize the number of genomes with a 5S/16s/23S gene into a table
(2) Grab the sequences and add them 2 a multifasta file for each of the 3 genes.

You need bedtools in your path (conda activate metagenomics).
usage:
parse_rRNA_gff_hits.py -i rRNA_gffs_dir -f bin_fasta_dir -o output_dir

"""


import sys, os, re
import subprocess
import fnmatch
import argparse
################
#Argument Parser
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_dir', required=True,
                    help='directory of rRNA gffs produced by barrnap')
parser.add_argument('-f', '--fasta_dir', required=True,
                    help='directory of the fasta files for each bin')
parser.add_argument('-o', '--output_dir', required=False,
                     help='directory to save the 4 produced files')
args = parser.parse_args()
################

################
#Functions
def call_bedtools_fasta(bed_line, fasta, bin, gene):
    p1  = subprocess.Popen(['echo', bed_line], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(['bedtools', 'getfasta', '-bed', 'stdin', '-fi', fasta],
        stdin=p1.stdout, stdout=subprocess.PIPE)
    seq = p2.communicate()
    seq = seq[0].decode("utf-8")
    seq = seq.split("\n")[1]
    seq = ">{0}_{1}\n{2}".format(gene, bin, seq)
    return seq
##################

###Output files
if not os.path.exists(args.output_dir):
    os.mkdir(os.path.abspath(args.output_dir))

rRNA_table = open(os.path.join(os.path.abspath(args.output_dir), "rRNA_table.txt"), "w")
TSU_seqs = open(os.path.join(os.path.abspath(args.output_dir), "5S_seqs.fasta"), "w")
SSU_seqs = open(os.path.join(os.path.abspath(args.output_dir), "16S_seqs.fasta"), "w")
LSU_seqs = open(os.path.join(os.path.abspath(args.output_dir), "23S_seqs.fasta"), "w")

####

bins=[bin for bin in os.listdir(args.fasta_dir) if fnmatch.fnmatch(bin, '*.fa')]

print("binname", "5S", "16S", "23S", sep="\t", file=rRNA_table)
for filename in os.listdir(args.input_dir):
    binname = os.path.splitext(filename)[0]
    binpath = os.path.join(args.input_dir, filename)
    with open(binpath) as f:
        working_bin_fasta = [bin for bin in bins if re.search(binname+".fa", bin)]
        working_bin_fasta = os.path.join(os.path.abspath(args.fasta_dir), working_bin_fasta[0])
        SSU=0
        LSU=0
        TSU=0
        for line in f:
            if not line[0] == "#":
                elems = line.split()
                gene_name=elems[8].split("=")[2]
                if gene_name == "5S":
                    TSU=TSU+1
                    seq = call_bedtools_fasta(line, working_bin_fasta, binname, "5S")
                    print(seq, file=TSU_seqs)
                elif gene_name == "16S":
                    SSU=SSU+1
                    seq = call_bedtools_fasta(line, working_bin_fasta, binname, "16S")
                    print(seq, file=SSU_seqs)
                elif gene_name == "23S":
                    LSU=LSU+1
                    seq = call_bedtools_fasta(line, working_bin_fasta, binname, "23S")
                    print(seq, file=LSU_seqs)
        print(binname, TSU, SSU, LSU, sep="\t", file=rRNA_table)

####
rRNA_table.close()
TSU_seqs.close()
SSU_seqs.close()
LSU_seqs.close()