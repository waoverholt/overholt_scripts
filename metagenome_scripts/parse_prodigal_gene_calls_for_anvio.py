#!/usr/bin/env python3

"""
Written by Will Overholt
03.02.2020

This script was written so I could use the same gene calls in the Probst overview files
as in Anvi'o.
It takes as input the default .faa file produced by prodigal.

Usage:
parse_prodigal_gene_calls_for_anvio.py -f gene_calls.faa -o external_gene_calls.txt
"""

import sys, os
import argparse

#Argument Parser
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fasta', required=True,
                    help='amino acid fasta file with the typical prodigal header produced from the assembly contigs / scaffolds')
parser.add_argument('-o', '--output_file', required=True,
                    help='output table for anvio')
args = parser.parse_args()

#Initiate outfile
OUT = open(args.output_file, "w")
print("gene_callers_id", "contig", "start", "stop", 
"direction", "partial", "source", "version", sep="\t", file=OUT)

#Process fasta file
with open(args.fasta) as fasta:
    iline=0
    for line in fasta:
        if line[0] == ">":
            fields = line.split(" # ")
            partial = fields[-1].split(";partial=")[1].split(";")[0]
            header = fields[0]
            gene_callers_id = iline
            contig = header[1:]
            contig = "_".join(contig.split("_")[:-1])
            start = int(fields[1])-1
            stop = fields[2]
            direction = {"1":"f", "-1":"r"}[fields[3]]
            partial = {True:"0",False:"1"}[partial == "00"]
            print(gene_callers_id, contig, start, stop, 
            direction, partial, "prodigal", "v2.6.3", sep="\t", file=OUT)
            iline=iline+1

OUT.close()