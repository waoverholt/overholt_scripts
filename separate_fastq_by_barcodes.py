#!/usr/bin/env python

"""
Written by Will Overholt
2015-06-09

This program will split a multiplexed fastq file into individual sample
files. It takes as input a fastq file (written for fastq.join.fastq),
a QIIME formatted mapping file, and a directory to save the files.

Each new file will be named based on the file name in the mapping file.
"""

import sys,os,re
import argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="fastq file to split", required=True)
parser.add_argument("-o", "--output_dir", help="output directory", required=True)
parser.add_argument("-m", "--mapping_fp", help="mapping file path with barcode and sample names", required=True)

args = parser.parse_args()

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

sample_barcodes = {}
with open(args.mapping_fp) as mapping:
    for line in mapping:
        if not line[0] == "#":
            line_elems = line.split()
            (key, val) = line_elems[0],line_elems[1]
            sample_barcodes[(key)] = val

for samp,barcode in sample_barcodes.items():
    f = open(os.path.join(args.output_dir, samp+".fastq"), 'w')
    for name, seq, qual in FastqGeneralIterator(open(args.input)):
        #header = name.split(":")
        header = name.split("=")
        #seq_barcode = header[9]
        seq_barcode = header[1]
        seq_barcode = seq_barcode.split(" ")
        if seq_barcode[0] == barcode:
            f.write("@{}\n{}\n+\n{}\n".format(name,seq,qual))
    f.close()

