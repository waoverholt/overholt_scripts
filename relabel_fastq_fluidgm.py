#!/usr/bin/env python
my_usage = """
Written by Will Overholt
6-29-2015

This is script was written to facillitate using QIIME's split_libraries_fastq.py function with GATech's illumina sequencer output. The raw files are placed in individual files with NO BARCODES.

This script assumes you have already merged and then joined all the R1 and R2 reads. The GATech sequence attaches a unique sample identifier to the end of the sequence fastq header.

E.g @M01562:92:000000000-AFV2M:1:1101:14041:1874 1:N:0:142
The 142 at the end of 1:N:0:142 is unique to this sample. It refers to the fluidigm barcode number FLD0142).

As input please provide:
    (1) the joined fastq file
    (2) A self made mapping file that is tab separated:
        FLD0142    samp1    barcode
        FLD0143    samp2    barcode
        FLD1120    wao_samp1_oil    barcode
    (3) A new file name to save the relabeled sequences

This script will rename each sequence according the the sample name provided in the second column of the mapping file. Sequences will be numbered consecutively starting at 0 for each sample
"""

import sys,os,re
import argparse

parser = argparse.ArgumentParser(usage = my_usage)
parser.add_argument("-i", "--input", help="concatenated fastq file that has all the sample files merged into one", required=True)
parser.add_argument("-m", "--mapping_fp", help="self made mapping file, uses the fluidigm barcodes to match sequences to sample names, should be tab delimited (fluidigm num \t sample_name \t barcode \n)", required=True)
parser.add_argument("-o", "--output_fp", help="output name new fastq file", required = True)

args = parser.parse_args()

#Get the unique identifier of the Illumina machine
#This will be used to parse the sequences and ensure we don't accidentally parse a quality line instead of a header line
f = open(args.input)
first_line = f.readline()
line_elements = first_line.split(":")
machine_name = line_elements[0]
f.close()

out = open(args.output_fp, 'w')

counter = 0
sample_list = []
with open(args.mapping_fp) as mapping:
    for line in mapping:
        line = line.rstrip()
        elems = line.split("\t")
        samp_id = re.search('[0-9]+', elems[0])
        samp_line = (samp_id.group(), elems[1], elems[2])
        sample_list.append(samp_line)
        
cur_sample = sample_list[0][1]        
with open(args.input) as inseqs:
    for line in inseqs:
        if re.match(machine_name, line):
            #I'm assuming the sample ID is the last element in the sequence name, and fields are separated by colons
            seq_name_elems = line.split(":")
            unique_samp_num = seq_name_elems[len(seq_name_elems)-1]
            for samp in sample_list:
                #print samp[0], "\t", unique_samp_num
                if int(unique_samp_num) == int(samp[0]):
                    new_samp_name = re.sub(unique_samp_num, samp[2], line)
                    out.write("{0}\n".format(new_samp_name))
                    if cur_sample == samp[1]:
                        counter = counter + 1
                    else:
                        counter = 0
                        cur_sample = samp[1]
        else:
            out.write(line)

out.close()

        
