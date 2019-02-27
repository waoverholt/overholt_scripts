#!/usr/bin/env python3
my_usage = """
Written by Will Overholt
7 June 2018

This script is supposed to grab 2 paired-end fastq files, identify the one with the fewest sequences, and subsample each to the same depth (specified in decimal format, 0.1 = 10%).

This script checks to make sure there are only 2 files per pair. If it finds more than 2 it means you need to adjust the file_regex variable for your specific sample naming scheme.

Script usage:
subsample_PE_fastq.py -i/--input-dir -o/--output-dir -f/--fraction-subsample
(1) Input directory containing all fastq paired-end reads for each sample
(2) Location to save subsampled fastq files
(3) Faction to subsample, e.g. 0.1 is 10%

This script subsamples WITHOUT replacement
"""
import sys, os, re, math, random
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description=my_usage)
parser.add_argument("-i", "--input_dir", help="input directory containing all forward and reverse reads for each sample", required=True)
parser.add_argument("-o", "--output_dir", help="output directory to save all randomly subsampled fastq files", required=True)
parser.add_argument("-f", "--fraction_subsample", help="decimal value for proportion to subsample to, 0.1 = 10% subsampled.", required=True)

args = parser.parse_args()

def blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b: break
        yield b

file_regex = re.compile(r"^(.*?)_.*")

INDIR = args.input_dir
OUTDIR = args.output_dir
PERC_SUB = args.fraction_subsample

perc_for_name = str(math.floor(float(PERC_SUB)*100))+"perc"
file_list = [f for f in os.listdir(INDIR)]
sample_names = [m.group(1) for f in file_list for m in [file_regex.search(f)] if m]
set_name = sorted(set(sample_names))

#Check only 2 files per pair
for uniq in set_name:
    r = re.compile('%s.*'%uniq)
    matches = [m.group() for m in [r.search(f) for f in file_list] if m]
    if len(matches) != 2:
        sys.exit("Mutiple fastq files [n = {0}, and should be exactly 2], match your input directory\nExiting the script and you'll need to rerun".format(len(matches)))
    #Grab matching files
    r1 = re.compile(".*_R1_.*")
    r2 = re.compile(".*_R2_.*")
    read1_match = [m.group() for m in [r1.search(elem) for elem in matches] if m]
    read2_match = [m.group() for m in [r2.search(elem) for elem in matches] if m]
    read1 = ''.join(read1_match)
    read2 = ''.join(read2_match)
    
    read1_fp = os.path.join(os.path.abspath(INDIR), read1)
    read2_fp = os.path.join(os.path.abspath(INDIR), read2)

    read1_basename = os.path.splitext(os.path.basename(read1_fp))[0]
    read2_basename = os.path.splitext(os.path.basename(read2_fp))[0]
    
    r1_len = 0
    r2_len = 0
    with open(read1_fp, "r", encoding="utf-8", errors="ignore") as f:
        r1_len = sum(bl.count("\n") for bl in blocks(f))

    with open(read2_fp, "r", encoding="utf-8", errors="ignore") as f:
        r2_len = sum(bl.count("\n") for bl in blocks(f))
    print("Forward File: {0}\tNumber of Lines: {1}".format(read1, r1_len))
    print("Reverse File: {0}\tNumber of Lines: {1}".format(read2, r2_len))
    min_len = min(r1_len, r2_len)
    sub_samp_num = math.floor(min_len * float(PERC_SUB) / 4)

    #generate seq subsample list
    random_list = random.sample(range(0, math.floor(min_len / 4)), sub_samp_num)
    print("Number of Sequences Sampled: {0}".format(len(random_list)))

    #select fasta lines mapping to the list
    outfile1 = os.path.join(os.path.abspath(OUTDIR), "{0}_{1}.fastq".format(read1_basename, perc_for_name))
    OUT1 = open(outfile1, "w")
    with open(read1_fp, "r") as f:
        for i,seq in enumerate(SeqIO.parse(f, "fastq")):
            #print(i, seq.id)
            if i in random_list:
                OUT1.write(seq.format("fastq"))
    OUT1.close()

    outfile2 = os.path.join(os.path.abspath(OUTDIR), "{0}_{1}.fastq".format(read2_basename, perc_for_name))
    OUT2 = open(outfile2, "w")
    with open(read2_fp, "r") as f:
        for i,seq in enumerate(SeqIO.parse(f, "fastq")):
            #print(i, seq.id)
            if i in random_list:
                OUT2.write(seq.format("fastq"))
    OUT2.close()
    print()



