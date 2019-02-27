#!/usr/bin/env python2.7
my_usage = """
Written by Will Overholt
01-14-2015

This script will wrap PEAR for use on a directory containing paired end reads.
It uses the default command and is not parallelized.

This script assumes you have copied all the sample files into one directory. It was built using FLUIDIGM names (FLOOO1...)

BE SUPER RIDICULOUSLY CAREFUL THAT YOUR SAMPLE NAMES CAN'T MATCH MULTIPLE FILE NAMES (uses REGEX matching)!!
I have added a step that makes sure only 2 fastq files match your sampleID, but this isn't full proof.

As input please provide:
    (1) folder containing all sequences
    (2) A file that has each sample listed
    (3) output directory

IMPORTANT
THERE IS A SECURITY ISSUE WHERE I USE SHELL=T TO CALL PEAR.

"""

import sys, os, re
import argparse
import subprocess as sub
import shlex
from my_useful_functions import *

parser = argparse.ArgumentParser(usage = my_usage)
parser.add_argument("-i", "--input", help="input directory that contains all forward & reverse reads for each sample", required=True)
parser.add_argument("-t", "--text_file", help="text file that contains each sample on an individual line", required=True)
parser.add_argument("-o", "--output", help="directory to save pear results")

args = parser.parse_args()

file_list = [f for f in os.listdir(args.input)]
#print file_list

make_sure_path_exists(args.output)

#Check that everything matches up first
with open(args.text_file) as mapping:
    for line in mapping:
        line = line.rstrip()
        #print line
        r = re.compile('%s.*'%line)
        matches = [m.group() for m in [r.search(f) for f in file_list] if m]
        #print matches
        if len(matches) != 2:
            sys.exit("Mutiple fastq files [n = {0}, and should be exactly 2], match your sample ID list\nExiting the script and you'll need to rerun".format(len(matches)))

with open(args.text_file) as mapping:
    for line in mapping:
        line = line.rstrip()
        #print line
        r = re.compile('%s.*'%line)
        matches = [m.group() for m in [r.search(f) for f in file_list] if m]
        r1 = re.compile(".*_R1_.*")
        read1_match = [m.group() for m in [r1.search(elem) for elem in matches] if m]
        read1 = ''.join(read1_match)
        r2 = re.compile(".*_R2_.*")
        read2_match = [m.group() for m in [r2.search(elem) for elem in matches] if m]
        read2 = ''.join(read2_match)


        read1 = os.path.join(os.path.abspath(args.input), read1)
        read2 = os.path.join(os.path.abspath(args.input), read2)
        output = os.path.join(os.path.abspath(args.output), line)
        #print read1, read2, output
        commands = '$HOME/program_files/pear-0.9.6-bin-64/pear-0.9.6-bin-64 -f {0} -r {1} -o {2}'.format(read1, read2, output)
        print commands
        inargs = shlex.split(commands)
        sub.call(commands, shell=True)
