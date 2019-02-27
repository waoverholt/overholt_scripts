#!/usr/bin/env python2.7
my_usage = """
This script looks through a directory of fastq files. It was written to copy those files and rename them to a user specified name. I had originally meant it to rename default fluidigm sample names. For example, it would copy FLD0001_S1_L001_R1_001.fastq and relabel it as Will01_.fastq.

It leaves in the trailing underscore(_) after the name so QIIME's multiple_split_libraries.py script recognizes your sample name based on the underscore.

You need to provide a mapping file, which is 2 columns separated by a tab. The first column should be the first element (separated by underscores, eg FLD0001) of the file name and the second column is the the name you want to change it to.

FLD0001\tWill01
FLD0002\tWill02
...


THIS SCRIPT WILL BREAK IF YOUR SAMPLE NAMES CONTAIN SPECIAL CHARACTERS (such as &, %, $, etc...). Basically anything that will break bash.
"""

import sys,os,re
import argparse
from my_useful_functions import *
import subprocess as sub

parser = argparse.ArgumentParser(usage = my_usage)
parser.add_argument("-i", "--input_dir", help="directory containing all fluidigm labeled sequence files", required=True)
parser.add_argument("-m", "--mapping_fp", help="self made mapping file, matches the current file names to a specified sample name, should be tab delimited (sample_file_name \t desired_sample_name \n)", required=True)
parser.add_argument("-o", "--output_dir", help="where to save the new files", required = True)

args = parser.parse_args()

make_sure_path_exists(args.output_dir)

for f in os.listdir(args.input_dir):
    sample_file_elements = f.split("_")
    sample_file_id = sample_file_elements[0]
    
    with open(args.mapping_fp) as mapping:
        for line in mapping:
            line = line.rstrip()
            elems = line.split("\t")
            fld_label = elems[0]
            samp_id = elems[1]

            if (fld_label == sample_file_id):
                in_file = os.path.join(os.path.abspath(args.input_dir), f)
                outfile = os.path.join(os.path.abspath(args.output_dir), "{0}_.fastq".format(samp_id))
                command = "cp {0} {1}".format(in_file, outfile)
                sub.call(command, shell = True)
                
                
