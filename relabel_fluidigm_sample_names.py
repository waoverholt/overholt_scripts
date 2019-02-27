#!/usr/bin/env python2.7
my_usage = """
This script renames a directory of files that include the fluidigm keyword in the file name. For example, it will find a file labeled FLD0001_S1_L001_R1_001.fastq and rename it Will01_.fastq.

Note, in the key step it uses regex matching to find the fluidigm keyword in the entire header. THIS SCRIPT WILL FUCK UP YOUR DAY IF YOUR FILES ARE LABELED FLD1, FLD10, FLD11 as it will match all of them with FLD1.

The user provides a tab separated two column file, the first containing the fluidigm keyword, the second containing the new name you would like.

Any special characters that break bash will break this script and give an error, eg. $, %, &, @, ...

"""

import sys,os,re
import argparse
from my_useful_functions import *
import subprocess as sub

parser = argparse.ArgumentParser(usage = my_usage)
parser.add_argument("-i", "--input_dir", help="directory containing all fluidigm labeled sequence files", required=True)
parser.add_argument("-m", "--mapping_fp", help="self made mapping file, matches the fluidigm names to a specified sample name, should be tab delimited (fluidigm name \t sample_name \n)", required=True)
parser.add_argument("-o", "--output_dir", help="where to save the new files", required = True)

args = parser.parse_args()

make_sure_path_exists(args.output_dir)

for f in os.listdir(args.input_dir):
    with open(args.mapping_fp) as mapping:
        for line in mapping:
            line = line.rstrip()
            elems = line.split(" ")
            fld_label = elems[0]
            samp_id = elems[1]

            if re.match(fld_label, f):
                in_file = os.path.join(os.path.abspath(args.input_dir), f)
                outfile = os.path.join(os.path.abspath(args.output_dir), "{0}_.fastq".format(samp_id))
                command = "cp {0} {1}".format(in_file, outfile)
                sub.call(command, shell = True)
                
                
