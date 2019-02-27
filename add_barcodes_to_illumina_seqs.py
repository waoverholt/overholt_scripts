#!/usr/bin/env python2.7

import sys, os, re
import argparse
from my_useful_functions import *

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="directory containing pear joined files", required=True)
parser.add_argument("-m", "--mapping_fp", help="self made mapping file, uses the fluidigm barcodes to match sequences to sample names, should be tab delimited (fluidigm num \t barcode \n)", required=True)
parser.add_argument("-o", "--output_fp", help="output directory", required = True)

args = parser.parse_args()
make_sure_path_exists(args.output_fp)

file_list = [f for f in os.listdir(args.input)]
#print file_list

with open(args.mapping_fp) as mapping:
    line = mapping.readline()
    line.rstrip()
    elems = line.split("\t")
    fld = elems[0]
    barcode = elems[1]
    reg = re.compile('%s.*'%fld)
    matches = [m.group() for m in [reg.match(f) for f in file_list] if m]
    match = ''.join(matches)
    #print match

    cur_file = os.path.join(os.path.abspath(args.input), match)
    with open(cur_file) as f:
        out = open(os.path.join(os.path.abspath(args.output_fp), "{0}_barcoded.fastq".format(fld)), 'w')
        for line in f:
            line = line.rstrip()
            if re.match("@M[0-9]+:.*", line):
                out.write("{0}#{1}".format(line,barcode))
            else:
                out.write(line)
        out.close()
