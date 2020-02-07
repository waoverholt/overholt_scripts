#!/usr/bin/env python3

import os, sys
import argparse

####################
#Argument Parser
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_dir', required=True,
                    help='directory of bin directories, each containing a "genes.faa" file, typically generated from checkM. This can be modified with -t')
parser.add_argument('-o', '--output_file', required=True,
                     help='concatenated file of all gene sequences from all bins, with their sequence names renamed')
parser.add_argument('-f', '--file_name', required=False, nargs='?', default="genes.faa", type=str,
                    help='The name of the protein sequence file found within each bin')
parser.add_argument('-d', '--delimiter_bin', required=False, nargs='?', default=".", type=str, 
                    help='The bin name delimiter. CheckM is usually a period while anvio used _')
                    
args = parser.parse_args()
#####################
o_bins = [d for d in os.listdir(args.input_dir) if "bin" in d]
o_bins.sort(key=lambda x: int(x.split(args.delimiter_bin)[1]))

outfile=open(args.output_file, 'w')
for bin in o_bins:
    with open(os.path.join(args.input_dir, bin, args.file_name)) as f:
        i = 1
        for line in f:
            if line[0] == ">":
                line_elems = line.rstrip().split(" # ")
                line_id = line_elems[0][1:]
                print(">" + bin + "_" + str(i) + "__" + line_id, file=outfile)
                i = i+1
            else:
                print(line, end='', file=outfile)

outfile.close()