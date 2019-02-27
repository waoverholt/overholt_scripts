#!/usr/bin/env python
import sys,os
import re
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", help="the otu file created using the make-otus.sh script found in htc-clust", required=True)
parser.add_argument("-o", "--output_file", help="the output file in qiime otu table format", required=True)
args = parser.parse_args()

INFILE = os.path.abspath(args.input_file)
OUTFILE = os.path.abspath(args.output_file)
out = open(OUTFILE, "w")

otus = defaultdict(list)
#test = ''
with open(INFILE) as f:
    for line in f:
        line = line.rstrip()
        if not re.match("#", line):
            if re.match(">", line):
                (otu, size) = line.split()
                otuid = re.search("[0-9]+", otu)
                #test += ("\n{0}\t".format(otuid.group()))
                out.write("\n{0}\t".format(otuid.group()))
            else:
                out.write("{0}\t".format(line))
                #test += ("{0}\t".format(line))
#print OUTFILE test
out.close()
