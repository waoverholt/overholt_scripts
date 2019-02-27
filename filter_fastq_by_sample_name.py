#!/usr/bin/env python
import sys
from itertools import islice
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import re

def isa_group_separator(line):
    return line=='@'

input_fp = sys.argv[1]
sample_list = sys.argv[2]
output_fp = sys.argv[3]

out = open(output_fp, "w")

mapping = [line.rstrip() for line in open(sample_list)]
map_set = set(mapping)

for name, seq, qual in FastqGeneralIterator(open(input_fp)):
    sample_match = re.match("([^_]+)", name)
    sample_id = sample_match.group(1)
    if sample_id in map_set:
        out.write("@"+name+'\n'+seq+'\n'+"+"+'\n'+qual+'\n')
out.close()



