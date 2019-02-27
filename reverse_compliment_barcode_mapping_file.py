#!/usr/bin/env python

import sys,os,re
from Bio.Seq import Seq

mapping_fp = sys.argv[1]
out = open(sys.argv[2], "w")

with open(mapping_fp) as mapping:
    for line in mapping:
        if line[0] == "#":
            out.write(line)
        else:
            line_contents = line.split("\t")        
            barcode_seq = Seq(line_contents[1])
            new_line = re.sub(line_contents[1], barcode_seq.reverse_complement().strip(), line)
            out.write("{0}\n".format(new_line))
            
        
