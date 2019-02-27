#!/usr/bin/env python

import os,sys,re
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_map", help="otu mapping file (from pick_otus.py)", required=True)
parser.add_argument("-s", "--list_of_otus", help="a list of OTUs that you want to filter the otu mapping file by", required=True)
parser.add_argument("-o", "--output_map", help="outfile file path")

args=parser.parse_args()

otu_map = [line.rstrip() for line in open(args.input_map, 'r')]

for otu in otu_map:
    print otu
    
"""
with open(sys.argv[2]) as fp:
    for line in fp:
        line = line.rstrip()
        for otu in otu_map:
            r = "^"+line+"\t"
            if re.match(r,otu):
                print otu
                break
        
"""        
