#!/usr/bin/env python3
"""
Written by Will Overholt
13.6.2018

This script parses a kegg.json file and saves it a tab delimited.
"""
import sys, os
import json

json_data = sys.argv[1]
outfile = sys.argv[2]
OUT = open(outfile, "w")
with open(json_data, "r") as f:
    data = json.load(f)
for A in data['children']:
    for B in A['children']:
        for C in B['children']:
            try:
                for D in C['children']:
                    for key, value in D.items():
                        elems = value.split("  ")
                        KO = elems[0]
                        gene = elems[1]
                        print(KO, gene, A['name'], B['name'], C['name'], sep="\t", file=OUT)
            except KeyError:
                next
OUT.close()
