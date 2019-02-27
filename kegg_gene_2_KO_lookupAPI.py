#!/usr/bin/env python3
import sys, os, re
from Bio.KEGG.REST import kegg_link

#test = kegg_link("ko", "sly:544114")
#print(test.read())

Kgene_file_from_Uniprot = sys.argv[1]
outfile = sys.argv[2]

OUT = open(outfile, "w")
num_processed = 0
with open(Kgene_file_from_Uniprot) as fp:
    for i, line in enumerate(fp):
        if i != -120:
            line = line.rstrip()
            elems = line.split("\t")
            KO_values = elems[1].split(":")
            tax = KO_values[1]
            num = KO_values[2]
            lookup = "{0}:{1}".format(tax, num)
            conn = kegg_link("ko", lookup)
            for line in conn.readlines():
                line = line.strip()
                if line :
                    ko = line.split("\t")
                    print(elems[0], ko[0], ko[1], sep="\t", file=OUT)
        if i % 10000 == 0:
            num_processed = i
            print("Number Processed: ", num_processed, sep="")
out.close()
