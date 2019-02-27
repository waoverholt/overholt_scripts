#!/usr/bin/env python2
"""
This script was written as a way to recover only the 16S portion of sequences that may have much longer sequences. It assumes that you have already blasted your long sequences against a reference full length 16S sequence only (the query - ie. E. coli). The blast results should be saved in an xml file. This script then parses that xml file to recover the subject sequence name, and alignment that has had its gaps removed.

I orginally wrote this sequence since many plant pathogen sequences contained the ITS gene & tRNA genes, and would not be aligned with PyNAST.
"""


import sys, os, re
from Bio.Blast import NCBIXML

IN = sys.argv[1]

result_handle = open(IN)

blast_records = NCBIXML.parse(result_handle)

for blast_record in blast_records:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            print(">{0}".format(alignment.hit_def))
            aligned_seq = hsp.sbjct
            rm_gap = re.sub("-", "", aligned_seq)
            print(rm_gap)


