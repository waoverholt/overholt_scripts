#!/usr/bin/env python3
"""
Written by Will Overholt
12.9.19

This script was written to parse the default output produced by cutadapt.
It works with the output from the cutadapt_pairedend_2step.sh raw_seqs/ > cutadapt_results.txt command.

"""
import os, sys, re

in_file = sys.argv[1]

if (len(sys.argv) > 2):
    out_file = sys.argv[2]
    OUT=open(out_file, "w")

fr = ""
rr = ""
total = ""
fra = ""
rra = ""
pr = ""        
if (len(sys.argv) > 2):
    print("R1_path", "R2_path", "Total Initial Reads", "ForwardReads_w_Adapters", "ReverseReads_w_Adapters", "PairedReads_w_Adapters", sep="\t", file=OUT)
else:
    print("R1_path", "R2_path", "Total Initial Reads", "ForwardReads_w_Adapters", "ReverseReads_w_Adapters", "PairedReads_w_Adapters", sep="\t")




with open(in_file, "r") as f:
    for line in f:
        if re.match(r'Command line parameters', line):
            elems = line.split(" ")
            fr = os.path.basename(elems[elems.index('-o')+1])
            rr = os.path.basename(elems[elems.index('-p')+1])

            #print(fr, rr)
        elif re.match(r'Total read pairs processed', line):
            elems = line.split()
            total = elems[4]

        elif re.match(r'  Read 1 with adapter', line):
            elems = line.split()
            #print(elems[4], elems[5])
            fra = " ".join([elems[4], elems[5]])

        elif re.match(r'  Read 2 with adapter', line):
            elems = line.split()
            #print(elems[4], elems[5])
            rra = " ".join([elems[4], elems[5]])

        elif re.match(r'Pairs written', line):
            elems = line.split()
            #print(elems[4], elems[5])
            pr = " ".join([elems[4], elems[5]])
            
            if (len(sys.argv) > 2):
                print(fr, rr, total, fra, rra, pr, sep="\t", file=OUT)
            else:
                print(fr, rr, total, fra, rra, pr, sep="\t")


if (len(sys.argv) > 2):
    OUT.close()