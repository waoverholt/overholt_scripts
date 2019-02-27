#!/usr/bin/env python

"""
Written by Will Overholt
8-8-2013

this script takes a qiime biom file and performs a transformation
options are:
(1) square root transformation (default)
(2) forth root transformation
(3) log(x+1) transformation

If there is taxonomy please pass the -x flag

"""

from biom.parse import parse_biom_table
import argparse
import subprocess
import tempfile
import os
import math
import re

#defining functions for transformation that will be used later
def square_root(x):
    return math.sqrt(x)
def forth_root(x):
    return math.sqrt(square_root(x))
def log(x):
    return math.log(x+1,10)

#set up arguments for script
#currently, -i for input, -o for output, -t for transformation
parser = argparse.ArgumentParser(description='This simple script does a few transformations on a biom style OTU table as produced by qiime. You can currently do a squareroot (sqrt), 4throot (4rt), or a log(x+1) (log) transformation to downweight abundant OTU. If there is taxonomy associated with the biom file please pass the -x flag')
parser.add_argument("-i", "--input", help="path to input biom file", required=True)
parser.add_argument("-o", "--output", help="output file path", required=True)
parser.add_argument("-t", "--transform", help="transformatin you wish to do, options are sqrt, 4rt, log", default='sqrt')
parser.add_argument("-x", "--taxonomy", help="does the OTU table contain a taxonomy file? Pass this flag if it does", action='store_true')

args = vars(parser.parse_args())

#set up variables based on users arguments
input = args['input']
output = args['output']
trans = args['transform']
taxonomy = args['taxonomy']

#bulky way to do this, but I convert the .biom file to a .txt file using a temp
#then perform the transformation and finally convert it back to a .biom file

temp_otu = tempfile.NamedTemporaryFile(delete=False)
temp_otu2 = tempfile.NamedTemporaryFile(delete=False)

#convert to .txt file
p = subprocess.call(['convert_biom.py', '-i', input, '-o', temp_otu.name, '-b', '--header_key', 'taxonomy'])
f = open(temp_otu.name)
f2 = open(temp_otu2.name, 'w')

#here's where the work actually gets done
for line in f:
    line.rstrip() #remove newlines / white space
    if re.search('#', line): #skip the headers
        f2.write(line)
    else:
        fields = line.split('\t') 
        sample = []
        sample.append(fields[0]) #don't do anything to OTU id
        num_e = len(fields)
        trans_value = None
        if taxonomy == True:
            for e in fields[1:num_e-1]:
                v = float(e)
                if trans == 'sqrt':
                    trans_value = square_root(v)
                elif trans == '4rt':
                    trans_value = forth_root(v)
                elif trans == 'log':
                    trans_value = log(v)
                sample.append(trans_value)
            sample.append('\t'.join(map(str,fields[num_e-1:num_e])).rstrip())

            f2.write('\t'.join(map(str,sample)))
            f2.write('\n')
        elif taxonomy == False:
            for e in fields[1:num_e]:
                v = float(e)
                if trans == 'sqrt':
                    trans_value = square_root(v)
                elif trans == '4rt':
                    trans_value = forth_root(v)
                elif trans == 'log':
                    trans_value = log(v)
                sample.append(trans_value)
            
            f2.write('\t'.join(map(str,sample)))
            f2.write('\n')
f.close()
f2.close()
p2 = subprocess.call(['convert_biom.py', '-i', temp_otu2.name, '-o', output, '--biom_table_type=otu table', '--process_obs_metadata', 'taxonomy'])

os.unlink(temp_otu.name)
os.unlink(temp_otu2.name)

