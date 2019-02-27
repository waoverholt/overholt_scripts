#!/usr/bin/env python
import sys, os, re

infile = sys.argv[1]
size = int(sys.argv[2])

i = 0
with open(infile) as f:
    count = 0
    for line in f:
        line = line.rstrip()
        if re.match('(^@O)|(^@M)', line):
            print line
            count += 1
            i += 1
        elif re.match('^(A|T|G|C)', line):
            print line[0:size+1]
            count += 1
            i += 1
        elif re.match('\+$', line):
            print line
            count += 1
            i += 1
        elif count >2:
            print line[0:size+1]
            count = 0 
            i += 1
