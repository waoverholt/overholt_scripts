#!/usr/bin/env python
import sys, os, re

infile = sys.argv[1]
size = int(sys.argv[2])

IN = open(infile)

while 1:
    line = IN.readline()
    line = line.rstrip()
    if len(line) == 0:
        break
    count = 0
    elif re.match('@M', line):
        print line
        count += 1
    elif re.match('(A|T|G|C)', line):
        print line[0:size+1]
        count += 1
    elif re.match('\+$', line):
        print line
        count += 1
    elif count == 3:
        print line[0:size+1]
        count = 0 
IN.close()
