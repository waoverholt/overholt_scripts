#!/usr/bin/env python

import os
import sys
import re

d = sys.argv[1]
path = os.listdir(d)

sample=[]
results=[]
data = []

f = open(os.path.join(d, path[0]))
next(f)         
for line in f:
    fields = line.split('\t')
    sample.append([fields[0], float(fields[1]), float(fields[2].rstrip())])
#print sample
f.close()
count = 0
for fname in path[1:]:
    f = open(os.path.join(d, fname), 'r')
    f.readline()
    for line in f:
        fields = line.split('\t')
        for row in sample:
            if re.match(fields[0], row[0]):
                sum1 = (row[1])+float(fields[1])
                sum2 = (row[2])+float(fields[2])
                results.append([row[0], sum1, sum2])
num = len(path)
final_results=[]

for row in results:
    avg1 = row[1] / float(num)
    avg2 = row[2] / float(num)
    final_results.append([row[0], avg1, avg2])

print final_results
#print sample
