#!/usr/env python3

import os, sys
from datetime import datetime

## needs to be sorted
"""
# convert from fastq to all info on single line, tab delim
cat <fastq> | paste - - - - > out.txt

# sort on the "start_time" field (happens to be column 6 for this dataset)
sort -t"=" -k6,6 out.txt > sorted_out.txt
"""

in_file = 'test_split/sorted_06031011_1line.txt'

start_time = '0'

breakpoints = [30,45,60,120,180,240,300,360]

for cut_time in breakpoints:
    out_fp = f"test_split/{cut_time}min.fastq"
    out_f = open(out_fp, 'w')

    with open(in_file, 'r') as f:
        for line in f:
            header, seq, div, qual = line.split("\t")
            time = header.split("=")[5]
            time = datetime.strptime(time,"%Y-%m-%dT%H:%M:%SZ")
            if start_time == '0':
                start_time = time
            time_diff = time - start_time
            min_diff = time_diff.total_seconds()/60

            if min_diff < cut_time:
                out_f.write(f"{header}\n{seq}\n{div}\n{qual}")
            else:
                print(f"finished {cut_time}")
                out_f.close()
                break

            #print(time_diff.total_seconds()/60)
        