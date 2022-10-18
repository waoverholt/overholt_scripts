#!/bin/bash

# Written by Will Overholt
# 2022-07-28

# Input: a multiline fastq file
# usage: calc_n50+reads+bases_fastq.sh <path to final_summary file>

fastq=$1

echo -e "n50 num_reads total_bases" > tmp_out 

awk '{getline;print length($0); s += length($1); getline; getline;} END{print "+"s}' $fastq | sort -gr | awk 'BEGIN{bp = 0; f = 0} {if(NR == 1){sub(/\+/, "", $1);s=$1} else{bp += $1; if(bp > s / 2 && f == 0) {n50 = $1; f=1}}} END {print n50, (NR - 1), s}' >> tmp_out

column -t tmp_out

rm tmp_out
