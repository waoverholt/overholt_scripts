#!/bin/bash

# Written by Will Overholt
# 2022-07-28

# Input: a multiline fastq file
# usage: script.sh <path to fastq> 

fastq=$1

# Change these to your use case
maxlen=7000 #assemblying plasmids
size=6000 #expecting 6kb
coverage=200

output=$(basename $fastq)
output="${output%.*}_downsampled.fastq"

num_bp=$(( $size*$coverage ))

cat $fastq | paste - - - - | awk '{printf "%d %s\n", length($2), $0}' | sort -n -k1,1 | sed -E -e 's/^[0-9]+ //' | awk -F "\t" 'BEGIN {total_len=0} length($2) < '"$maxlen"' {print $0; total_len+=length($2); if (total_len > '"$num_bp"') {exit}}' | sed -E "s|\t|\n|g" > $output

