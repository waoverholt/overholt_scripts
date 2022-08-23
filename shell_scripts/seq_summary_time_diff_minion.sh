#!/bin/bash

# Written by Will Overholt
# 2022-08-09

# Input: "sequencing_summary.*.txt" generated from minION run or from running guppy on a set of fast5 files
# usage: seq_summary_time_diff_minion.sh <path to summary file>

seq_output=$1

cut -f 7 "$seq_output" | \
  sed 1d | \
  sort -n | \
  sed -n '1p;$p' | \
  paste - - | \
  awk 'BEGIN {printf "Run Time\n"}; {sec=($2-$1); min=sec/60; hour=min/60}; END {printf "seconds = %s\nminutes = %s\nhours = %s\n", sec, min, hour}'

