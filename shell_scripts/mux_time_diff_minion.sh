#!/bin/bash

# Written by Will Overholt
# 2022-08-09

# Input: "mux_scan_data_XXXXX_XXX.csv" generated from minION run
# usage: mux_time_diff_minion.sh <path to mux file>

mux_output=$1

cut -d, -f 4 "$mux_output" | \
  sed 1d | \
  sort -n | \
  sed -n '1p;$p' | \
  paste - - | \
  awk 'BEGIN {printf "Run Time\n"}; {sec=($2-$1); min=sec/60; hour=min/60}; END {printf "seconds = %s\nminutes = %s\nhours = %s\n", sec, min, hour}'

