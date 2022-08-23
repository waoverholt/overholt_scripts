#!/bin/bash

# Written by Will Overholt
# 2022-07-12

# Input: "final_summary_XXXXX_XXX.txt" generated from minION run
# usage: time_diff_minion.sh <path to final_summary file>

final_report=$1

start_time=$(grep "started" $final_report | cut -f2 -d= | sed -E "s/\-[0-9]{2}:[0-9]{2}$//")
end_time=$(grep "acquisition_stopped" $final_report | cut -f2 -d= | sed -E "s/\-[0-9]{2}:[0-9]{2}$//")


difference_s=$(( $(date -d "$end_time" "+%s") - $(date -d "$start_time" "+%s") ))

echo "difference in minutes:"
echo "$difference_s/60" | bc 
echo "difference in hours:"
echo "scale=2; $difference_s/3600" | bc
