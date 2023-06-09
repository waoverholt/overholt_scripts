#!/bin/bash

# Written by Will Overholt
# 2022-07-28

# Input: genome assembly file (fasta)
# usage: script.sh <path to fastq> 

fasta=$1

#mamba activate kraken2
KRAKENDB=/home/sharedFolder/referenceData/minikraken2_v2_8GB_201904_UPDATE/

##Run kraken2
INPATH="$(dirname $fasta)"
kraken2 --db $KRAKENDB/ --report "$INPATH/kraken2_report.txt" $fasta > "$INPATH/kraken2.out"

echo
echo
echo
##Calculate assembly length
echo "The assembly length was:"
cat $fasta | paste - - | cut -f2 | awk 'BEGIN {sum=0} {sum+=length($0)} END {print sum}'

##Calculate %human (9606)
echo "The number of reads mapping to human was:"
cut -f3 "$INPATH/kraken2.out" | grep 9606 | wc -l

echo "The bp mapping to human was:"
awk 'BEING {sum=0} $3 == 9606 {sum+=$4;} END {print sum}' "$INPATH/kraken2.out"

