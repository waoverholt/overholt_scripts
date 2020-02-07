#!/usr/bin/env python3

"""
Written by Will Overholt
4.11.19

This script was written to process a collection of outputs from multiple ROCker models for multiple samples.
It takes as input a directory, that contains directories named after each rocker gene model used.
Here is an example directory tree.

  |-ROCker
  |  |-RpoB.100
        |-T1C1.output
        |-T1C2.output
        |-T1C3.output
        |-T1O1.output
  |  |-AmoA_A_v2.125
  |  |-NosZ_v2.125
  |  |-NifH.100
  |  |-AmoA_B_v2.125
  |  |-NifH.150
  |  |-NorB.125
  |  |-Hao.125
  |  |-NirK.100

The input files were created with the scripts:
/nv/hp10/woverholt3/job_scripts/metaT_scripts/rocker_scripts
submit_multiple_rocker.sh
multiple_qsub_rocker.pbs
"""
#Imports
import os, sys
import argparse

#Argument Parser
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_dir', required=True,
                    help='directory containing directories named after each rocker model used')
parser.add_argument('-o', '--output_file', required=True,
                    help='output table summarize results for each genexsample')
args = parser.parse_args()

#Main script
indir = args.input_dir
outfile = args.output_file

OUT=open(outfile, "w")
print("GeneModel\tSampleID\tCounts", file=OUT)
#loop over each gene model
for gene_dir in os.listdir(indir):
    if os.path.isdir(os.path.join(indir,gene_dir)):
        #rocker directories are labeled "gene_name.model_size", so I keep only the "gene_name" portion
        gene_id = gene_dir.split(".")[0]
        #loop over each sample within the gene model
        for sample in os.listdir(os.path.join(indir, gene_dir)):
            if os.path.isfile(os.path.join(indir, gene_dir, sample)):
                #my rocker output names are "sample_name.output"
                sampleid = sample.split(".")[0]
                with open(os.path.join(indir, gene_dir, sample)) as f:
                    num_lines = sum(1 for line in f)
                    print(gene_id, sampleid, num_lines, sep="\t", file=OUT)
OUT.close()