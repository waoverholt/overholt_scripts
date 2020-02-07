#!/usr/bin/env python3

"""
Written by Will Overholt
05.02.2020

This script was written so I could use the same gene calls in the Probst overview files
as in Anvi'o.
It takes as input the default .faa file produced by prodigal.

Usage:
parse_prodigal_gene_calls_for_anvio.py -f gene_calls.faa -o external_gene_calls.txt
"""

import sys, os
import argparse
import numpy as np
import pandas as pd


#Argument Parser
parser = argparse.ArgumentParser()
parser.add_argument('-p', '--prodigal_genes', required=True,
                    help='The result of the script: "parse_prodigal_gene_calls_for_anvio.py"')     
parser.add_argument('-t', '--taxonomy', required=True,
                    help='The probst file *.scaff2tax.txt /n scaffold/taxonomy_string')
parser.add_argument('-o', '--output_file', required=True,
                    help='output table for anvio')
args = parser.parse_args()


#Create dictionary of the taxonomy
tax_dict = {}
with open(args.taxonomy) as tax:
    for line in tax:
        line = line.rstrip()
        scaf, tax_str = line.split("\t")
        tax_dict[scaf] = tax_str


#Parse the gene file

OUT = open(args.output_file, 'w')
print("gene_callers_id", "t_domain", "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species", sep="\t", file=OUT)
with open(args.prodigal_genes) as genes:
	header = genes.readline()
	for line in genes:
		fields = line.split("\t")
		tax_str = tax_dict[fields[1]]
		tax_fields = tax_str.split(";")
		fixed_tax = []
		for i in range(7):
			try:
				fixed_tax.append(tax_fields[i])
			except IndexError:
				fixed_tax.append("unclassified")
		print(fields[0], "\t".join(fixed_tax), sep="\t", file=OUT)
		#print(fields[0], tax_str, sep="\t")

OUT.close()
