#! /usr/bin/env python3
# %%
import json
import os,sys,re
import argparse
import pandas as pd
import random
from collections import defaultdict
import matplotlib.pyplot as plt

my_usage = """
This script was written to supplement Miguel's script JPlace.to_iToL.rb script. It is meant to take a .jplace file produced by RaxML-EPA.

The idea is to make a good tree with long sequences, then map illumina OTUs onto the tree. This script will ID the OTUs that map to the nodes on the phylogenetic tree so you can inflate the OTUs back into counts.

The end goal is to use the iToL (or ggtree) to annotate the reference tree with relative abundance of reads associated with the nodes.

It takes as input (-i) a .jplace file, a counts file, and output (-o) a file name.

I've only tested this on one dataset. The OTUIDs are called denovo[0-9]+ and the BIOM table has the same row names. I use pandas to split the OTU table based on sequences that are associated with the same node.

Updated to be more flexible for the OTUIDs
User asked to provide the text string format of the otuIDs:
e.g. "ASV_1, "ASV_2", "ASV_3,...ASV_999" -> "ASV_"
"denovo0, denovo1, ..., denovon" -> "denovo"
"""

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help=".jplace file from raxml-epa", required=True)
parser.add_argument("-o", "--output", help="output file name to store results")
parser.add_argument("-c", "--counts", help="species counts table that is tab separated")
parser.add_argument("-s", "--otu_string_format", help = "text format for the otuids, include all characters up to the number")

args = parser.parse_args()

otu_match = re.compile(re.escape(args.otu_string_format) + r'[0-9]+')

#Delete output file if it already exists
try:
    os.remove(args.output)
except OSError:
    pass

with open(args.input) as jplace:
    json_obj = json.load(jplace)

## This is stepping through the reference tree and pulling out the internal nodes
Tree = json_obj['tree']
Tree = re.split(r"\(|\)|,", Tree)
node_names = []
for char in Tree:
    if char:
        m = re.match("(.*):[0-9]+\.[0-9]+({[0-9]+})", char)
        if m:
            if len(m.groups()) > 1:
                node_names.append([m.group(2), m.group(1)])
            elif len(m.groups()) == 1:
                node_names.append([m.group(1)])
        else:
            next
    
node_dict = defaultdict(list)
for item in json_obj["placements"]:
    node_key = item['p'][0][0]
    otu_value = item['n'][0]
    m = re.search(otu_match, otu_value)
    node_dict[node_key].append(m.group())


#Subsetting the OTU table and summarizing node results
otu_df = pd.read_table(args.counts, sep="\t", skiprows=(0), header=(1), index_col=(0))

new_df = pd.DataFrame(index = node_dict.keys(), columns = list(otu_df))
#print otu_df[otu_df.index.isin(['denovo0', 'denovo11'])]

for node, otu in node_dict.items():
    line = "{0}\t{1}\n".format(node, otu)
    #print(otu)
    #print(otu_df[otu_df.index.isin(otu)])
    #print(otu_df.sum())
    new_df.loc[node] = otu_df[otu_df.index.isin(otu)].sum()

new_index_names = []
for index in new_df.index:
    for nodeid in node_names:
        newid = re.match("\{([0-9]+)\}", nodeid[0])
        if int(newid.group(1)) == int(index):
            if nodeid[1]:
                #new_df.index.rename = nodeid[1]
                new_index_names.append(nodeid[1])
            else:
                new_index_names.append("{{{0}}}".format(index))


new_df.index = new_index_names
new_df = new_df.assign(PiePlace=0.5)
new_df = new_df.assign(RowSum=new_df.sum(axis=1))
col_list = new_df.columns.tolist()
col_list = col_list[-2:] + col_list[:-2]
new_df = new_df[col_list]

#format database for iTol
f = open(args.output, 'a')
f.write("DATASET_PIECHART\n")
f.write("SEPARATOR TAB\n")
f.write("DATASET_LABEL\tReadPlacement\n")
f.write("COLOR\t#1f2122\n")
f.write("FIELD_LABELS\t")
f.write("\t".join(map(str, list(new_df)[2:])))
f.write("\n")
f.write("#I am being a bit lazy here, but I haven't spent the time to generate 'good' colors for the samples. I recommend using the webiste http://tools.medialab.sciences-po.fr/iwanthue/examples.php or http://phrogz.net/css/distinct-colors.html to get a list of colors that you like. They should be tab separated.\n")
f.write("FIELD_COLORS\n")
f.write("DATA\n")
new_df.to_csv(f, sep="\t", header = False)
f.close()

# %%
