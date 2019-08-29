#!/usr/bin/env python3
help_message = """
Written by Will Overholt
28 August 2019

Purpose: Get a taxonomic consensus from a list of Uniprot protein accession numbers.
It requires (2) input files.
The first is a user provided text file, with each line representing a single or semi-colon separated group of accession numbers:
A0A0T2YIR2
A0A3D5VRX5
A0A0D7KBM4;A0A2W5IXE0

The second file is an Accession number to taxonomy file that is generated from the first file using uniprot (this can be automated in the future if necessary).
To generate this file:
1. Go to https://www.uniprot.org/uploadlists/
2. Provide the first file as input & hit Submit
3. Click the "Columns" button above the table.
4. The order of columns should to be: Entry // Taxonomic Lineage (SUPERKINGDOM) // Taxonomic Lineage (PHYLUM) // Taxonomic Lineage (CLASS) // 
Taxonomic Lineage (ORDER) // Taxonomic Lineage (FAMILY) // Taxonomic Lineage (GENUS) // Taxonomic Lineage (SPECIES)
(Note: The species column (or any) can be removed without issue. The Entry column needs to be present)
5. Hit Save & OK if there is an error (shouldn't matter)
6. Click Download
7. Change the Format to Tab-separated
8. Change to "Uncompressed"
9. Save

If no output file is provided, the results will print to the terminal.

usage:
python protein_consensus_taxonomy.py -i <user_protein_groups> -t <uniprot_taxomony_file>
                               -o <optional: output txt file>
"""

####################
#Imports
import os,sys,re
import argparse

####################
#Argument Parser
parser = argparse.ArgumentParser(help_message)
parser.add_argument('-i', '--input_protein_groups', required=True,
                    help='text file with each line representing a single or semi-colon separated group of accession numbers')
parser.add_argument('-t', '--taxonomy', required=True,
                    help='downloaded uniprot taxonomy file')
parser.add_argument('-o', '--output_file', required=False,
                     help='a tab delimited file of the results, otherwise results will be printed to terminal')
                
args = parser.parse_args()

if args.output_file:
    outFile = open(args.output_file, "w")

#####################
#Main Script

#Loop over the uniprot taxonomy and generate a dictionary 
# of the taxonomic groups for each accession
taxDict = {}
with open(args.taxonomy, "r") as tax:
    for line in tax:
        elems = line.rstrip().split("\t")
        taxDict[elems[0]] = elems[1:]

#Loop over each protein / protein group
with open(args.input_protein_groups, "r") as f:
    missing_accessions = []
    for line in f:
        IDs = line.rstrip().split(";")
        # If there is only a single protein accession number in this line
        if len(IDs) == 1:
            if args.output_file:
                print(IDs[0], ",".join(taxDict[IDs[0]]),sep="\t", file=outFile)
            else:
                print(IDs[0], ",".join(taxDict[IDs[0]]), sep="\t")
        
        #If there are multiple accession nubmers per line
        #Now we need to find the consensus taxnomoy string
        else:
            #This part borrowed heavily from the script by Tony Walterst
            #https://gist.github.com/walterst/bd69a19e75748f79efeb
            
            #List of the taxonomy strings of all accession numbers of the group
            common_acc = []
            for ACC in IDs:
                try:
                    common_acc.append(taxDict[ACC])
                except KeyError:
                    #print("Missing an Accession number")
                    missing_accessions.append(ACC)
            #Dynamically determine how many taxonomic groups we have
            levels = len(common_acc[0])
            taxa_string = []
            for n in range(levels):
                cur_taxa = []
                #loop over each of the taxonomic levels of each of the accession numbers
                for tax in common_acc:
                    cur_taxa.append(tax[n])
                #If these strings are identical then we have a match at that level (match = T)
                match = all(x==cur_taxa[0] for x in cur_taxa)
                if match:
                    taxa_string.append(cur_taxa[0])
                else:
                    taxa_string.append("")
            #Format the consensus taxonomy for each group
            final_taxa_string = ",".join(taxa_string)
            if args.output_file:
                print(";".join(IDs), final_taxa_string, sep="\t", file=outFile)
            else:
                print(";".join(IDs), final_taxa_string, sep="\t")

    if len(missing_accessions) > 0:
        print("\nThe following accession numbers were missing from your uniprot taxonomy file and")
        print("therefore were not used to determine the consensus. You may want to fix your input files to fix this.")
        for ACC in missing_accessions:
            print(ACC)

outFile.close()