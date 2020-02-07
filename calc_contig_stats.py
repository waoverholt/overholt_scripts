#!/usr/bin/env python

import sys,os,re
import argparse
import subprocess as sub

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", help="fasta formatted file containing the contigs", required=True)
parser.add_argument("-f", "--forward_reads", help="forward reads used to generate the assembly, this will be used to calculate percent assembled", required=False)
parser.add_argument("-r", "--reverse_reads", help="reverse reads used to generate the assembly, this will be used to calculate percent assembled", required=False)
parser.add_argument("-p", "--percent_aligned", help="Calculate the number of reads that map back to the alignment. This uses bowtie2 and can take some time. It gives an approximation of how representative the contigs are of the raw reads. Of course, it requires that bowtie2 is installed and in your path.", required=False, action='store_true')

args = parser.parse_args()

if args.raw_reads:
    args.percent_aligned = True

#Remove spaces between fasta sequences (to determine contig sequence lengths)
def convert_to_single_fasta():
    command = '''awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'''
    OUTF = open(os.path.abspath("{0}.single".format(args.input_file)), "w")
    with open(args.input_file) as f:
        curr_seq = ""
        for l in f:
            if l[0] == ">":
                if curr_seq != "":
                    OUTF.write(curr_seq+"\n")
                curr_seq = ""
                OUTF.write(l)
            else:
                curr_seq = curr_seq + l.rstrip()
    OUTF.close()
    args.input_file = os.path.abspath("{0}.single".format(args.input_file))

#Check if the input file is in multi or single-lined fasta format
"""
This is a VERY simple check and may give errors. All it does is read in the first 3 lines 
and makes sure the 3rd line starts with a fasta header (e.g. ">"). I can see this easily
going wrong in the future.
"""
def check_fasta_format(fname):
    lines = []
    with open(fname) as f:
        lines.extend(f.readline() for i in xrange(3))
    if lines[2][0] != ">":
        print ""
        print "Detected file is in multi-line fasta format."
        print "Converting file to single-line fasta format."
        convert_to_single_fasta()
        print "Finished conversion, calculating stats."
        print ""

#Count the number of sequences in the file (based on the fasta headers)        
def num_seqs(fname):
    with open(fname) as f:
        count=0
        for l in f:
            if l[0] == ">":
                count = count+1
        return count
        
#calculate N50, average length, number >1kb
def calc_N50(fname):
    #defining variables to be determined
    n50=0
    average_contig_size=0
    contig_over_1kb=0
    #set up list of contig sizes & sort size
    seq_length_array = []
    with open(fname) as f:
        for l in f:
            if l[0] == ">":
                pass
            else:
                seq_len = len(l)-1
                seq_length_array.append(seq_len)
    seq_length_array.sort(key=int)
    largest_contig = seq_length_array[len(seq_length_array)-1]
    total_size = sum(seq_length_array)
    average_contig_size = total_size / len(seq_length_array)
    contig_over_1kb = [n for n in seq_length_array if n > 999]
    num_contig_over_1kb = len(contig_over_1kb)
    #run loop to calculate the n50
    running_total = 0
    n50 = 0
    for i, elem in enumerate(seq_length_array):
        running_total = running_total + elem
        if running_total < (total_size / 2):
            pass
        elif running_total == (total_size / 2):
            print elem
            n50 = ((seq_length_array[i] + seq_length_array[i+1])/2)
        else:
            n50 = elem
            break
    return(n50, largest_contig, average_contig_size, num_contig_over_1kb)

def run_bowtie2_scaffold(contig):
    print "Building the scaffold file for Bowtie2."
    outfile = os.path.abspath(contig)
    outfile = "{0}.bowtie".format(outfile)
    sub.call(["bowtie2-build", args.input_file, outfile], stdout=sub.PIPE, stderr=sub.PIPE)

def run_bowtie2(contig, foward, reverse):
    print "Mapping raw reads using Bowtie2, this may take awhile..."
    outfile = os.path.abspath(contig)
    index_file = "{0}.bowtie".format(outfile)
    outfile = "{0}.sam".format(outfile)
    output = sub.check_output(["bowtie2", "-x", index_file, "-1", foward, "-2", reverse, "-S", outfile], stderr=sub.STDOUT)
    percent_aligned = 0
    for line in output.split(os.linesep):
        if re.search("overall alignment", line):
            line_fields = line.split()
            percent_aligned = line_fields[0]
            return(percent_aligned)

def clean_up():
    #remove the bowtie & the sam files
    filelist = [f for f in os.listdir(".") if f.endswith(".bt2") ]
    for f in filelist:
        os.remove(f)
    os.remove("{0}.sam".format(args.input_file))

    
##Capturing all the outputs
check_fasta_format(args.input_file)
num_seqs = num_seqs(args.input_file)
seq_lens = calc_N50(args.input_file)
if args.percent_aligned:
    run_bowtie2_scaffold(args.input_file)
    percent_assem = run_bowtie2(args.input_file, args.foward_reads, args.reverse_reads)
    clean_up()


row_format = "{:<30}" * 2
print
print row_format.format(*["Number of Contigs", num_seqs])
print row_format.format(*["N50", seq_lens[0]])
print row_format.format(*["Largest Contig", seq_lens[1]])
print row_format.format(*["Average Contig Size", seq_lens[2]])
print row_format.format(*["Number of Contigs >1kb", seq_lens[3]])
if args.percent_aligned:
    print row_format.format(*["Percent of Reads used", percent_assem])

