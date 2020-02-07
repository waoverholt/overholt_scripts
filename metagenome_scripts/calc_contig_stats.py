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

if args.forward_reads:
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
      
def calc_N50(seq_len_array, total_size):
    running_total = 0
    n50 = 0
    for i, elem in enumerate(seq_len_array):
        running_total = running_total + elem
        if running_total < (total_size / 2):
            pass
        elif running_total == (total_size / 2):
            print elem
            n50 = ((seq_len_array[i] + seq_len_array[i+1])/2)
        else:
            n50 = elem
            break
    return(n50)

#calculate N50, average length, number >1kb
def calc_stats(fname):
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
    num_seqs = len(seq_length_array)
    largest_contig = seq_length_array[len(seq_length_array)-1]
    total_size = sum(seq_length_array)
    average_contig_size = total_size / len(seq_length_array)
    n50 = calc_N50(seq_length_array, total_size)

    #contigs > 1kb
    contig_over_1kb = [n for n in seq_length_array if n > 999]
    num_contig_over_1kb = len(contig_over_1kb)
    len_1kb = sum(contig_over_1kb)
    avg_1kb = len_1kb / num_contig_over_1kb
    n50_1kb = calc_N50(contig_over_1kb, len_1kb)
    #contigs > 3kb
    contig_over_3kb = [n for n in seq_length_array if n > 2999]
    num_contig_over_3kb = len(contig_over_3kb)
    len_3kb = sum(contig_over_3kb)
    avg_3kb = len_3kb / num_contig_over_3kb
    n50_3kb = calc_N50(contig_over_3kb, len_3kb)

    return(
        num_seqs, n50, total_size, average_contig_size, 
        num_contig_over_1kb, n50_1kb, len_1kb, avg_1kb,
        num_contig_over_3kb, n50_3kb, len_3kb, avg_3kb,
        largest_contig,
        )

def run_bowtie2_scaffold(contig):
    print "Building the scaffold file for Bowtie2."
    outfile = os.path.abspath(contig)
    outfile = "{0}.bowtie".format(outfile)
    sub.call(["bowtie2-build", args.input_file, outfile], stdout=sub.PIPE, stderr=sub.PIPE)

def run_bowtie2(contig, forward, reverse):
    print "Mapping raw reads using Bowtie2, this may take awhile..."
    outfile = os.path.abspath(contig)
    index_file = "{0}.bowtie".format(outfile)
    outfile = "{0}.sam".format(outfile)
    output = sub.check_output(["bowtie2", "-x", index_file, "-1", forward, "-2", reverse, "-S", outfile], stderr=sub.STDOUT)
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
seq_lens = calc_stats(args.input_file)
if args.percent_aligned:
    run_bowtie2_scaffold(args.input_file)
    percent_assem = run_bowtie2(args.input_file, args.forward_reads, args.reverse_reads)
    clean_up()


row_format = "{:<30}" * 2
print
print "All Contigs"
print row_format.format(*["Number of Contigs", seq_lens[0]])
print row_format.format(*["N50", seq_lens[1]])
print row_format.format(*["Total Assembly Size", seq_lens[2]])
print row_format.format(*["Average Contig Size", seq_lens[3]])
print
print "Contigs > 1kb"
print row_format.format(*["Number of Contigs >1kb", seq_lens[4]])
print row_format.format(*["N50 of >1kb", seq_lens[5]])
print row_format.format(*["Total Assembly Size > 1kb", seq_lens[6]])
print row_format.format(*["Average Contig Size > 1kb", seq_lens[7]])
print
print "Contigs > 3kb"
print row_format.format(*["Number of Contigs >3kb", seq_lens[8]])
print row_format.format(*["N50 of >3kb", seq_lens[9]])
print row_format.format(*["Total Assembly Size > 3kb", seq_lens[10]])
print row_format.format(*["Average Contig Size > 3kb", seq_lens[11]])
print
print row_format.format(*["Largest Contig", seq_lens[12]])

if args.percent_aligned:
    print row_format.format(*["Percent of Reads used", percent_assem])
