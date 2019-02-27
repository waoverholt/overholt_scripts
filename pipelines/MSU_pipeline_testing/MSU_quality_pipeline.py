#!/usr/bin/env python
import re
import sys
import argparse
import os
from ruffus import *
import subprocess as sub
import shlex
import shutil
from Bio import SeqIO
import textwrap
#------------------------------------------------------#

#Get commandline options
my_usage = """
Written by Will Overholt (waoverholt@gatech.edu)
14 Jan 2014

This pipeline was written to facillitate format conversion and initial
quality filtering(using the superior usearch version 7.0 quality filter
settings). It takes raw sequence data returned from the MSU sequencing 
facility. This format has 2 files (forward, reverse reads) for each 
sample, with the barcode in the header of each read
(delimited with colons).

Dependencies:
(1) You must have the QIIME core installation in your path.
(2) You must have usearch version 7.0 in your path. Note it 
    should be called usearch1.7.0 (email me to bitch about my lazyness)
(3) You need to have both ruffus and Biopython installed and 
    in your pythonpath.

This script does NOT check for chimeras. That step will be performed in 
the second pipeline (for OTU clustering using usearch)
"""


parser = argparse.ArgumentParser(usage = my_usage)
parser.add_argument("-i", "--input_seqs_dir", help="the directory storing the raw sequence files from MSU", required=True)
parser.add_argument("-o", "--output_dir", help="the directory all results will be stored", required=True)
parser.add_argument("-m", "--mapping_file", help="mapping file linking barcodes to the sample names", required=True)
parser.add_argument("-n", "--seqs_to_split", help="number specifying how many sequences to add to each file when splitting larger files. 500000 is a good number" , default=500000, type=int, required=False)
parser.add_argument("-c", "--cores", help="maximum number of cores to use when running tasks in parallel", default=1, type=int, required=False)
args = parser.parse_args()

#-------------------------------------------------------#
IN_DIR = os.path.abspath(args.input_seqs_dir)
WORK_DIR = os.path.abspath(args.output_dir)
#-------------------------------------------------------#


#-------------------------------------------------------#
#functions used in the pipelines
def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise

def read_and_write_lines(output_file, input_dir, filename):
    f = open(os.path.join(input_dir, filename), 'r')
    lines = f.readlines()
    for line in lines:
        output_file.write(line)
    f.close

def check_installed(program):
    try:
        devnull = open(os.devnull)
        sub.Popen([program, '--version'], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:

        if e.errno == os.errno.ENOENT:
            return sys.exit("\n\nThe pre-requisite "+ program + " is not in your path. Is it installed?\n")
    return True

def batch_iterator(iterator, batch_size) :
    """Returns lists of length batch_size.
 
    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.
 
    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch

def split_large_file(infile, size, prefix, seqformat, suffix, outdir):
    make_sure_path_exists(outdir)
    record_iter = SeqIO.parse(open(infile), seqformat)
    for i, batch in enumerate(batch_iterator(record_iter, size)):
        filename = "%s_%i.%s.%s" % (prefix, i+1, seqformat, suffix)
        handle = open(os.path.join(WORK_DIR, outdir, filename), "w")
        count = SeqIO.write(batch, handle, seqformat)
        handle.close()        

#-------------------------------------------------------#


#-------------------------------------------------------#
#THE FIRST STEP IN THE PIPELINE
#first need to create joined read files
@files(IN_DIR, [os.path.join(WORK_DIR, 'concat_R1_seqs.fastq'), os.path.join(WORK_DIR, 'concat_R2_seqs.fastq')])
def convert_MSU_seqs_to_qiime(outfile, blank): #needed the blank input because I didn't like the double outfile list
    make_sure_path_exists(WORK_DIR)
    concat_R1 = open(os.path.join(WORK_DIR, 'concat_R1_seqs.fastq'), "a")
    concat_R2 = open(os.path.join(WORK_DIR, 'concat_R2_seqs.fastq'), "a")
    for filename in os.listdir(IN_DIR):
        if re.search("_R1_", filename):
            read_and_write_lines(concat_R1, IN_DIR, filename)
            read_and_write_lines(concat_R2, IN_DIR, filename)
    concat_R1.close()
    concat_R2.close()

#-------------------------------------------------------#
#Make the barcodes file
@files(convert_MSU_seqs_to_qiime, os.path.join(WORK_DIR, 'barcodes.fastq'))
def make_barcode_files(infile, outfile):
    outdir = os.path.join(WORK_DIR, 'temp')
    make_sure_path_exists(outdir)
    check_installed("extract_barcodes.py")
    command = 'extract_barcodes.py -f {0} -c barcode_in_label --char_delineator ":" --bc1_len 12 -o {1}'.format(infile[0], outdir)
    in_args = shlex.split(command)
    qiime = sub.Popen(in_args, shell=False)
    qiime.wait()
    os.rename(os.path.join(WORK_DIR, 'temp', 'barcodes.fastq'), os.path.join(WORK_DIR, 'barcodes.fastq'))
    shutil.rmtree(os.path.join(WORK_DIR, 'temp'))

#-------------------------------------------------------#
#join paired end reads
@files(make_barcode_files, os.path.join(WORK_DIR, 'join_paired_ends'))
def join_paired_ends(infile, outfile):
    read1 = os.path.join(WORK_DIR, 'concat_R1_seqs.fastq')
    read2 = os.path.join(WORK_DIR, 'concat_R2_seqs.fastq')
    barcodes = os.path.join(WORK_DIR, 'barcodes.fastq')
    output = os.path.join(WORK_DIR, 'join_paired_ends')
    command = 'join_paired_ends.py -f {0} -r {1} -b {2} -o {3}'.format(read1, read2, barcodes, output)
    in_args = shlex.split(command)

    qiime = sub.Popen(in_args, shell=False)
    qiime.wait()

#-------------------------------------------------------#
#demultiplex with QIIME - no error checking
@files(join_paired_ends, os.path.join(WORK_DIR, 'qiime_deplex'))
def deplex_qiime(in_dir, out_dir):
    mapping = os.path.abspath(args.mapping_file)
    barcodes = os.path.join(WORK_DIR, 'join_paired_ends', 'fastqjoin.join_barcodes.fastq')
    inseqs = os.path.join(WORK_DIR, 'join_paired_ends', 'fastqjoin.join.fastq')
    output = os.path.join(WORK_DIR, 'qiime_deplex')
    command = 'split_libraries_fastq.py -q 0 --store_demultiplexed_fastq -m {0} --barcode_type golay_12 -b {1} -i {2} -o {3}'.format(mapping, barcodes, inseqs, output)
    in_args = shlex.split(command)
    
    qiime = sub.Popen(in_args, shell=False)
    qiime.wait()

#-------------------------------------------------------#
#quality filter with usearch7.0
@split(deplex_qiime, '%s/qiime_deplex/split_fastq_seqs/*.chunks' % WORK_DIR)
def split_fastq(input_file, output_files):
    check_installed('usearch1.7.0')
    split_large_file(os.path.join(input_file, "seqs.fastq"), args.seqs_to_split, "seqs", "fastq", "chunks", os.path.join(WORK_DIR, "qiime_deplex", "split_fastq_seqs"))
    
@transform(split_fastq, suffix(".chunks"), ".fasta")
def run_usearch_qual(infile, outfile):
    check_installed('usearch1.7.0')
    command = 'usearch1.7.0 -fastq_filter {0} -fastq_maxee 0.5 -fastaout {1}'.format(infile, outfile)
    inargs = shlex.split(command)
    program = sub.Popen(inargs)
    program.wait()

@merge(run_usearch_qual, 'qc_seqs.fasta')
def merge_usearch_qual(infiles, summary_file):
    output_file = open(os.path.join(WORK_DIR, "qiime_deplex", summary_file), "a")
    for filename in infiles:
        read_and_write_lines(output_file, infiles, filename)

#Just make sure entire pipeline is run - rerunning will show you everything is up to date
@posttask(touch_file(os.path.join(WORK_DIR, "MSU_pipeline_now_complete")))
@follows(merge_usearch_qual)
@files(None, "MSU_pipeline_now_complete")
def pipeline_finished(input_file, output_file):
    pass

#pipeline_printout_graph(os.path.join(WORK_DIR, "flowchart.jpg"), "jpg", [pipeline_finished], no_key_legend=True)
pipeline_run([pipeline_finished])

#os.remove(os.path.join(WORK_DIR, "MSU_pipeline_now_complete"))

pipeline_run()
