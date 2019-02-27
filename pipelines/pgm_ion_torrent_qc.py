#!/usr/bin/env python

import os,sys,argparse,re
from ruffus import *
import subprocess as sub
import shlex
import shutil
from my_useful_functions import *
import tempfile
from filter_fasta_by_length import *

"""
Remove all spaces in names before you begin.
ls *.fastq | perl -ne 'use File::Copy;chomp;$old=$_;s/\s+/_/g;move($old,$_);'
"""
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_seqs_dir", help="the directory containing the ion torrent .fastq files from the sequencing facility", required=True)
#parser.add_argument("-m", "--sample_name_file", help="a file with each sample listed on its own line", required=True)
parser.add_argument("-o", "--output_dir", help="a file name that will contain the quality controlled sequences", required=True)
parser.add_argument("-p", "--phred", help="the phred score you wish to set as a minimum phred score to keep", required=False, default=14, type=int)
args = parser.parse_args()

###############

IN_DIR = os.path.abspath(args.input_seqs_dir)
tempfile.tempdir = IN_DIR
concat_seqs = os.path.join(os.path.abspath(args.output_dir),'concat_seqs.fna')
concat_hist = os.path.join(os.path.abspath(args.output_dir),'concat_histograms.fna')
concat_log = os.path.join(os.path.abspath(args.output_dir),'concat_logs.txt') 

make_sure_path_exists(os.path.abspath(args.output_dir))


for filename in os.listdir(IN_DIR):
    if re.search("fastq", filename):
        
        temp_dir = tempfile.mkdtemp()
        #print temp_dir
        sample_name = os.path.splitext(os.path.basename(filename))[0]
        temp_mapping = tempfile.NamedTemporaryFile(mode='w+t')
        temp_name = temp_mapping.name
        temp_mapping.write("#SampleID	BarcodeSequence	LinkerPrimerSequence	Description\n")
        temp_mapping.write(sample_name+"\n")
        temp_mapping.seek(0)
        command = 'split_libraries_fastq.py -i {0} -m {1} -o {2} --barcode_type not-barcoded --phred_offset 33 --sample_ids {3}'.format(os.path.join(IN_DIR,filename),temp_name,temp_dir,sample_name)
        in_args = shlex.split(command)
        #print in_args
        qiime = sub.Popen(in_args, shell=False)
        qiime.wait()
        #print filename
        #print sample_name
        for f in os.listdir(temp_dir):
            f = os.path.join(temp_dir,f)
            if re.search("fna", f):
                with open(concat_seqs, 'a+') as outfile:
                    with open(f) as infile:
                        for line in infile:
                            outfile.write(line)
            elif re.search("histogram", f):
                with open(concat_hist, 'a+') as outfile:
                    with open(f) as infile:
                        for line in infile:
                            outfile.write(line)
            elif re.search("log", f):
                with open(concat_log, 'a+') as outfile:
                    with open(f) as infile:
                        for line in infile:
                            outfile.write(line)

            temp_mapping.close()  
        shutil.rmtree(temp_dir)

final_out_fp = os.path.join(os.path.abspath(args.output_dir), 'filtered_concat_seqs.fna')
filter_seq_by_length(concat_seqs, 150, 1000, final_out_fp)
