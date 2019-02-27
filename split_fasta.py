#!/usr/bin/env python

"""
Written by Will A. Overholt
10/24/14

"""

from my_useful_functions import batch_iterator,split_large_file
import sys, os
import shlex
import subprocess as sub
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_seqs", help="parent sequence file to split", required=True)
parser.add_argument("-o", "--output_dir", help="the directory all split files  will be stored", required=True)
parser.add_argument("-n", "--seqs_to_split", help="number specifying how many sequences to add to each file when splitting larger files (needs to be in a multiple of 4. 1200000 is a good number (400000 sequences)", type=int, required=False)
parser.add_argument("-f", "--number_of_files", help="number of file you want to be created", required=False, type=int)

args = parser.parse_args()

if args.number_of_files is not None and args.seqs_to_split is not None:
    sys.exit("\n\nYou need to specify either -n or -f, not both\n\n")
elif args.number_of_files is None and args.seqs_to_split is None:
    sys.exit("\n\nYou need to define either the number of sequences you want per file, or total number of files\n\n")

if args.seqs_to_split is not None:
    if (args.seqs_to_split % 2 != 0):
        sys.exit("\n\nYou need to specify a value that is a multiple of 2 to prevent breaking up sequences\n\n")


in_base = os.path.splitext(args.input_seqs)[0]

if not os.path.exists(args.output_dir):
    os.mkdir(args.output_dir)

if args.number_of_files is not None:
    num_lines = sum(1 for line in open(args.input_seqs))
    batch_size = num_lines / args.number_of_files
    if batch_size % 2 == 1:
        batch_size = int(round(batch_size + 1, 0))
        command = 'split -l {0} -d {1} {2}/seqs_'.format(batch_size, args.input_seqs, args.output_dir)
        in_args = shlex.split(command)
        split = sub.Popen(in_args, shell=False)
        split.wait()
    elif batch_size % 2 == 0:
        command = 'split -l {0} -d {1} {2}/seqs_'.format(batch_size, args.input_seqs, args.output_dir)
        in_args = shlex.split(command)
        split = sub.Popen(in_args, shell=False)
        split.wait()


elif args.seqs_to_split is not None:
    command = 'split -l {0} -d {1} {2}/seqs_'.format(args.seqs_to_split * 2, args.input_seqs, args.output_dir)
    in_args = shlex.split(command)
    split = sub.Popen(in_args, shell=False)
    split.wait()
