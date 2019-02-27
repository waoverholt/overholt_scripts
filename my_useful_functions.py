#!/usr/bin/env python
import sys
import os
from Bio import SeqIO
import subprocess as sub
import shlex
import re


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
        handle = open(os.path.join(outdir, filename), "w")
        count = SeqIO.write(batch, handle, seqformat)
        handle.close()

def run_sys_command(command):
    inargs = shlex.split(command)
    program_call = sub.Popen(inargs)
    program_call.wait()

def label_fasta_as_OTU_X(infile, outfile):
    o = open(outfile, "w")
    i = 0
    with open(infile) as fp:
        for line in fp:
            if line[0] == '>':
                o.write('>OTU_{0}\n'.format(i))
                i += 1
            else:
                o.write(line)
    o.close()

def add_cluster_sizes_to_rep_set(infile, mapping_file, outfile):
    o = open(outfile,"w")
    mapping = open(mapping_file, "r")
    mapping_contents = mapping.readlines()
    mapping_sorted = sorted(mapping_contents)
    #print mapping_sorted
    mapping.close()
    #for element in mapping_sorted:
    #    print element

    i = 0
    with open(infile, "r") as fp:
        for line in fp:
            line = line.rstrip()
            if line[0] == '>':
                seq_name = line.split()
                mapping_id = mapping_sorted[i].split('\t')
                mod_seq_name = re.search(">(.*)", seq_name[0])
                if not mod_seq_name.group(1) == mapping_id[0]:
                    sys.exit("Mapping file and sequences are out of order\n");
                else:
                    size = (len(mapping_id)-1)
                    output_str = "{0}_{1};size={2};".format(seq_name[0], seq_name[1], size)
                    o.write(output_str+'\n')
                    i += 1
            else:
                o.write(line+'\n')
    o.close()
