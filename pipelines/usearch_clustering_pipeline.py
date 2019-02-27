#!/usr/bin/env python
import os
import sys
import argparse
import datetime
import logging
from ruffus import *
import shlex
import re
import subprocess as sub
from my_useful_functions import * #need to have the python script "my_useful_functions.py" in your pythonpath

#Get commandline options
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_seqs", help="fasta file of sequences that have previously been quality filtered", required=True)
parser.add_argument("-o", "--output_dir", help="output directory to store all results from pipeline functions", required=True)
parser.add_argument("-c", "--cores", help="specify the maximum number of cores to use", default=1, type=int, required=False)
parser.add_argument("-n", "--seqs_to_split", help="number specifying how many sequences are added to each file when splitting a larger file. 500000 is a good number", default=500000, type=int, required=False)
parser.add_argument("-r", "--ref_seqs", help="If given, this pipeline will test for chimeras against the reference sequence file given. If not, there will be no reference based chimera detection.", required=False)
parser.add_argument("-t", "--test", help="test the script, output the tasks that are uptodate and those that need to be run", required=False, action="store_true")


args = parser.parse_args()

WORK_DIR = os.path.abspath(args.output_dir)
make_sure_path_exists(WORK_DIR)
input_seqs_basename = os.path.splitext(args.input_seqs)[0]
logging.basicConfig(filename=(os.path.join(WORK_DIR, "usearch_clustering_pipeline.log")), level=logging.INFO,format="")
cur_time = str(datetime.datetime.now())
logging.info("\n#Started Running usearch clustering pipeline\n{0}\n\n".format(cur_time))

#First dereplicate using qiime's prefix-suffix otu picker
@files(os.path.abspath(args.input_seqs), os.path.join(WORK_DIR, '{0}_usearch-labeled.fna'.format(input_seqs_basename)))
def relabel_input_seqs(infile, outfile):
    o = open(outfile, "w+")
    with open(infile) as fp:
        for line in fp:
            line = line.rstrip()
            if line[0] == ">":
                first = re.search(">(.*)_[0-9]+\s", line)
                o.write(line+';'+'barcodelabel='+first.group(1)+';\n')
            else:
                o.write(line+'\n')
    o.close()
    logging.info("#Relabeled the input sequences to work with usearch algorithm. Added ;barcodelabel='sample name'; to each sequence\n\n")

@follows(relabel_input_seqs)
@files(relabel_input_seqs, os.path.join(WORK_DIR, "pre-suf"))
def prefix_suffix_dereplicate(in_file, outdir):
    check_installed("pick_otus.py")
    command = 'pick_otus.py -m prefix_suffix -p 500 -i {0} -o {1}'.format(in_file, outdir)
    
    run_sys_command(command)

    logging.info("#Ran prefix suffix otu picking using qiime. The command was:\n{0}\n\n".format(command))


#pick representative sequences from prefix suffix simplification
@files(prefix_suffix_dereplicate, os.path.join(WORK_DIR, "pre_clustered_seqs.fna"))
def pick_rep_set(in_dir, outfile):
    basename = os.path.splitext(args.input_seqs)[0]
    in_file = os.path.join(WORK_DIR, "pre-suf", basename+"_usearch-labeled_otus.txt")
    command = 'pick_rep_set.py -i {0} -f {1} -o {2}'.format(in_file, os.path.abspath(args.input_seqs), outfile)
    run_sys_command(command)

    logging.info("#Ran qiime's pick rep set on the prefix-suffix otu map\n{0}\n\n".format(command))


#add sizes from the dereplicate step to the header file
@files(pick_rep_set, os.path.join(WORK_DIR, "pre_clustered_seqs_sized.fna"))
def add_cluster_sizes(infile, outfile):
    basename = os.path.splitext(args.input_seqs)[0]
    mapping_file = os.path.join(WORK_DIR, "pre-suf", basename+"_usearch-labeled_otus.txt")

    add_cluster_sizes_to_rep_set(infile, mapping_file, outfile)
    
    logging.info("#Added cluster sizes to each prefix-suffix otu header.\n\n")


#sort by cluster size (use ~/overholt_scripts/sort_by_cluster_abundance.py if usearch runs out of memory)    
@files(add_cluster_sizes, os.path.join(WORK_DIR, "pre_clustered_seqs_sized_sorted.fa"))
def sort_by_sizes(infile, outfile):
    check_installed("usearch1.7.0")
    command = 'usearch1.7.0 -sortbysize {0} -minsize 1 -output {1}'.format(infile, outfile)
    
    run_sys_command(command)
    
    logging.info("#Sorted the pre-suf rep seqs by cluster size (so largest clusters act as the seeds\n{0}\n\n".format(command))

#precluster at 99% using cluster_smallmem from usearch (values were recommended from QIIME pick_otus usearch61)
@files(sort_by_sizes, os.path.join(WORK_DIR, "99_centroids.fa"))
def pre_clust99(infile, outfile):
    log_file = os.path.join(WORK_DIR, "usearch_99_cluster.log")
    command = 'usearch1.7.0 -cluster_smallmem {0} -maxaccepts 1 -maxrejects 32 -id 0.97 -minseqlength 64 -wordlength 8 -centroids {1} -strand both -log {2} -usersort'.format(infile, outfile, log_file)
    
    run_sys_command(command)
    
    logging.info("#Pick otus at 99% similarity using cluster_smallmem\n{0}\n\n".format(command))
              

@files(pre_clust99, os.path.join(WORK_DIR, "usearch_otus_rep_97.fna"))
def usearch_clustering(infile, outfile):
    #includes denove chimera filtering
    log_file = os.path.join(WORK_DIR, "usearch_97_cluster.log")
    command = 'usearch1.7.0 -cluster_otus {0} -otus {1} -log {2} -otuid 0.97'.format(infile, outfile, log_file)
    run_sys_command(command)

    logging.info("#Cluster otus at 97% using usearch7.0\n{0}\n\n".format(command))

@files(usearch_clustering, os.path.join(WORK_DIR, "final_otu_rep_set.fna"))
def relabel_otus(infile, outfile):
    if args.ref_seqs == None:
        print "\nReference based chimera checking skipped because no reference file provided\n"
        logging.info("Reference based chimera checking skipped because no reference file provided\n\n")
        label_fasta_as_OTU_X(infile, outfile)
    else:   
        nochim_fp= os.path.join(WORK_DIR, "usearch_otus_rep_97_nochim.fna")
        refs = os.path.abspath(args.ref_seqs)
        command = 'usearch1.7.0 -uchime_ref {0} -db {1} -nonchimeras {2} -strand plus -threads {3}'.format(infile, refs, nochim_fp, args.cores) 
        run_sys_command(command)

        label_fasta_as_OTU_X(nochim_fp, outfile)

        logging.info("#Removed chimeras using reference based chimera removal against db: {0}\n{1}\n\n".format(refs, command))
        
@follows(relabel_otus)
@split(relabel_input_seqs, '%s/split_input_fasta/*.chunks' % WORK_DIR)
def split_input_fasta(infile, output_files):
    output_dir = os.path.join(WORK_DIR, "split_input_fasta")
    split_large_file(infile, args.seqs_to_split, input_seqs_basename, "fasta", "chunks", output_dir) 

@transform(split_input_fasta, suffix(".chunks"), ".uc")
def map_seqs_to_otus(infile, outfile):
    command = 'usearch1.7.0 -usearch_global {1} -db {0} -strand both -id 0.97 -uc {2}'.format(os.path.join(WORK_DIR, "final_otu_rep_set.fna"), infile, outfile)
    run_sys_command(command)
    
    logging.info("#Split the input sequences into multiple files and then mapped these back to the usearch identified OTUs\n{0}\n\n".format(command))

@merge(map_seqs_to_otus, os.path.join(WORK_DIR, "full_otu_mapping.uc"))
def merge_mapping_results(infiles, summary_file): #script messes up if you are rerunning it multiple times (starts overcounting sequences)
    try:
        with open(summary_file) as summary_file:
            summary_file.close()
            os.remove(summary_file)

        o = open(summary_file, "a")
        for filename in infiles:
            with open(filename, "r") as fp:
                for line in fp:
                    o.write(line)
        o.close()

    except IOError:
        o = open(summary_file, "a")
        for filename in infiles:
            with open(filename, "r") as fp:
                for line in fp:
                    o.write(line)
        o.close()

@files(merge_mapping_results, os.path.join(WORK_DIR, "otu_table.txt"))
def convert_uc2otu(infile, outfile):
    path_to_uc2otutab = "/data/home/woverholt3/bin/uc2otutab.py"
    if os.path.isfile(path_to_uc2otutab) == True:
        command = 'python {0} {1}'.format(path_to_uc2otutab, infile)
        inargs = shlex.split(command)
        program = sub.Popen(inargs, stdout=open(outfile, 'w'))
        program.wait()


        logging.info("#Convert the usearch mapping file into an otu table\n{0}\n\n".format(command))

    elif os.path.isfile(path_to_uc2otutab) == False:
        sys.exit("\n\n{0} does not exist. Please change line 140 in this script (path_to_uc2otutab variable) to reflect the absolute path on where you installed the script uc2otutab.py that can be downloaded from http://drive5.com/python/\n".format(path_to_uc2otutab))


    
@files(convert_uc2otu, os.path.join(WORK_DIR, "otu_table.biom"))
def convert_to_biom(infile, outfile):
    command = 'biom convert -i {0} -o {1} --table-type="otu table"'.format(infile, outfile)
    run_sys_command(command)

    logging.info("#Convert the otu text file into a .biom file\n{0}\n\n".format(command))

if args.test == False:
    #pipeline_run([relabel_input_seqs], multiprocess=6)
    pipeline_run([convert_to_biom], multiprocess=6)
else: 
    pipeline_printout(sys.stdout, [convert_to_biom], verbose = 4)

#pipeline_printout_graph(os.path.join(WORK_DIR, "flowchart.jpg"), "jpg", [convert_to_biom], no_key_legend=True)
