#!/usr/bin/perl

use strict;
use warnings;

my $fp = "/data/qiime_files/split_ANL-12-20-12/check-samples";

system("filter_fasta.py -f ~/data/qiime_files/split_ANL-12-20-12/check-sample/seqs_fixed.fna --sample_id_fp ~/data/qiime_files/split_ANL-12-20-12/check-samples/analysis_with_denovo_clust/parallel_ref_otus/Meta_only.txt -o ~/data/qiime_files/split_ANL-12-20-12/check_samples/meta_only_seqs.fna")

system("parallel_pick_otus_uclust_ref.py -i $fp/meta_only_seqs.fna -o $fp/meta_only_pf_script -r $GG_otus -s 0.6 -z -O 8")

system("filter_fasta.py -f $fp/meta_only_seqs.fna -s $fp/meta_only_pf_script/meta_only_seqs_failures.txt -n -o $fp/meta_seqs_pf.fna")
