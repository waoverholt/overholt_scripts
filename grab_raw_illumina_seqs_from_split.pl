#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Getopt::Long;
use File::Temp;

my $split_lib_seqs;
my $raw_seqs;
my $out_file;

GetOptions(
    "i=s" => \$raw_seqs,
    "r=s" => \$split_lib_seqs,
    "o=s" => \$out_file) or die "missing or incorrect options\n";

my $seqio_obj = Bio::SeqIO->new(-file => $split_lib_seqs, -format => "fasta");


my $tmp = File::Temp->new(
    TEMPLATE => "tempXXXX",
    DIR => "/data/home/woverholt3/Desktop/",
    suffix => ".fna",
    UNLINK => 1,
    ) or die "Cannot create temp file\n";
my $filename = $tmp->filename;

#my $filename = "/data/home/woverholt3/Desktop/tmp_file";

open my $temp, ">>", $filename or die "Could not open $filename: $!";

while (my $seq = $seqio_obj->next_seq) {
    my $full_name = $seq->desc;
    (my $id = $full_name) =~ s/ orig_bc.*//g;
    print { $temp } "$id\n" or die "Cannot write to $filename\n";
}

`filter_fasta.py -f $raw_seqs -o $out_file -s $filename`
