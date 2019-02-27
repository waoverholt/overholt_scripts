#!/usr/bin/perl

=pod
written by Will Overholt
10-21-2013

This script is meant to count the number of sequence reads that are related to a cultured isolate.
As input, 

(1) run blastn: "blastn -query ~/Projects/sanger_sequencing_raw/101813_Xiaoxu-isolate
s/all_isolates.fasta -db blast_db_rep_set/rep-set-db -outfmt 7 -num_threads 6 > blast_results.txt"

(2) Filter this output to only include matches > some cutoff (94, or 97%): perl -ane 'if ($_ =~ m/#/) {print}; if ($F[2] > 94) {print $_}' blast_results.txt > blast_results_filtered.txt

(3) Run this script passing the otu_mapping file from qiime (final_otu_map_mc2)

It'll count the number of raw sequences that are within the similarity cut-off to your isolate

=cut
use warnings;
use strict;

#set up the user input, mapping file should be first, then the blast results (outfmt 7)
my $otu_map = $ARGV[0];
my $blast_results = $ARGV[1];

open(MAP, $otu_map);

my %otu_map;

#read the mapping file into a hash of arrays
while (<MAP>) {
    my @lines = split('\t', $_);
    $otu_map{$lines[0]} = [ @lines[1..$#lines] ]; #skip the first element
}

close(MAP);

open(SEQ_IDS, $blast_results);


my %isolates;

#read the blast results into a hash of arrays
#each isolate sequence will be its own key, all the results are part of the key
while (<SEQ_IDS>) {
    unless ($_ =~ m/#/) {
	my @fields = split('\t', $_);
	if (exists $isolates{$fields[0]}) {
	    push (@{ $isolates{$fields[0]} }, $fields[1]);
	}
	else {
	    push (@{$isolates{$fields[0]}}, $fields[1]);
	}
    }
}

#search your mapping file for the raw sequences (instead of OTUs)

foreach my $group (keys %isolates) {
    #print "$group: @{$isolates{$group}} \n";
    my $count = 0;
    foreach my $lines (@{ $isolates{$group} }) {
	my $otu_to_lookup = $lines;
	if (exists $otu_map{$otu_to_lookup}) {
	    my $num_elements = @{$otu_map{$otu_to_lookup}};
	    #print "\t".$lines."\t".$num_elements."\n";
	    $count += $num_elements;
	}
    }
    print "$group\t$count\n";
}
 
close(SEQ_IDS);

