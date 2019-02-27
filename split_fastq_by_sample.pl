#!/usr/bin/env perl
use warnings;
use strict;

my $sample_list = $ARGV[0];
my $seqs = $ARGV[1];
my $out_dir = $ARGV[2];

`rm -r $out_dir`;
`mkdir $out_dir`;

open(NAMES, $sample_list);

my @sample_names = <NAMES>;
close(NAMES);

foreach (@sample_names) {
    chomp;
    unless ($_ =~ m/#/) {
	my $sample_name = $_;
	open(RESULTS, ">>", $out_dir."/".$sample_name.".fastq");
	my $regex = "^\@".$sample_name;
	#my $regex = "^\@Azivy-6-3";
	#print `echo 'grep -A 3 "$regex" $seqs'`;
	my @results =`grep -A 3 $regex $seqs`;
	foreach (@results) {
	    print RESULTS $_;
	}
    }
}

