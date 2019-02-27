#!/usr/bin/perl

=pod
Written by Will Overholt
10-23-13

This script takes a GFF file from a genome and searches it for a list of key words
Provide a directory of GFFs and it will sequentially search through them

Currently, you can specify 2 files with a list of REGEX key words to search

Script will output the GFF_file_name followed by the number of hits of each gene.

It was written to find Aliphatic and aromatic degradation genes:
~/Projects/Beach_sands/genomes

=cut


use strict;
use warnings;

my $gff_dir = $ARGV[0] or die "Please enter a path to a gff file.\n";
my $aliphatics = $ARGV[1] or die "Please enter a file with gene names.\n";
my $aromatics = $ARGV[2];

opendir(D, $gff_dir) or die "Can't open directory specified, please make sure the path is correct\n";
my @file_list = grep !/^\.\.?$/, readdir(D);

open(AROMATIC, $aromatics);
my @aromatics = <AROMATIC>;
close(AROMATIC);

open(ALIPHATIC, $aliphatics);
my @aliphatics = <ALIPHATIC>;
close(ALIPHATIC);

open(OUT, ">>", "/data/home/woverholt3/Desktop/test.txt");
    
foreach my $file (@file_list) {
    #print "$gff_dir".$file."\n";
    open(GFF, "$gff_dir".$file);
    my $count_aliphatic=0;
    my $count_aromatic=0;
    while (<GFF>) {
	foreach my $aliphatic (@aliphatics) {
	    chomp $aliphatic;
	    if ($_ =~ m/$aliphatic/) {
		my @fields = split('\t', $_);
		print OUT $file."\t"."aliphatic"."\t".$fields[3]."\t".$fields[8];
		$count_aliphatic++;
	    }
	}
	foreach my $aromatic (@aromatics) {
	    chomp $aromatic;
	    if ($_ =~ m/$aromatic/) {
		my @fields = split('\t', $_);
		print OUT $file."\t"."aromatc"."\t".$fields[3]."\t".$fields[8];
		$count_aromatic++;
	    }
	}
    }
    print $file."\t".$count_aliphatic."\t".$count_aromatic."\n";
    close(GFF);
}
close(OUT);

