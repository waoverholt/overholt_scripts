#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;

my $seqs = $ARGV[0];
my $nodes = $ARGV[1];
my $out_dir = $ARGV[2];
my $basename=basename($seqs);
(my $file_prefix = $basename) =~ s/[\.].*//g;

unless (-d "$out_dir") {
    `mkdir $out_dir`;
}    
else {
    print "Warning, directory already exists.\nDo you want to continue(Y/N)\n";
    my $usr_ans = <STDIN>;
    if ($usr_ans =~ m/(N|n)/) {
	exit;
    }
}

my $seq_count = qx(wc -l $seqs);
$seq_count =~ s/(^[0-9]+).*/$1/;

#print $seq_count;

my $i = 1;
my $j = 0;
my $count = 1;
#print $j+($seq_count/$nodes)."\n";

open(SEQS, $seqs);
open(OUT, ">", "$out_dir\/".$file_prefix."_$count"."txt");

while (<SEQS>) {
    if ($i < ($j+$seq_count/$nodes)) {
	print OUT $_;
	$i++;
    }
    else {
	unless ($_ =~ m/\>/) {
	    print OUT $_;
	    $j = $i;
	    $count++;
	    $i++;
	    close(OUT);
	    open(OUT, ">", "$out_dir\/".$file_prefix."_$count".".txt");
	}
	else {
	    print OUT $_;
	    $i++;
	}
    }
}
