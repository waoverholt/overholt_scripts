#!/usr/bin/perl
use warnings;
use strict;

my $otu_to_keep = $ARGV[0];
my $otu_map = $ARGV[1];

open(OTUS, $otu_to_keep);
open(MAP, $otu_map);

my @otus = <OTUS>;
my @map = <MAP>;

close(OTUS);
close(MAP);

foreach (@otus) {
    unless ($_ =~ m/#.*/) {
	#print $_;
	my @matches = grep { /^$_\t.*/ } @map;
	foreach (@matches) {
	    print $_;
	}
    }
}


