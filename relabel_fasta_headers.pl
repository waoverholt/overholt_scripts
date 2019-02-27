#!/usr/bin/env perl
use strict;
use warnings;

my $in = $ARGV[0];
my $label = $ARGV[1];

my $i = 0;
open(IN, $in) or die "Can't open input file\n";

while (<IN>) {
    if ($_ =~ m/>/) {
	$i++;
	print ">".$label."_".$i."\n";
    }
    else {
	print $_;
    }
}
