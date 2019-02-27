#!/usr/bin/perl

use strict;
use warnings;

# Author: Eric Johnston
# Date: 10/21/2014


open(IN,$ARGV[0]);

my $UniprotID;
my $GO_annotation;


while(<IN>){
	chomp;
	my @col=split("\t", $_);
	if ($_ =~ m/^ID/) { 
		($UniprotID) = $_ =~ m/^ID\ \ \ ([A-Z;_;0-9;a-z]+)\ ?/;
		next;
	}
	if ($_ =~ m/^DR\ \ \ GO;\ GO:[0-9]+;\ P/) { 
		print $UniprotID."\t".."\n";
		($GO_annotation) = $_ =~ m/^DR\ \ \ GO;\ GO:[0-9]+;\ P:(.+)?;/;
		print $UniprotID."\t".$GO_annotation."\n";
		next;
	}			
}

close(IN);

