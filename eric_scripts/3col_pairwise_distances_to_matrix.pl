#!/usr/bin/perl


use warnings;
use strict;


# Author: Eric Johnston
# Date: 10/21/2016


open(IN1,$ARGV[0]); # pairwise distances as 3-columns


my %sample_ID_hash;
my %pairwise_distance_hash;
while(<IN1>){
	chomp;
	my @cols = split("\t",$_);
	$sample_ID_hash{$cols[0]}=1;
	$sample_ID_hash{$cols[1]}=1;
	$pairwise_distance_hash{$cols[0]."__".$cols[1]}=$cols[2];
}
close(IN1);
for my $keys0 (sort keys %sample_ID_hash) {
	print "\t".$keys0;
}
print "\n";

for my $keys1 (sort keys %sample_ID_hash) {
	print $keys1;
	for my $keys2 (sort keys %sample_ID_hash) {
		print "\t".$pairwise_distance_hash{$keys1."__".$keys2};
	}
	print "\n";
}


