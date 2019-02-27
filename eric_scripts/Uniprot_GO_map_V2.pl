#!/usr/bin/perl


use warnings;


# Author: Eric Johnston
# Date: 07/16/2014


open(IN1,$ARGV[0]); ###Tab-separated counts (matrix) table
open(IN2,$ARGV[1]); ###Output from 'Uniprot_GO_mapping_F.pl' or 'Uniprot_GO_mapping_P.pl'


my $first_line = <IN1>;
# print "\t".$first_line;

my %hash_ID;

while (<IN1>) {
	chomp;
	my @col = split("\t",$_);
	my ($SP_ID) = $col[0] =~ m/sp\|.+\|(.+)/;
	$col[0] = $SP_ID;
	my $storage_array = join("\t",@col);
	$hash_ID{$SP_ID} = $storage_array;
}
# foreach my $keys (keys %hash_ID) {
# 	print $hash_ID{$keys}."\n";
# 	

close(IN1);

while (<IN2>) {
	chomp;
	my @col = split("\t",$_);
	if (exists $hash_ID{$col[0]}) {
		print $col[1]."\t".$hash_ID{$col[0]}."\n";
	}	
}
	
close(IN2);