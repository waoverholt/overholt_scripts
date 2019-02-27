#!/usr/bin/perl


use warnings;
use strict;
use List::MoreUtils 'pairwise';


# Author: Eric Johnston
# Date: 07/16/2014


open(IN1,$ARGV[0]); ###Sorted output from 'Uniprot_GO_consolidating.pl'


my @col;
my $old_PF_name;
my $counts;
my @counts_array;
my @existing_array;
my @sum_array;

###########
my $initial_line = <IN1>;
@col = split("\t", $initial_line);
$old_PF_name = $col[0];
($counts) = $initial_line =~ m/.+?\t.+?\t(.+)/;
@counts_array = split("\t", $counts);
@existing_array = @counts_array;
###########



while (<IN1>) {
	chomp;
	@col = split("\t", $_);
	($counts) = $_ =~ m/.+?\t.+?\t(.+)/;
	@counts_array = split("\t", $counts);
	my $new_PF_name = $col[0];
	
	if ($new_PF_name eq $old_PF_name) {
		@sum_array = pairwise { $a + $b } @counts_array, @existing_array;
		@existing_array = @sum_array;
	} else {
		print $old_PF_name."\t".join("\t",@existing_array)."\n";
		@existing_array = @counts_array;
	}		
	$old_PF_name = $new_PF_name;
}	
	
print $old_PF_name."\t".join("\t",@existing_array)."\n";
@existing_array = @counts_array;
	
close(IN1);	
	
	
	
	