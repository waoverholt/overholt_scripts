#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Text::CSV;

=pod

Written by: Will Overholt
4-16-2013

Script takes the raw .csv file from the bioanalyzer as input
It goes through and converts to bioanalyzer output format to genemapper's input
which is required by T-REX (http://trex.biohpc.org/index.aspx).

Input = bioanalyzer.csv file
output_dir = directory you want to save the file

Script outputs:
[1] arisa_out.txt
        - this is the converted fragment file
        - columns are: (1) dye,peak# [just used blue as default] (2) sample name; 
          (3) peak size; (4) peak height; (5) peak area; (6) default 1 [required by T-REX, not used]

        - to modify the output columns change the array values in line 71

[2] arisa_label.txt
        - T-REX requires a labeling file
        - this file only contains the sample name by default
        - user will have to add meta data (see T-REX documentation)


=cut

#input files
my $file = $ARGV[0];
my $out_dir = $ARGV[1];

my $csv = Text::CSV->new();
open (CSV, "<$file") or die $!;

#files the script creates as output
open (RESULTS, ">$out_dir/arisa_out.txt");
open (LABEL, ">$out_dir/arisa_label.txt");

print LABEL "FileName\n";

my @array;
my $name;
my $i;
my $j;

while(<CSV>) { #read each line from the CSV file
    unless ($_ =~ m/^Name/ || $_ =~ m/^Size/) { #screws up if I don't skip these lines (symbols I think?)
	if ($csv->parse($_)) {
	    my @columns = $csv->fields();
	    if ($columns[0] =~ m/Sample Name/) {
		foreach (@array) {
		    print RESULTS $_."\n";
		}
		@array=();
		$j++;
		$name = $columns[1];
		print LABEL "$name\n";
		$i = 0;
	    }
	    elsif ($columns[0] =~ m/[0-9]+/) {
		$i++;
		my $size = $columns[0];
		$size =~ s/,//g;
		push(@array, "B,$i\t$name\t$size\t$columns[6]\t$columns[9]\t1");
	    }		
	}
	else {
	    my $err = $csv->error_input;
	    print "Failed to parse lin: $err";
	}
    }
}
foreach (@array) {
    print RESULTS $_."\n";
}
close CSV;

