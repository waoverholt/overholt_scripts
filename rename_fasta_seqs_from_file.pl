#!/usr/bin/perl

use warnings;
use strict;

use Bio::SeqIO;

my $FASTA = $ARGV[0];
my $re_name = $ARGV[1];

my $seqio_obj = Bio::SeqIO->new(-file => $FASTA, -format => 'fasta');

open(RE_NAME, "$re_name");

my @a = <RE_NAME>;

while (my $seq_obj= $seqio_obj->next_seq) {
    my $seq_name = $seq_obj->primary_id();
    my $seq = $seq_obj->seq();
    unless ($seq_name =~ m/_<|>[0-9]+/) {
	foreach (@a) {
	    chomp;
	    my @fields = split("\t", $_);
	    if ($seq_name =~ m/$fields[0]/) {
		print ">".$fields[1]."\n";
		print $seq."\n";
	    }
	}
    }
}


close(RE_NAME);
