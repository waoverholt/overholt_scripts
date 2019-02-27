#!/usr/bin/perl

=pod

written by Will Overholt
4-3-2013

-----------------------------------------------------
PURPOSE OF SCRIPT
-----------------
Sometimes Stefan doesn't send us trimmed sequences
This script takes the .seq, .qual, and .TXT files
that are generated from phred and from ABI and uses
them to trim out the crap.

It takes as input the directory containing all the
individual sequences you are interested in trimming. It
also takes a description as to where and what the resulting
output files will be called
(e.g. ./this_scrip /path/to/input/dir /path/to/output/output_identifier_string)

When this script was written each individual sequence
had 3 files associated with it. Therefore the first step
was to make three concatenated multifasta files.

Then I call Lucy (http://lucy.sourceforge.net/) and use it
to identify bad regions.

Finally I trim the shit and output the final file (XXX_trimmed.fasta).

=cut

use warnings;
use strict;
use Bio::SeqIO; #you'll need to have the bioperl core modules installed


my $dir = $ARGV[0]; #directory containing all the sequences you want to trim
my $out = $ARGV[1]; #path and unique identifier you want to use for the sequences

opendir(DIR, $dir); 
my @files = readdir(DIR);
closedir(DIR);

my @good_files;
my @fasta;
my @qual;
my @second;

#remove Linux's default directors (. and ..)
foreach (@files) { 
    if ($_ =~ m/[A-z]/) {
	push(@good_files, $_);
    }
}

#set up 3 arrays for each file type
foreach (@good_files) {
    if ($_ =~ m/seq/) {
	push(@fasta, $_);
    }
    elsif ($_ =~ m/qual/) {
	push(@qual, $_);
    }
    elsif ($_ =~ m/TXT/) {
	push(@second, $_);
    }
}
@fasta = sort(@fasta);
@qual = sort(@qual);
@second = sort(@second);

open(FASTA, ">>", $out.".fasta");
open(QUAL, ">>", $out.".qual");
open(SECOND, ">>", $out.".TXT");

#for whatever reason the fasta sequences didn't have a header
#We first add the sequence name to the output file, then
#add the sequence
foreach (@fasta) {
    my @fields = split('\.', $_);
    my $seq_name = $fields[0];
    print FASTA ">".$seq_name."\n";
    open(FILE, $dir.$_);
    my @seq = <FILE>;
    print FASTA @seq;
    close(FILE);
}
#write the qual files to a new file
foreach (@qual) {
    open(FILE2, $dir.$_);
    my @quality = <FILE2>;
    print QUAL @quality;
    close(FILE2);
}
#write the TXT files to a new file
foreach (@second) {
    open(FILE3, $dir.$_);
    my @seq2 = <FILE3>;
    print SECOND @seq2;
    close(FILE3);
}

close(FASTA);
close(QUAL);
close(SECOND);

#run lucy using default parameters
`lucy $out.fasta $out.qual $out.TXT -output $dir"lucy_out.fasta" $dir"lucy_out.qual"`;

#remove the damn spaces lucy adds to the names so Bio::SeqIO can parse the sequence headers
`sed -i 's/ /_/g' $dir"lucy_out.fasta"`;

open(FINAL, ">", $out."_trimmed.fasta");

#trim the sequences, since lucy is a little bitch and doesn't trim them herself
my $new_fasta = $dir."lucy_out.fasta";
my $seqio_obj = Bio::SeqIO->new( -file => $new_fasta, -format => "fasta");
while (my $seq_obj = $seqio_obj->next_seq) {
    my $seq_name = $seq_obj->primary_id;
    my @fields = split("_", $seq_name);
    my $start = $fields[4];
    my $end = $fields[5];
    print FINAL ">".$fields[0]."\n";
    print FINAL $seq_obj->subseq($start, $end)."\n";
}

