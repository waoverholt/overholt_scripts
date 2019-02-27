#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::Tools::Run::StandAloneBlastPlus;
use Getopt::Long;
use File::Spec;
use File::Basename;

=pod 

written by Will Overholt
for questions/comments please contact Will at waoverholt@gatech.edu
Spring 2013
adapted from a script written by Dr. Konstantinidis

Purpose of script:
This script is used to calculate AAI (average amino acid identity) or ANI (average nucleotide identity)
from a folder containing predicted gene files
Each predicted gene file (can be translated) should be in a multi-fasta format, and each genome
should have its own multifasta file in the in_dir.
Script will output a list of the comparisons in the terminal as well as write a distance matrix
to the specified output file

this script can be used along with parse_genes_from_draftgenome.pl, translate_fasta_seqs.pl, find_orthologs.pl, and get_seqs_of_ortho_genes.pl

#########################################################
#be sure blast+ is in your path
#########################################################

Options:
-i input directory containing predicted gene sequences for each genome
-p what type of sequences are you comparing? protein (p) or nucleic acid (n)
-o output file path (where do you want to save the generated distance matrix)
-t flag if you want to calculate AAI from DNA sequences
-c number of cores blast will use when determining orthologs

=cut

#genome1 and genome2 labels are somewhat arbitrary since only reciprocal blsat hits are reported
my $in_dir; #make sure the path has a trailing black slash (I'm too lazy to account for this right now)
my $program; #current options are p (protein) or n (nucleotide)
my $output_file;
my $cores=1;
my $translate=0; #you can calculate AAI from a list of predicted DNA sequences, puts translated sequences in directory titled "translated_seqs"

GetOptions(
    "i=s" => \$in_dir,
    "p=s" => \$program,
    "o=s" => \$output_file,
    "t" => \$translate,
    "c=i" => \$cores) or die "missing or incorrect options\n";

#initially uses bioperl to translate list of DNA sequences for each genome
#then will execute blastp to calculate AAI

if ($translate == 1) {
    $program="p";
    my $dir_gene_predictions=$in_dir;
    #print $in_dir."\n";

    my @dirs = File::Spec->splitdir($in_dir);
    @dirs=splice(@dirs, 0, (@dirs - 2));
    my $dir_faa_out = $dirs[0];
    for (my $j = 1; $j < @dirs; $j++) {
	$dir_faa_out=$dir_faa_out."/".$_;
    }

    $dir_faa_out=$dir_faa_out."/translated_seqs/";
    `mkdir $dir_faa_out`;

    opendir(D, $in_dir) or die "Can't open directory specified";
    my @list = grep !/^\.\.?$/, readdir(D);
    closedir(D);
    my @sorted_list = sort(@list);
    
    foreach (@sorted_list) {
	open(OUT, ">", "$dir_faa_out/$_".".faa");
	my $file = $_;
	my $seqio_obj = Bio::SeqIO->new(-file => "$dir_gene_predictions/$file", -format => "fasta");
	while (my $seq_obj = $seqio_obj->next_seq) {
	    my $prot_obj = $seq_obj->translate;
	    print OUT ">".$seq_obj->primary_id."\n";
	    print OUT $prot_obj->seq()."\n";
	}
	
	close(OUT);
    }
    $in_dir=$dir_faa_out;
}


#Read the input files from the specificed directory and defines the distance matrix for the output file
opendir(D, $in_dir) or die "Can't open directory specified, please make sure the parth is correct";

#This reads the list of files in the directory and filters out the default dirs (. and ..)
my @list = grep !/^\.\.?$/, readdir(D);
closedir(D);

#sort the list of files
my @sorted_list = sort(@list);


#setting up the pairwise comparisons, direct comparison will be printed in the terminal, output file is a distance matrix

my $i;
my $dist->[0] = [@sorted_list]; #populate distance matrix with genomes as column labels
unshift(@{$dist->[0]}, 0);

for ($i=0; $i < @sorted_list; $i++) {
    my $filename_ref = "$in_dir".$sorted_list[$i];
    my @temp_array;
    
    for (my $j = 0; $j <= $i; $j++) {
	push(@temp_array, 0); #fill in 0s where you are not directly comparing the 2 genomes
    }
    push(@{$dist->[$i+1]}, $sorted_list[$i]); #add row labels
    push(@{$dist->[$i+1]}, @temp_array); #add correct number of 0s to set up adding ANI values to correct location

    #setting up the pairwise comparions using the input files from the input directory
    for (my $k = $i+1; $k < @sorted_list; $k++) {
	my $filename_db = "$in_dir".$sorted_list[$k];
	
        #create the blast database for the one way comparison
	my $fac = Bio::Tools::Run::StandAloneBlastPlus->new(
	    -db_dir => $in_dir,
	    -db_data => $filename_db,
	    -create => 1,
	    );
	#the database will be tempory and will be deleted at the end of each loop
	my $dbname = $fac->db;
	
	#check to see if you want to compare amino acids or DNA bases
	if ($program eq "p") {
	    my @blast1_results=();
	    my @blast2_results=();
	    #set up blast run
	    my $result = $fac->blastp(
		#-outfile => $in_dir."blast_results",
		-db_data => $dbname,
		-query => $filename_ref,
		-db_make_args => [ '-outfmt' => '6 qseqid sseqid length pid' ],
		-method_args => [ '-num_threads' => $cores,
				  '-evalue' => 0.00001,
				  '-num_descriptions' => 1,
				  '-num_alignments' => 1,
		]
		);

	    #print results from blast
	    while (my $blast_hit = $fac->next_result) {
		while (my $best_hit = $blast_hit->next_hit) {
		    while (my $hsp = $best_hit->next_hsp) {
			if ($hsp->length('total') / $result->query_length > 0.7 && $hsp->percent_identity > 30) {
			    push(@blast1_results, $blast_hit->query_name."\t".$best_hit->name."\t".$hsp->length('total')."\t".$hsp->percent_identity);
			    last;
			    
			}
		    }
		}
	    }
	    $fac->cleanup; #deletes temporary blastdb files
	    
	    #now we set up the reciprocal blast to get only true orthologs (remove possible paralogs)
	    my $fac2 = Bio::Tools::Run::StandAloneBlastPlus->new(
		-db_dir => $in_dir,
		-db_data => $filename_ref,
		-create => 1,
		);
	    
	    my $dbname2 = $fac2->db;
	    my $result2 = $fac2->blastp(
		-db_data => $dbname2,
		-query => $filename_db,
		-db_make_args => [ '-outfmt' => '6 qseqid sseqid length pid' ],
		-method_args => [ '-num_threads' => $cores,
				  '-evalue' => 0.00001,
				  '-num_descriptions' => 1,
				  '-num_alignments' => 1,
		]);
	    #collect blast results from the reciprocal blast
	    while (my $blast_hit2 = $fac2->next_result) {
		while (my $best_hit2 = $blast_hit2->next_hit) {
		    while (my $hsp2 = $best_hit2->next_hsp) {
			if ($hsp2->length('total') / $result2->query_length > 0.7 && $hsp2->percent_identity > 30) {
			    push(@blast2_results, $blast_hit2->query_name."\t".$best_hit2->name."\t".$hsp2->length('total')."\t".$hsp2->percent_identity);
			    last;
			    
			}
		    }
		}
	    }
	    #remove strings introduced by blast
	    foreach (@blast1_results) {
		$_ =~ s/lcl\|//g;
	    }
	    foreach (@blast2_results) {
		$_ =~ s/lcl\|//g;
	    }
	    #sort each results by the same sequence (query of first blast run)
	    @blast1_results = sort { (split '\t', $a)[0] cmp (split '\t', $b)[0] } @blast1_results;
	    @blast2_results = sort { (split '\t', $a)[1] cmp (split '\t', $b)[1] } @blast2_results;
	    
	    my $sum_ANI=0;
	    my $counter=0;
	    my $num_ANI=0;
	    foreach (@blast1_results) {

		my $ref_genome_name = fileparse($filename_ref, qr/\.[^.]*/);
		my $db_genome_name = fileparse($filename_db, qr/\.[^.]*/);

		open(HOMOLOG_OUT, ">", $in_dir."../".$ref_genome_name."_vs_".$db_genome_name.".txt");
		for (my $i = $counter; $i < @blast2_results; $i++) {
		    my @lines1 = split('\t', $_);
		    my @lines2 = split('\t', $blast2_results[$i]);
		    if  ($lines1[0] eq $lines2[1]) { #only calculate ANI from orthologs
			print HOMOLOG_OUT join("\t",@lines1)."\n";
			$sum_ANI += $lines1[3];
			$counter=$i;
			$num_ANI++;
			last;
		    }
		}
		close(HOMOLOG_OUT)
	    }
	    my $ANI = $sum_ANI / $num_ANI;
	    print "$sorted_list[$i]\t$sorted_list[$k]\t$ANI\t$num_ANI\n";
	    
	    $fac2->cleanup;

	    #populate 2-way ANI values into the distance matrix structure
	    push(@{$dist->[$i+1]}, $ANI);

	}

	#crappy, lazy code, but I just duplicated the above section for blastn

	elsif ($program eq "n") {
	    my @blast1_results=();
	    my @blast2_results=();
	    my $result = $fac->blastn(
		#-outfile => $in_dir."blast_results",
		-db_data => $dbname,
		-query => $filename_ref,
		-method_args => [ '-num_threads' => $cores,
				  '-evalue' => 0.00001,
				  '-num_descriptions' => 1,
				  '-num_alignments' => 1,
		]
		);
	    while (my $blast_hit = $fac->next_result) {
		while (my $best_hit = $blast_hit->next_hit) {
		    while (my $hsp = $best_hit->next_hsp) {
			if ($hsp->length('total') / $result->query_length > 0.7 && $hsp->percent_identity > 30) {
			    push(@blast1_results, $blast_hit->query_name."\t".$best_hit->name."\t".$hsp->length('total')."\t".$hsp->percent_identity);
			    last;
			}
		    }
		}
	    }
	    
	
	    $fac->cleanup;
	    
	    my $fac2 = Bio::Tools::Run::StandAloneBlastPlus->new(
		-db_dir => $in_dir,
		-db_data => $filename_ref,
		-create => 1,
		);
	    
	    my $dbname2 = $fac2->db;
	        
	    my $result2 = $fac2->blastn(
		-db_data => $dbname2,
		-query => $filename_db,
		-method_args => [ '-num_threads' => $cores,
				  '-evalue' => 0.00001,
				  '-num_descriptions' => 1,
				  '-num_alignments' => 1,
		]);
	    while (my $blast_hit2 = $fac2->next_result) {
		while (my $best_hit2 = $blast_hit2->next_hit) {
		    while (my $hsp2 = $best_hit2->next_hsp) {
			if ($hsp2->length('total') / $result2->query_length > 0.7 && $hsp2->percent_identity > 30) {
			    push(@blast2_results, $blast_hit2->query_name."\t".$best_hit2->name."\t".$hsp2->length('total')."\t".$hsp2->percent_identity);
			    last;
			}
		    }
		}
	    }

	    foreach (@blast1_results) {
		$_ =~ s/lcl\|//g;
	    }
	    foreach (@blast2_results) {
		$_ =~ s/lcl\|//g;
	    }
	    @blast1_results = sort { (split '\t', $a)[0] cmp (split '\t', $b)[0] } @blast1_results;
	    @blast2_results = sort { (split '\t', $a)[1] cmp (split '\t', $b)[1] } @blast2_results;


	    my $sum_ANI=0;
	    my $counter=0;
	    my $num_ANI=0;
	    foreach (@blast1_results) {
		for (my $i = $counter; $i < @blast2_results; $i++) {
		    my @lines1 = split('\t', $_);
		    my @lines2 = split('\t', $blast2_results[$i]);
		    if  ($lines1[0] eq $lines2[1]) {
			$sum_ANI += $lines1[3];
			$counter=$i;
			$num_ANI++;
			last;
		    }
		}
	    }
	    my $ANI = ($sum_ANI / $num_ANI);
	    print "$sorted_list[$i]\t$sorted_list[$k]\t$ANI\t$num_ANI\n";
	    #print "$filename_ref\t$filename_db\t$ANI\t$num_ANI\n";
	    $fac2->cleanup;
	    push(@{$dist->[$i+1]}, $ANI);
	}
    }
}
pop(@{$dist});

open(OUTPUT, ">", $output_file);

foreach (@{$dist}) {
    foreach (@{$_}) {
	print OUTPUT $_."\t";
    }
    print OUTPUT  "\n";
}
close(OUTPUT);
