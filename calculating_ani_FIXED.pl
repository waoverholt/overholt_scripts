#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Bio::SeqIO;

=pod

Written by Will Overholt
For questions or comments please contact me at waoverholt@gatech.edu
Spring 2014
adapted from previous scripts, "find_orthologs.pl" and "calculating_ani_v2.pl"

Purpose of Script:
This script is used to calculate either ANI (average nucleotide identify) or AAI (average amino acid identify from a folder containing predicted genes sequences or predicted protein sequences.

***************************************************
For AAI all files should have a .faa extension!!!!
For ANI all files should have a .ffn extension!!!!
***************************************************

Each predicted gene file should be in multi-line fasta format, and each genome should have its own file in the input directory (in_dir).
Script will output a list of the pairwise comparisons (_3col.txt) and a distance matrix (_dm.txt) into the input directory.

It was written to work with blast 2.2.27+

Options:
-i input directory containing predicted gene sequences for each genome
-p what type of sequences are you comparing? protein (p) or nucleic acid (n)
-o output file name
-c number of cores blast will use when determining orthologs

=cut


#first read directory of .faa files
#define number of cores for blastp to use (default = 1);
my $cores=1;
my $in_dir; #directory of sequences to compare
my $output_file; #output file name (do not have an extension)

GetOptions(
    "i=s" => \$in_dir,
    "o=s" => \$output_file,
    "c=i" => \$cores) or die "missing or incorrect options specified\n";

$in_dir =~ s/\/*$//;
$output_file =~ s/\.[^.]+$//;


my $blast_version = `blastp -version | head -1`;
$blast_version =~ s/\x0a//g;

unless ($blast_version =~ m/blastp: 2.2.27+/) {
    print "You do not have blast installed in your path\nor you are using a different version from 2.2.27+\nIf you wish to proceed, please re-run after commenting out lines 14-17 in the perl script";
    exit;
}


opendir(D, $in_dir) or die "Can't open directory specified";
my @list = grep !/^\.\.?$/, readdir(D);
closedir(D);

#################################################

#translate the gene sequences to determine true orthologs
`mkdir $in_dir/faa_out`;
my @ffn_list;
foreach (@list) {
    if ($_ =~ m/\.*\.ffn/) {
	push(@ffn_list, $_);
    }
}
my @sorted_list = sort(@ffn_list);
foreach (@sorted_list) {
    open(OUT, ">", "$in_dir/faa_out/$_".".faa");
    my $file = $_;
    my $seqio_obj = Bio::SeqIO->new(-file => "$in_dir/$file", -format => "fasta");
    while (my $seq_obj = $seqio_obj->next_seq) {
	my $prot_obj = $seq_obj->translate;
	print OUT ">".$seq_obj->primary_id."\n";
	print OUT $prot_obj->seq()."\n";
    }
    close(OUT)
}

my $in_dir_faa="$in_dir/faa_out";

######################################
#Run blastp to find orthologs from predicted sequences

opendir(D, $in_dir_faa) or die "Can't open directory specified";
my @list1 = grep !/^\.\.?$/, readdir(D);
closedir(D);


my @aa_list;
foreach (@list1) {
    if ($_ =~ m/.*\.faa/) {
	push(@aa_list, $_);
    }
}

my @sorted_list1 = sort(@aa_list);


my $i;
my $dist->[0] = [@sorted_list1];
unshift(@{$dist->[0]}, 0);
my @final_results;

for ($i = 0; $i < @sorted_list1; $i++) {
    my $filename_ref = "$in_dir_faa/".$sorted_list1[$i];
    my $filename_ref_short = $sorted_list1[$i];
    my @temp_array;
    
    for (my $j = 0; $j <= $i; $j++) {
	push(@temp_array, 0);
    }
    push(@{$dist->[$i+1]}, $sorted_list1[$i]); #add row labels
    push(@{$dist->[$i+1]}, @temp_array); #add correct number of zeros to set up adding ANI values to correct location
    
    for (my $k = $i+1; $k < @sorted_list1; $k++) {
	my $filename_db = "$in_dir_faa/".$sorted_list1[$k];
	my $filename_db_short = $sorted_list1[$k];
	
	#set up reference database
	system(
	    "if [ ! -d $in_dir_faa/temp_blast_dbs ];
        then
            mkdir $in_dir_faa/temp_blast_dbs
        fi");
	
	system("makeblastdb -in $filename_db -input_type 'fasta' -dbtype 'prot' -out $in_dir_faa/temp_blast_dbs/$filename_db_short");
	
	system(
	    "if [ ! -d $in_dir_faa/temp_blast_results ]; 
        then
            mkdir $in_dir_faa/temp_blast_results
        fi");
	
	#run first blast
	system("blastp -query $filename_ref -db $in_dir_faa/temp_blast_dbs/$filename_db_short -out $in_dir_faa/temp_blast_results/$filename_ref_short -max_target_seqs 1 -outfmt '6 qseqid sseqid length slen pident' -num_threads $cores");
	#set up reciprocal blast
	system("makeblastdb -in $filename_ref -input_type 'fasta' -dbtype 'prot' -out $in_dir_faa/temp_blast_dbs/$filename_ref_short");
	
	system("blastp -query $filename_db -db $in_dir_faa/temp_blast_dbs/$filename_ref_short -out $in_dir_faa/temp_blast_results/$filename_db_short -max_target_seqs 1 -outfmt '6 sseqid qseqid pident' -num_threads $cores");
	print "\n\n";

	#keep only reciprocal hits
	open(B_RESULTS, "$in_dir_faa/temp_blast_results/$filename_ref_short");
	open(KEEP, ">", "$in_dir_faa/temp_blast_results/$filename_ref_short"."_vs_"."$filename_db_short.good_results");
	open(B_RECIP, "$in_dir_faa/temp_blast_results/$filename_db_short");
	
	my @good_results;
	while(<B_RESULTS>) {
	    my @fields = split('\t', $_);
	    if ($fields[2] / $fields[3] > 0.2 && $fields[4] > 0.2) { #### adjust sensitivity here
		push(@good_results, $_);
		print KEEP $_;
	    }
	}
	
	my @recip_results = <B_RECIP>;
	@good_results = sort @good_results;
	@recip_results = sort @recip_results;
	
	my $count //=0;
	my $number_hits = 0;
	my $sum_ANI = 0;
	my @recip_good_blasts;
	foreach (@good_results) {
	    chomp;
	    #print $_."\n";
	    my @fields = split('\t', $_);
	    my $cur_file_gene_name = $fields[0];
	    my $ref_gene_name = $fields[1];
	    my $percent_id = $fields[4];
	    
	    for (my $iii = $count; $iii < @recip_results;) {
		my @fields_2 = split('\t', $recip_results[$iii]);
		my $cur_gene_name2 = $fields_2[0];
		my $ref_gene_name2 = $fields_2[1];
		if ($cur_file_gene_name eq $cur_gene_name2 && $ref_gene_name eq $ref_gene_name2) {
		    #print $cur_file_gene_name."\t".$ref_gene_name."\t".$percent_id."\n";
		    $sum_ANI += $percent_id;
		    $number_hits++;
		    $iii++;
		    $count = $iii;
		    last;
		}
		else {
		    $iii++;
		}
	    }
	}
	push(@final_results, $filename_ref_short."\t".$filename_db_short."\t".$number_hits."\t".($sum_ANI / $number_hits)."\n");
	
	push(@{$dist->[$i+1]},($sum_ANI / $number_hits));
    }
}

open(COL_RESULTS, ">", $in_dir."/".$output_file."_aai_3col.txt");
open(DM, ">", $in_dir."/".$output_file."_aai_dm.txt");

foreach (@final_results) {
    print COL_RESULTS $_
}

foreach (@{$dist}) {
    foreach (@{$_}) {
	print DM $_."\t";
    }
    print DM "\n";
}
close(COL_RESULTS);
close(DM);

#########################################
#using predicted orthologs, calculate ANI from original gene sequences


opendir(D, "$in_dir_faa/temp_blast_results") or die "Can't open directory $in_dir_faa/temp_blast_results created by translating the gene sequences\n $!\n";
my @list2 = grep !/^\.\.?$/, readdir(D);
closedir(D);

my @blastp_orths;
foreach (@list2) {
    if ($_ =~ m/\.good_results/) {
	push(@blastp_orths, $_)
    }
}

foreach (@blastp_orths) {
    open(IN_RESULTS, "$in_dir/faa_out/temp_blast_results/$_");

    `mkdir $in_dir/temp_ortholog_files/`;
    `mkdir $in_dir/temp_ortholog_files/gene_seqs/`;
    my $file_header = $_;
    my @file_names = split('_vs_',$file_header);
    foreach (@file_names) {
	$_ =~ s/\.good_results//;
	$_ =~ s/\.faa//;
    }
    open(FILE1, ">", "$in_dir/temp_ortholog_files/$file_names[0].seq-list");
    open(FILE2, ">", "$in_dir/temp_ortholog_files/$file_names[1].seq-list");

    my @blastp_initial_results = <IN_RESULTS>;
    close(IN_RESULTS);
    print "\n";
    foreach (@blastp_initial_results) {
	chomp;
	my @fields = split('\t', $_);
	my @first_genes_to_find = $fields[0];
	my @second_genes_to_find = $fields[1];

	print FILE1 "$_\n" for (@first_genes_to_find);
	print FILE2 "$_\n" for (@second_genes_to_find);
    }

    `filter_fasta.py -f $in_dir/$file_names[0] -o $in_dir/temp_ortholog_files/gene_seqs/$file_names[0].ffn -s $in_dir/temp_ortholog_files/$file_names[0].seq-list`;
    `filter_fasta.py -f $in_dir/$file_names[1] -o $in_dir/temp_ortholog_files/gene_seqs/$file_names[1].ffn -s $in_dir/temp_ortholog_files/$file_names[1].seq-list`;

    
    close(FILE1);
    close(FILE2);
}


my $in_dir_ffn = "$in_dir/temp_ortholog_files/gene_seqs";
opendir(D, "$in_dir_ffn") or die "Can't open directory created by translating the gene sequences\n";
my @list3 = grep !/^\.\.?$/, readdir(D);
closedir(D);


my @ffn_list1;
foreach (@list3) {
    if ($_ =~ m/.*\.ffn/) {
	push(@ffn_list1, $_);
    }
}

my @sorted_list2 = sort(@ffn_list1);

my $i2;
my $dist2->[0] = [@sorted_list2];
unshift(@{$dist2->[0]}, 0);
my @final_results2;

for ($i2 = 0; $i2 < @sorted_list2; $i2++) {
    my $filename_ref = "$in_dir_ffn/".$sorted_list2[$i2];
    my $filename_ref_short = $sorted_list2[$i2];
    my @temp_array;
    
    for (my $j = 0; $j <= $i2; $j++) {
	push(@temp_array, 0);
    }
    push(@{$dist2->[$i2+1]}, $sorted_list2[$i2]); #add row labels
    push(@{$dist2->[$i2+1]}, @temp_array); #add correct number of zeros to set up adding ANI values to correct location
    
    for (my $k = $i2+1; $k < @sorted_list2; $k++) {
	my $filename_db = "$in_dir_ffn/".$sorted_list2[$k];
	my $filename_db_short = $sorted_list2[$k];
	
	#set up reference database
	system(
	    "if [ ! -d $in_dir_ffn/temp_blast_dbs ];
        then
            mkdir $in_dir_ffn/temp_blast_dbs
        fi");
	
	system("makeblastdb -in $filename_db -input_type 'fasta' -dbtype 'nucl' -out $in_dir_ffn/temp_blast_dbs/$filename_db_short");
	
	system(
	    "if [ ! -d $in_dir_ffn/temp_blast_results ]; 
        then
            mkdir $in_dir_ffn/temp_blast_results
        fi");
	
	#run first blast
	system("blastn -task dc-megablast -query $filename_ref -db $in_dir_ffn/temp_blast_dbs/$filename_db_short -out $in_dir_ffn/temp_blast_results/$filename_ref_short -max_target_seqs 1 -outfmt '6 qseqid sseqid length slen pident' -num_threads $cores");
	#set up reciprocal blast
	system("makeblastdb -in $filename_ref -input_type 'fasta' -dbtype 'nucl' -out $in_dir_ffn/temp_blast_dbs/$filename_ref_short");
	
	system("blastn -task dc-megablast -query $filename_db -db $in_dir_ffn/temp_blast_dbs/$filename_ref_short -out $in_dir_ffn/temp_blast_results/$filename_db_short -max_target_seqs 1 -outfmt '6 sseqid qseqid pident' -num_threads $cores");
	print "\n\n";
	#keep only reciprocal hits
	open(B_RESULTS, "$in_dir_ffn/temp_blast_results/$filename_ref_short");
	open(KEEP, ">", "$in_dir_ffn/temp_blast_results/$i.good_results");
	open(B_RECIP, "$in_dir_ffn/temp_blast_results/$filename_db_short");
	
	my @good_results;
	while(<B_RESULTS>) {
	    my @fields = split('\t', $_);
	    if ($fields[2] / $fields[3] > 0 && $fields[4] > 0) { ##### adjust sensitivity here
		push(@good_results, $_);
		print KEEP $_;
	    }
	}
	
	my @recip_results = <B_RECIP>;
	@good_results = sort @good_results;
	@recip_results = sort @recip_results;
	
	my $count //=0;
	my $number_hits = 0;
	my $sum_ANI = 0;
	my @recip_good_blasts;
	foreach (@good_results) {
	    chomp;
	    #print $_."\n";
	    my @fields = split('\t', $_);
	    my $cur_file_gene_name = $fields[0];
	    my $ref_gene_name = $fields[1];
	    my $percent_id = $fields[4];
	    
	    for (my $iii = $count; $iii < @recip_results;) {
		my @fields_2 = split('\t', $recip_results[$iii]);
		my $cur_gene_name2 = $fields_2[0];
		my $ref_gene_name2 = $fields_2[1];
		    if ($cur_file_gene_name eq $cur_gene_name2 && $ref_gene_name eq $ref_gene_name2) {
			#print $cur_file_gene_name."\t".$ref_gene_name."\t".$percent_id."\n";
			$sum_ANI += $percent_id;
			$number_hits++;
			$iii++;
			$count = $iii;
			last;
		    }
		    else {
			$iii++;
		    }
		}
	    }
	push(@final_results2, $filename_ref_short."\t".$filename_db_short."\t".$number_hits."\t".($sum_ANI / $number_hits)."\n");
	
	push(@{$dist2->[$i2+1]},($sum_ANI / $number_hits));
    }
}

open(COL_RESULTS, ">", $in_dir."/".$output_file."_ani_3col.txt");
open(DM, ">", $in_dir."/".$output_file."_ani_dm.txt");

foreach (@final_results2) {
    print COL_RESULTS $_;
    #print $_."\n";
}

foreach (@{$dist2}) {
    foreach (@{$_}) {
	print DM $_."\t";
    }
    print DM "\n";
}
close(COL_RESULTS);
close(DM);

