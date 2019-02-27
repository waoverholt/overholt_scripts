Will Overholt
12-18-13

#convert to individual read files
ls /nv/hp10/kmarks6/data/raw_data/raw_seqs/ | grep R1 | xargs -I file cat /nv/hp10/kmarks6/data/raw_data/raw_seqs/file >> MSU_concat_raw_reads_R1.fastq

ls /nv/hp10/kmarks6/data/raw_data/raw_seqs/ | grep R2 | xargs -I file cat /nv/hp10/kmarks6/data/raw_data/raw_seqs/file >> MSU_concat_raw_reads_R2.fastq

#join paired end reads using qiime
join_paired_ends.py -f MSU_concat_raw_reads_R1.fastq -r MSU_concat_raw_reads_R2.fastq -b MSU_concat_barcodes.fastq

#make barcode file
extract_barcodes.py -f fastqjoin.join.fastq -c barcode_in_label --char_delineator ':' --bc1_len 12 -o barcodes

#deplex with qiime (no quality filter)
split_libraries_fastq.py -q 0 --store_demultiplexed_fastq -m ../mapping_corrected.txt --barcode_type golay_12 -b barcodes.fastq -i fastqjoin.join.fastq -o qiime_deplex
*note barcodes are not reserved in MSU seqs

#quality filter with usearch1.7.0
ls | xargs -n 1 -P 4 -I file usearch1.7.0 -fastq_filter file -fastq_maxee 0.5 -fastaout file.out -eeout

#merge together
cat seqs.1.fastq.out seqs.2.fastq.out seqs.3.fastq.out seqs.4.fastq.out >> ../qc_seqs.fna

#run chimera checking through usearch61 (used job script)
msub job_scripts/qiime_usearch61_chimera.pbs
export INPUT=$HOME/data/RAW_data_files/MSU_raw_seqs_11-15-13/join-paired/qiime_deplex/qc_seqs.fna
export OUTPUT=$HOME/data/RAW_dat_files/MSU_raw_seqs_11-15-13/join-paired/qiime_deplex/usearch61_chimera_checking
identify_chimeric_seqs.py -i $INPUT -r $GG_otus -m usearch61 --split_by_sampleid -o $OUTPUT

#discard chimeras
filter_fasta.py -f qc_seqs.fna -o qc_seqs_no_chim.fna -s usearch61_chimera_checking/non_chimeras.txt

#dereplicate using prefix suffix filter
pick_otus.py -m prefix_suffix -i qc_seqs_no_chim.fna -o prefix_picked_otus

#pick rep set of deplicated seqs
pick_rep_set.py -f ../qc_seqs_no_chim.fna -i qc_seqs_no_chim_otus.txt -o qc_seqs_no_chim_pre-suf-rep.fna

#cluster using usearch61 [need to think of way to use usearch1.7??]
pick_otus.py -m usearch61 -o usearch61_otus -s 0.97 -z -i qc_seqs_no_chim_pre-suf-rep.fna -k

#next steps
1. merge otu maps from pre-cluster and usearch61 cluster
2. go back, use discarded chimeras file (qc_seqs_no_chim.fna) and pull out only deep sea. Finish this process with the other 2 seq files from deep sea. Include Liu's data.
