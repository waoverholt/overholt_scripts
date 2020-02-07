#!/bin/bash 

#Directory containing F & R reads of files to assemble
INDEX=$HOME/Projects/GATech/Dissertation/metaG/06.annotation/ROCker/T222C1_amoA/amoA_assembly_T222C1/bowtie2_contig_index
INDIR=~/Projects/GATech/Dissertation/metaG/metaG_concat_lanes_paired_fastq/
OUTDIR=~/Projects/GATech/Dissertation/metaG/06.annotation/ROCker/T222C1_amoA/bowtie2_mapping/
PROCS=2
FILES=($INDIR/*)

if [ ! -d "$OUTDIR" ]; then
    mkdir "$OUTDIR";
fi

if [ ! -d "$INDIR" ]; then
    echo "ERROR"
    echo "$INDIR doesn't exist, please check your path and try again"
    exit 1
fi


for FILE in ${FILES[*]};
do
    ext=$(basename $FILE);
    temp_file_name=${ext%.*};
    file_name=${temp_file_name%.*};
    list_sample_name+=($file_name);
done

unique_list=($(echo "${list_sample_name[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

FORWARD=""
REVERSE=""

count=0
for sampleid in ${unique_list[@]};
do
    for FILE in ${FILES[@]};
    do
	if `echo $FILE | grep -q $sampleid`
	then
	    if `echo $FILE | grep -q "R1"`
	    then
		FORWARD=$FILE
		count=$((count+1))
	    fi

	    if `echo $FILE | grep -q "R2"`
	    then
		REVERSE=$FILE
		count=$((count+1))
	    fi
	    if [ "$count" -eq "2" ]
	    then
		echo "Going to run bwa-mem on the following files:"
		echo $FORWARD
		echo $REVERSE
		echo "using the index: $INDEX"
		echo
		samp_name=$(basename $sampleid)
		samp_name=${samp_name%%_*}
		OUTFILE=$samp_name".sam"
		OUTFILE="$OUTDIR$OUTFILE"	
		#echo $FORWARD, $REVERSE, $INDEX
		bowtie2 -f -x $INDEX -1 $FORWARD -2 $REVERSE -S $OUTFILE -p $PROCS --no-unal
		samtools view -F 4 -bS $OUTFILE -o $OUTFILE.bam
		anvi-init-bam $OUTFILE.bam -o $OUTFILE.anvi.bam
		count=0
	    fi
	fi
    done
done

