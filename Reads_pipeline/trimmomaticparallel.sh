#!/bin/bash

trimmomaticparallel () {
	IN=$1
	OUT=$2
	FASTQ_EXTENSION=$3
	sample=$4

	# Testing existence of paired-end fastq
	if [ -f ${IN}${sample}_1.${FASTQ_EXTENSION} ] && [ -f ${IN}${sample}_2.${FASTQ_EXTENSION} ]; then
		echo "Sample ${sample}: Paired-end fastq files exist!" >> "${OUT}/trimmomatic_run.log"
		
		read1="${sample}_1.${FASTQ_EXTENSION}"
		read2="${sample}_2.${FASTQ_EXTENSION}"

		fastqc -o $OUT/input_quality \
			$IN/$read1 \
			$IN/$read2

		read1_pair="${sample}_pair_1.${FASTQ_EXTENSION}"
		read2_pair="${sample}_pair_2.${FASTQ_EXTENSION}"
		read1_unpair="${sample}_unpair_1.${FASTQ_EXTENSION}"
		read2_unpair="${sample}_unpair_2.${FASTQ_EXTENSION}"

		java -jar /opt/progs/Trimmomatic/trimmomatic.jar PE \
		 	$IN/$read1 $IN/$read2 \
		 	$OUT/pair/$read1_pair $OUT/unpair/$read1_unpair \
		 	$OUT/pair/$read2_pair $OUT/unpair/$read2_unpair \
		 	LEADING:20 TRAILING:20 SLIDINGWINDOW:15:20

		fastqc -o $OUT/pair/quality \
			$OUT/pair/$read1_pair \
			$OUT/pair/$read2_pair

		fastqc -o $OUT/unpair/quality \
			$OUT/unpair/$read1_unpair \
			$OUT/unpair/$read2_unpair
		
	elif [ -f ${IN}${sample}.${FASTQ_EXTENSION} ]; then
		echo "Sample ${sample}: Single read fastq file exist! (At least one paired-end file not found)." >> "${OUT}/trimmomatic_run.log"

		read="${sample}.${FASTQ_EXTENSION}"

		fastqc -o $OUT/input_quality \
			$IN/$read

		java -jar /opt/progs/Trimmomatic/trimmomatic.jar SE \
		 	$IN/$read \
		 	$OUT/single_read/$read \
		 	LEADING:20 TRAILING:20 SLIDINGWINDOW:15:20

		fastqc -o $OUT/single_read/quality \
			$OUT/single_read/$read

	else
		echo "Sample ${sample}: Error: no fastq file found (paired-end or single read)." >> "${OUT}/trimmomatic_run.log"
	fi
}

trimmomaticparallel $1 $2 $3 $4
