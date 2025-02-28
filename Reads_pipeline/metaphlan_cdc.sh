#!/bin/bash

# docker run --name "NG_metaphlan" -v "/home/nicolas.godron/:/home/nicolas.godron/" cdc.local/metaphlan:4.0.2 /home/nicolas.godron/IAME/NiGo_PhD/MArising/scripts/metaphlan_cdc.sh

THREADS=24
FASTQ_EXTENSION="fastq.gz"

DB_FOLDER="/home/nicolas.godron/THESE-KLA/2024_Final/NG_allQC/2.2-metaphlan/metaDB/metaphlan_databases2/"

IN="/home/nicolas.godron/IAME/NiGo_PhD/MArising/run_2/trimmomatic/"
IN_SINGLE="${IN}single_read/"
IN_PAIR="${IN}pair/"
ID_LIST="${IN}../ID_list.txt"

OUT="/home/nicolas.godron/IAME/NiGo_PhD/MArising/run_2/metaphlan/"
LOG="${OUT}../metaphlan_run.log"

cd $IN_PAIR
# Renaming strains
mv 394-2_staphylococcus-aureus.R1.fastq.gz 394-2_staphylococcus-aureus_1.fastq.gz
mv 394-2_staphylococcus-aureus.R2.fastq.gz 394-2_staphylococcus-aureus_2.fastq.gz

# Making list of (unique) sample IDs
# touch $ID_LIST
# for fastq_file in $(ls)
# do
# 	echo $fastq_file | sed -E 's/(_1)?(_2)?.fastq.gz//g' >> $ID_LIST
# done

cd $IN_SINGLE
for fastq_file in $(ls)
do
	echo $fastq_file | sed -E 's/(_1)?(_2)?.fastq.gz//g' >> $ID_LIST
done
cat $ID_LIST | sort | uniq | tee $ID_LIST


mkdir -p $OUT/bowtie2/ $OUT/metaphlan_profiles/

time metaphlan --install --index mpa_vJun23_CHOCOPhlAnSGB_202403 --bowtie2db $DB_FOLDER
echo "time metaphlan --install --index mpa_vJun23_CHOCOPhlAnSGB_202403 --bowtie2db $DB_FOLDER" > $LOG

while read sample
do
	echo $sample >> $LOG
	# Testing existence of paired-end fastq
	# if [ -f ${IN_PAIR}${sample}_1.${FASTQ_EXTENSION} ] && [ -f ${IN_PAIR}${sample}_2.${FASTQ_EXTENSION} ]; then
	# 	echo "Paired-end fastq files exist!" >> $LOG

	# 	time metaphlan "${IN_PAIR}${sample}_1.${FASTQ_EXTENSION}","${IN_PAIR}${sample}_2.${FASTQ_EXTENSION}" \
	# 	--input_type fastq \
	# 	--index mpa_vJun23_CHOCOPhlAnSGB_202403 \
	# 	--bowtie2db $DB_FOLDER \
	# 	--bowtie2out ${OUT}/bowtie2/${sample}.bowtie2.bz2 \
	# 	--nproc $THREADS \
	# 	-o ${OUT}/metaphlan_profiles/${sample}_metaphlan.txt
	
	if [ -f ${IN_SINGLE}${sample}.${FASTQ_EXTENSION} ]; then
		echo "Single read fastq file exist! (At least one paired-end file not found)." >> $LOG
		time metaphlan "${IN_SINGLE}${sample}.${FASTQ_EXTENSION}" \
		--input_type fastq \
		--index mpa_vJun23_CHOCOPhlAnSGB_202403 \
		--bowtie2db $DB_FOLDER \
		--bowtie2out ${OUT}/bowtie2/${sample}.bowtie2.bz2 \
		--nproc $THREADS \
		-o ${OUT}/metaphlan_profiles/${sample}_metaphlan.txt
	else
		echo "Error: no fastq file found (paired-end or single read). Maybe a single and lonely paired-end file (forward or reverse) is present." >> $LOG
	fi
done < $ID_LIST

