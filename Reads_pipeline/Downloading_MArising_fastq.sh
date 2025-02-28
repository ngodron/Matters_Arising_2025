#!/bin/bash
main_folder="/home/nicolas/NiGo_PhD/MArising/"
SRS_ID_file=${main_folder}"/data/SRA_confirmed.IDs"
ENA_URL_file=${main_folder}"/data/SRA_confirmed.info"

cd ${main_folder}/data
touch "$ENA_URL_file"

while read SRA
do
	echo $SRA
	ena_url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRA}&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,sra_ftp&format=tsv&download=true"
	curl $ena_url >> "$ENA_URL_file" 
done < ${SRS_ID_file}

cat $ENA_URL_file | grep "fastq.gz" | tr ";" "	" >"${ENA_URL_file}2"

awk 'BEGIN { FS = "\t"; OFS = "\n" } { for (i = 1; i <= NF; i++) \
 { while (match($i, /ftp\.[^ ]*\.fastq\.gz/)) { print substr($i, RSTART, RLENGTH); \
  $i = substr($i, RSTART + RLENGTH) } } }' "${ENA_URL_file}2" > "${ENA_URL_file}.URLs"

cd ${main_folder}
mkdir -p ./fastq && cd ./fastq

while read url
do
	echo $url
	wget $url
done < "${ENA_URL_file}.URLs"

echo "All done sir! (No guaranteed results, MIT license kind-of stuff)"
