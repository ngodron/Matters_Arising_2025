#!/bin/bash

# MAIN_DIR="/home/nicolas/Work/CatiBioMed/NiGo_PhD/MArising/"
# docker run -v $MAIN_DIR:/main cdc.local/quality_parallel /main/scripts/trimming.sh

THREADS=12
FASTQ_EXTENSION="fastq.gz"
PARALLEL_PATH="/parallel-20240722/src/parallel"

IN="/main/fastq/"
ID_LIST="/main/output/IDs/ID_list.txt"

OUT="/main/run_2/trimmomatic/"

mkdir -p $OUT/input_quality
mkdir -p $OUT/pair/quality $OUT/unpair/quality $OUT/single_read/quality

echo "Trimming MArising strains, run start : $(date)" > "${OUT}/trimmomatic_run.log"

awk -v infolder="$IN" -v outfolder="$OUT" -v fastqext="$FASTQ_EXTENSION" \
	'{print infolder,outfolder,fastqext,$0}' "$ID_LIST" | $PARALLEL_PATH -j $THREADS 'echo {1} {2} {3} {4}; /main/scripts/trimmomaticparallel.sh {1} {2} {3} {4}'

echo "Trimming MArising strains, run end : $(date)" >> "${OUT}/trimmomatic_run.log"
