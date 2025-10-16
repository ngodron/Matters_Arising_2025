#!/usr/bin/bash

SAM=$1
threads=$2
BAM=`echo $1 | sed 's/sam$/bam/i'`
echo $BAM

if [ -z $2 ]
then
 threads=2
fi 

#id=$3
#samtools view -@ $threads -bS -h $SAM | bamtools sort > /tmp/TMP$id.bam
#bamaddrg -b /tmp/TMP$id.bam -s $id > $BAM

samtools view -@ $threads -bS -h $SAM | bamtools sort > $BAM
bamtools index -in $BAM

