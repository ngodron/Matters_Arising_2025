mkdir -p 03.MAPPED_IS1
mkdir -p 04.MAPPING/DB
mkdir -p 05.COVERAGE
mkdir -p 06.KRAKEN
mkdir -p 07.PLOTS

# To run after collection of genomic assemblies and SRA reads.
#

for genome in `ls -1 02.ASSEMBLIES/*.fna`; do
 name=`basename $genome .fna`

 # search for IS1 in contigs
 blastn -subject $genome -query IS1_Sequences.fasta -outfmt 6 -perc_identity 90 -out 03.MAPPED_IS1/$name.IS1.out

 # create bowtie database for mapping
 bowtie2-build $genome 04.MAPPING/DB/$name

 # combine all reads in a single file
 gzip -d 00.READS/$name/*.gz && cat 00.READS/$name/*.fastq > 00.READS/$name/$name.allreads.fastq

 
 # perform mapping
 bowtie2 -x 04.MAPPING/DB/$name --local -p 20 -U 00.READS/$name/$name.allreads.fastq -S 04.MAPPING/$name.mapped.sam --score-min L,10,0.8 >> 04.MAPPING/stats_Mapping.$name.log 2>> 04.MAPPING/stats_Mapping.$name.log
 sam2bam 04.MAPPING/$name.mapped.sam $name 10
 samtools coverage 04.MAPPING/$name.mapped.bam --ff UNMAP | sort -n -k 3 > 05.COVERAGE/$name.coverage.out
 
 # run kraken (done on Migale) and extract important information
 kraken2 --db /db/outils/kraken2-2025/bacteria_k25 --use-names --threads 10 02.ASSEMBLIES/$name.fna | sed \"s/\\t/_/g\" > 06.KRAKEN/$name.kk.out
 cut -d "_" -f 1,2,3,4 06.KRAKEN/$name.kk.out > 06.KRAKEN/$name.kkshort.out
 
 # display plot
 echo --------------------------------------------
 echo create plot for $name
 echo --------------------------------------------
 echo ############################################
 Rscript barplotCoverage_IS.R $name 05.COVERAGE/$name.coverage.out 03.MAPPED_IS1/$name.IS1.out 06.KRAKEN/$name.kkshort.out 07.PLOTS/$name.barplot.png
 
done


