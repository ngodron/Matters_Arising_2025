# Matters_Arising_2025
Bioinformatics scripts to analyze data from Ellabaan et al., Nature Comms (2021)

## Reads_pipeline
### Downloading FASTQ files
- Script: Downloading_MArising_fastq.sh, calls fastq-dl

### Trimming with Trimmomatic
- Parameters:
LEADING:20 TRAILING:20 SLIDINGWINDOW:15:20

- Scripts: trimming.sh which calls trimmomaticparallel.sh

### Taxonomic assignation with MetaPhlAn 4.0.2
- Database: mpa_vJun23_CHOCOPhlAnSGB_202403
- Scripts: metaphlan_cdc.sh, merge_metaphlan_tables.py (from MetaPhlAn authors)

### Analysis
- Script: plots_family.R

## Assembly_pipeline
### ARG & MGE queries
- Snakemake script: SNAKEcontasearch-onAssembly

### Analysis
- Script: plot_assembled_contigs.R

## Supplementary_data
### Common to both pipelines (metadata & download)
- SRAlist.txt
List of SRA genomes in confirmation step of original study.

- SRAARGMapping.tab
MGE-ARG pairings found in these SRA genomes.

### Assembly_pipeline
- Suppl_Table_1.txt
Count of ARGs, Taxonomic family and SRA ID.

### Reads_pipeline
- All_run_info.tsv
SRA run metadate, to allow download and metadata parsing.

- Distribution_table.tsv
Aggregated results from MetaPhlAn taxonomic assignation at the family level. 

- Threshold_family.tsv
Presence/absence matrix at the family level (threshold for presence: > 5% abundance) and correspondance with SRA metadata (binned, taxonomic mismatch corresponds to < 5% of expected family).

- plots_family_RsessionInfo.txt
Info about the version of R (4.3.3) packages for reproducibility purposes.
