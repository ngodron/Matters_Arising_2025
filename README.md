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
- All_run_info.tsv:

- Distribution_table.tsv:

- SRAARGMapping.tab:

- SRAlist.txt:

- Suppl_Table_1.txt:

- Threshold_family.tsv:

- plots_family_RsessionInfo.txt:

