# Matters_Arising_2025
Bioinformatics scripts to analyze data from Ellabaan et al., Nature Comms (2021)

## Reads_pipeline
### Downloading FASTQ files:
- Script: Downloading_MArising_fastq.sh, calls fastq-dl

### Trimming: Trimmomatic
- Parameters:
LEADING:20 TRAILING:20 SLIDINGWINDOW:15:20

- Scripts: trimming.sh which calls trimmomaticparallel.sh

### Taxonomic assignation: MetaPhlAn 4.0.2
- Database: mpa_vJun23_CHOCOPhlAnSGB_202403
- Scripts: metaphlan_cdc.sh, merge_metaphlan_tables.py (from MetaPhlAn authors)

### Analysis:
- Script: plots_family.R
