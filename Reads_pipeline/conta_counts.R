library(tidyverse)
# 7 jul. 2025
# Author: N. Godron (nicolas.godron {a} inserm.fr)
# Review round 1, with regard to remarks of Reviewer #2 on contamination threshold.
# Counting contaminated samples per family under MISAG high and medium contamination threshold.


distribution_conta <- read_tsv("~/Downloads/Distribution_table.tsv")
distribution_conta$MISAG_high <- distribution_conta$sra_abundance > 95
distribution_conta$MISAG_medium <- distribution_conta$sra_abundance > 90

## Histogram of SRA abundance distribution
hist(distribution_conta$sra_abundance, breaks = 100)
abline(v = 90, col = "blue")
abline(v = 95, col = "red")


conta_counts <- table(distribution_conta$MISAG_high,
                                 distribution_conta$family)
conta_counts <- rbind(conta_counts,
                     table(distribution_conta$MISAG_medium,
                           distribution_conta$family))
conta_counts <- data.frame(conta_counts)
row.names(conta_counts) <- c("cont_high-5", "uncont_high-5", "count_med-10", "uncount_med-10")

# Enterobacteriaceae
fisher.test(matrix(conta_counts$Enterobacteriaceae, nrow = 2))
# p-value = 0.5719

# Staphylococcaceae
fisher.test(matrix(conta_counts$Staphylococcaceae, nrow = 2))
# p-value = 0.01766

# Streptococcaceae
fisher.test(matrix(conta_counts$Streptococcaceae, nrow = 2))
# p-value = 0.414

# Total
fisher.test(matrix(c(134,47,114,67), nrow =2))
# p-value = 0.03133