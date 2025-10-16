library(tidyverse)

main_dir <- getwd()

GUNC_df <-
  read_tsv(file = paste0(main_dir,
                         "/output/GUNC/GUNC.progenomes_2.1.maxCSS_level.tsv"))

quantile(GUNC_df$contamination_portion,
         probs = seq(0,1,0.1))


MPhLan_df <-
  read_tsv(file = paste0(main_dir,
                         "/Github/Matters_Arising_2025/Supplementary_data/Distribution_table.tsv"))


merged_df <- full_join(GUNC_df, MPhLan_df, by = join_by (genome == Name))

Entero_df <- merged_df[merged_df$family == "Enterobacteriaceae", ]
Staph_df <- merged_df[merged_df$family == "Staphylococcaceae", ]
Strept_df <- merged_df[merged_df$family == "Streptococcaceae", ]

### ECDF values ----
cdf_GUNC_Entero <-
  ecdf(Entero_df$contamination_portion[order(Entero_df$contamination_portion)] * 100)
cdf_GUNC_Staph <-
  ecdf(Staph_df$contamination_portion[order(Staph_df$contamination_portion)] * 100)
cdf_GUNC_Strept <-
  ecdf(Strept_df$contamination_portion[order(Strept_df$contamination_portion)] * 100)

cdf_MPhLan_Entero <-
  ecdf(Entero_df$sra_abundance[order(Entero_df$sra_abundance)])
cdf_MPhLan_Staph <-
  ecdf(Staph_df$sra_abundance[order(Staph_df$sra_abundance)])
cdf_MPhLan_Strept <-
  ecdf(Strept_df$sra_abundance[order(Strept_df$sra_abundance)])

# Sample counts
table(MPhLan_df$family)
table(merged_df[! is.na(merged_df$contamination_portion), "family"])

### GUNC ECDF plots ----
par(mfrow=c(2,3))
plot(cdf_GUNC_Entero, main = "GUNC Enterobacteriaceae ECDF (contigs, n=50)",
     xlab = "GUNC-estimated contamination values", ylab = "ECDF",
     xlim = c(0,100), xaxt= 'n')
axis(side = 1, at = seq(0,100,5))
abline(v=c(5,10,15), col = c("red","orange","green4"), lty = 3, lwd = 3)

plot(cdf_GUNC_Staph, main = "GUNC Staphylococcaceae ECDF (contigs, n=90)",
     xlab = "GUNC-estimated contamination values", ylab = "ECDF",
     xlim = c(0,100), xaxt= 'n')
axis(side = 1, at = seq(0,100,5))
abline(v=c(5,10,15), col = c("red","orange","green4"), lty = 3, lwd = 3)

plot(cdf_GUNC_Strept, main = "GUNC Streptococcaceae ECDF (contigs, n=23)",
     xlab = "GUNC-estimated contamination values", ylab = "ECDF",
     xlim = c(0,100), xaxt= 'n')
axis(side = 1, at = seq(0,100,5))
abline(v=c(5,10,15), col = c("red","orange","green4"), lty = 3, lwd = 3)


### MetaPhLan ECDF ----
plot(cdf_MPhLan_Entero, main = "MetaPhLan Enterobacteriaceae ECDF (reads, n=61)",
     xlab = "SRA concordance values", ylab = "ECDF",
     xlim = c(0,100), xaxt= 'n')
axis(side = 1, at = seq(0,100,5))
abline(v=c(85,90,95), col = c("green4","orange","red"), lty = 3, lwd = 3)

plot(cdf_MPhLan_Staph, main = "MetaPhLan Staphylococcaceae ECDF (reads, n=97)",
     xlab = "SRA concordance values", ylab = "ECDF",
     xlim = c(0,100), xaxt= 'n')
axis(side = 1, at = seq(0,100,5))
abline(v=c(85,90,95), col = c("green4","orange","red"), lty = 3, lwd = 3)

plot(cdf_MPhLan_Strept, main = "MetaPhLan Streptococcaceae ECDF (reads, n=23)",
     xlab = "SRA concordance values", ylab = "ECDF",
     xlim = c(0,100), xaxt= 'n')
axis(side = 1, at = seq(0,100,5))
abline(v=c(85,90,95), col = c("green4","orange","red"), lty = 3, lwd = 3)

# dev.off()

cor(x = merged_df$contamination_portion,
    y = merged_df$sra_abundance,
    use = "complete.obs",
    method = "pearson")
# [1] 0.1927582

cor(x = Entero_df$contamination_portion,
    y = Entero_df$sra_abundance,
    use = "complete.obs",
    method = "pearson")

cor(x = Staph_df$contamination_portion,
    y = Staph_df$sra_abundance,
    use = "complete.obs",
    method = "pearson")

cor(x = Strept_df$contamination_portion,
    y = Strept_df$sra_abundance,
    use = "complete.obs",
    method = "pearson")

plop <- as.vector(MPhLan_df[MPhLan_df$family == "Staphylococcaceae", "sra_abundance"])

ggplot(MPhLan_df) +
  geom_density(aes(x = sra_abundance,
                   fill = family),
               color = "black",
               position = "identity",
               alpha = 0.5) +
  geom_vline(aes(xintercept=95), color="red") +
  geom_vline(aes(xintercept=90), color="orange") +
  geom_vline(aes(xintercept=85), color="green4") +
  ggtitle("Density of MetaPhLan-to-SRA concordance (reads, n=181)") +
  xlab("MetaPhLan abundance of the SRA metadata family (%)") +
  scale_x_continuous(breaks = seq(0,100,5))

# GUNC_density_df <- 
#   merged_df[complete.cases(merged_df[, c("contamination_portion", "family")]),
#                                      c("contamination_portion", "family")]
# ggplot(GUNC_density_df) +
#   geom_density(aes(x = (contamination_portion*100),
#                    fill = family),
#                color = "black",
#                position = "identity",
#                alpha = 0.5) +
#   geom_vline(aes(xintercept=5), color="red") +
#   geom_vline(aes(xintercept=10), color="orange") +
#   geom_vline(aes(xintercept=15), color="green4") +
#   ggtitle("Density of the distribution\nof GUNC-estimated contamination (contigs)") +
#   xlab("GUNC-estimated contamination")
