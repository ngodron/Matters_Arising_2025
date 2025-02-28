# 16 sept. 2024
# Author: N. Godron (nicolas.godron {a} inserm.fr)
# Visuals for metaphlan interpretation in light of SRA metadata.

# Libraries ----
library(tidyverse)
library(taxize)
library(ggbeeswarm)
library(ComplexUpset)

# Parameters ----
# Expressed in % of relative abundance (as in metaphlan output)
cont_threshold <- 5

# Inputs ----
main_dir <- "/home/nicolas/Work/CatiBioMed/NiGo_PhD/MArising/"

## metadata ----
metadata <-
  read_tsv(paste0(main_dir, "output/metadata/All_run_info.tsv"))

too_many_columns <- problems(metadata)
# Upon examination, 6 more columns are found in 65 metadata rows
# The columns are not of interest to this work and are thus discarded.
# They are: library_prep_location,	rna_prep_3_protocol, ph,
#   sequencing_longitude,	tissue_type, and	isolation_source
rm(too_many_columns)

## metaphlan ----
contamination_df <-
  read_tsv(paste0(main_dir,
                  "run_2/abundance_table.tsv"),
           comment = "#")

# Filtering, parsing & querying taxonomy ----
## metadata, GNR ----
metadata_df <-
  metadata |>
  select(c("run_accession", "sample_alias", "scientific_name",
           "tax_id", "tax_lineage",
           "sample_accession", "submission_accession", 
           "study_accession", "experiment_accession",
           "library_selection", "fastq_bytes", "description",
           "read_count", "base_count", "library_strategy", "library_layout",
           "tag", "first_created"))

# Two samples are in the table twice due to aliases (cf. metadata_df)
duplicates <- data.frame(contamination_df$clade_name,
                         contamination_df$ERR594079, contamination_df$"0788",
                         contamination_df$ERR2204715, contamination_df$"7518")
( duplicates <- duplicates[rowSums(duplicates[2:5] != 0) > 0, ] )
# Filtering out these three sequences
contamination_df <- contamination_df |>
  select(! c("0788", "7518")) |>
  # Renaming a strain with a super long name with its ERR ID instead
  rename(ERR1549309 = `Isolat27_IonXpress_012_R_2012_11_16_09_40_45_user_SN1-88_Auto_user_SN1-88_89`)


# Global Names Resolver (GNR) service from the "Encyclopedia of Life"
resolved_names <- taxize::gnr_resolve(sci = metadata_df$scientific_name,
                                      # data source n°12 is the "Encyclopedia of Life"
                                      data_source_ids = 12,
                                      resolve_once = TRUE,
                                      best_match = TRUE)

names(resolved_names)[c(1,3)] <- c("scientific_name", "resolved_name")

metadata_df <-
  left_join(metadata_df,
            resolved_names[c("scientific_name", "resolved_name")],
            by = "scientific_name") |>
  rowwise() |>
  mutate(genus_species = word(resolved_name,
                              start = 1,
                              end = 2,
                              sep = " ")) |>
  mutate(sra_genus = word(resolved_name,
                          start = 1,
                          end = 1,
                          sep = " "))

# Special case for resolved name of length 1 word (i.e. only Genus)
metadata_df$genus_species[is.na(metadata_df$genus_species)] <-
  metadata_df$resolved_name[is.na(metadata_df$genus_species)]

### Fetching family ----
# The step below can take a minute
unique_metadata_family <- tax_name(unique(metadata_df$genus_species),
                                   get = c("family","phylum"))
metadata_df <- left_join(x = metadata_df,
                         y = unique_metadata_family[2:4],
                         by = join_by(genus_species == query) )

## metaphlan ----
### Family ----
contamination_family <-
  contamination_df[grepl(x = contamination_df$clade_name,
                         pattern = "f__"), ]
contamination_family <-
  contamination_family[! grepl(x = contamination_family$clade_name,
                               pattern = "g__"), ]

contamination_family <- contamination_family |>
  mutate(phylum = str_replace(string = clade_name,
                              pattern = ".*p__(.*)\\|c__.*",
                              replacement = "\\1")) |>
  relocate(phylum,
           .after = "clade_name")

# Removing prefix from species name
contamination_family$clade_name <-
  str_replace_all(string = contamination_family$clade_name,
                  pattern = ".*f__",
                  replacement = "")

clade_name <- contamination_family$clade_name

## Test condition: "> contamination threshold?" ----
threshold_family <-
  data.frame(
    ifelse(test = contamination_family[3:ncol(contamination_family)] > cont_threshold,
           yes = 1,
           no = 0))

threshold_family <- data.frame(t(threshold_family))
colnames(threshold_family) <- clade_name

## Joining metadata and metaphlan ----
metaphlan_sra_ref <- t(contamination_family)
colnames(metaphlan_sra_ref) <- clade_name

# Left join (keeping all IDs in metaphlan_sra_ref) by row names
metaphlan_sra_ref <-
  merge(data.frame(metaphlan_sra_ref),
        data.frame(metadata_df[, c("run_accession", "family")]),
        # Here, 0 refers to row.names()
        by.x = 0,
        by.y = "run_accession",
        all.x = TRUE)

names(metaphlan_sra_ref)[1] <- "Name"

# How abundant is the family found in SRA metadata?
metaphlan_sra_ref$sra_abundance <-
  apply(X = metaphlan_sra_ref,
        MARGIN = 1,
        FUN = function(x) round(as.numeric(x[(x["family"])]),
                                digits = 2))

# For family not found in result table
metaphlan_sra_ref$sra_abundance[is.na(metaphlan_sra_ref$sra_abundance)] <- 0

# Left join to add info on the abundance of the family found in SRA metadata
threshold_family$Name <- row.names(threshold_family)
threshold_family <- merge(threshold_family,
                          metaphlan_sra_ref[, c("Name", "sra_abundance", "family")],
                          by = "Name",
                          # Left join:
                          all.x = TRUE)

# Removing columns of families without any sample surpassing threshold
# Non-numeric columns are kept.
threshold_family <- threshold_family[colSums(threshold_family == 0) < nrow(threshold_family)]
threshold_family$family <- as.factor(threshold_family$family)

threshold_family <- threshold_family |>
  mutate(SRA_Agreement = case_when(between(sra_abundance, 95, 100) ~ "95-100%",
                                   between(sra_abundance, 80, 95) ~ "80-95%",
                                   between(sra_abundance, 5, 80) ~ "5-80%",
                                   between(sra_abundance, 0, 5) ~ "Taxonomic mismatch",
  ))

family_and_phylum <- contamination_family[contamination_family$clade_name %in% colnames(threshold_family),
                                          c("clade_name", "phylum")]
colnames(family_and_phylum) <- c("set", "phylum")
# Ordering by number of families in phylum and then by number of samples per family
family_and_phylum$least_families_in_phylum <- match(family_and_phylum$phylum, names(sort(table(family_and_phylum$phylum))))
family_and_phylum$samples <- colSums(threshold_family != 0)[family_and_phylum$set]
family_and_phylum <-
  family_and_phylum[with(family_and_phylum, order(least_families_in_phylum, samples)),]

# Plots ----
## metaphlan-SRA congruency ----
violin_df <- metaphlan_sra_ref[! metaphlan_sra_ref$Name %in% c("clade_name", "phylum"),]
violin_df <- violin_df[order(violin_df$family), ]
violin_plot <- ggplot(data = violin_df) +
  aes(x = "",
      y = sra_abundance) +
  # geom_violin(adjust = 1e-1,
  #             alpha = 0.2,
  #             fill = "#aaaaaa") +
  ggbeeswarm::geom_beeswarm(aes(x = family,
                                y = sra_abundance,
                                colour = family),
                            method = "center",
                            cex = 2,
                            side = 0L,
                            dodge.width = 1) +
  scale_colour_manual(values = c("#ffb3f2", "#cea30d", "#003a9e")) +

  xlab("") +
  ylab("Observed relative abundance (%)\nof the expected family") +
  ylim(c(0, 100)) +
  scale_y_continuous(breaks = c(0,5,25,50,75,95,100)) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) +
  annotate("rect", fill = "#64ab05",
           xmin = 0.4, xmax = 0.5,
           ymin = 95, ymax = 100) +
  annotate("rect", fill = "#f3d819",
           xmin = 0.4, xmax = 0.5,
           ymin = 80, ymax = 95) +
  annotate("rect", fill = "#c86c44",
           xmin = 0.4, xmax = 0.5,
           ymin = 5, ymax = 80) +
  annotate("rect", fill = "#9c2a17",
           xmin = 0.4, xmax = 0.5,
           ymin = 0, ymax = 5)
violin_plot

## metaphlan UpSet ----
upset_plot <- ComplexUpset::upset(
  ### Basic input & parameters ----
  data = threshold_family,
  intersect = family_and_phylum$set,
  height_ratio = 0.8,
  # name refers to the main x label (below the intersection grid) 
  name = "",
  # Removing empty intersections
  min_degree = 1,
  sort_sets = FALSE,
  sort_intersections_by = c('cardinality', 'degree', 'ratio'),
  
  ### Intersections ----
  base_annotations = list(
    "Intersection size" = intersection_size(
      text_mapping = aes(
        label = !!get_size_mode("exclusive_intersection"),
        # Shenanigan to place label on top of bar, there must exist a better way.
        y = ifelse(!!get_size_mode("exclusive_intersection"),
                   !!get_size_mode("exclusive_intersection") + .5,
                   0)),
      text_colors = c(on_background = "black",
                      on_bar = "black"),
      text = list(size = 5),
      mapping = aes(fill = SRA_Agreement)) +
      scale_fill_manual(values = c("95-100%"="#64ab05",
                                   "80-95%"="#f3d819",
                                   "5-80%"="#c86c44",
                                   "0-5%"="#9c2a17")) +
      ylab("") +
      theme(axis.title =
              element_text(size = 14,
                           face = "bold",
                           hjust = 1),
            legend.position = "none")
  ),
  
  ### Graphical adjustments ----
  set_sizes = (
    upset_set_size(
      geom = geom_bar(
        mapping = aes(fill = family),
        width = 0.8))) +
    ylim(c(0,120)) +
    scale_fill_manual(values= c("Staphylococcaceae"="#cea30d",
                                "Enterobacteriaceae"="#ffb3f2",
                                "Streptococcaceae"="#003a9e")) +
  # mapping = aes(y = after_stat(count)))) +
    ylab(paste0("Count of samples with > ",
                cont_threshold, "%\nobserved abundance of bacterial family")) +
    theme(axis.ticks.x = element_line(),
          axis.title = element_text(size = 11),
          legend.position = "none" ),
  
  themes = upset_modify_themes(to_update = list(
    "intersections_matrix" = theme(text = element_text(size = 15)))
  )
)
upset_plot

# Saving plots ----
png(filename = paste0(main_dir, "/paper/Distribution_observed_abundance_v3.png"),
    width = 1600, height = 1200, res = 300,
    pointsize = 2)
violin_plot
dev.off()

png(filename = paste0(main_dir, "/paper/UpSet_plot_v3.png"),
    width = 3400, height = 2400, res = 300)
upset_plot
dev.off()

# Legend info for Staph. examples ----
# Examples of expected Staphylococcaceae with problematic contamination
Staph_examples_df <-
  contamination_df[c("clade_name", "SRR1522630", "ERR212931")]
Staph_examples_df <-
  Staph_examples_df[rowSums(Staph_examples_df[2:3]) > 0,]

# Species-level info
Staph_examples_df <-
  Staph_examples_df[grepl(x = Staph_examples_df$clade_name,
                          pattern = "s__"), ]
Staph_examples_df <-
  Staph_examples_df[! grepl(x = Staph_examples_df$clade_name,
                            pattern = "t__"), ]
Staph_examples_df$clade_name <-
  str_replace_all(string = Staph_examples_df$clade_name,
                  pattern = ".*s__",
                  replacement = "")
View(Staph_examples_df[order(Staph_examples_df$clade_name), ])

# sessionInfo() ----
sessionInfo()
