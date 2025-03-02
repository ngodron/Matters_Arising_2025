# 16 sept. 2024
# Author: N. Godron (nicolas.godron {a} inserm.fr)
# Visuals for metaphlan interpretation in light of SRA metadata.

# Libraries ----
library(tidyverse)
library(taxizedb)
library(ggbeeswarm)
library(ComplexUpset)

# Parameters ----
# Expressed in % of relative abundance (as in metaphlan output)
cont_threshold <- 5

# Inputs ----
main_dir <- "/home/nicolas/IAME/Ellabaan_Response/"

## metadata ----
metadata_df <-
  read_tsv(paste0(main_dir, "All_run_info.tsv"))

too_many_columns <- problems(metadata_df)
# Upon examination, 6 more columns are found in 65 metadata rows
# The columns are not of interest to this work and are thus discarded.
# They are: library_prep_location,	rna_prep_3_protocol, ph,
#   sequencing_longitude,	tissue_type, and	isolation_source
rm(too_many_columns)

## metaphlan ----
contamination_df <-
  read_tsv(paste0(main_dir,
                  "20250203_All_metaphlan_profiles.tsv"),
           comment = "#")

# Filtering, parsing & querying taxonomy ----
## metadata, GNR ----
metadata_df <-
  metadata_df |>
  select(c("run_accession", "sample_alias", "scientific_name",
           "tax_id", "tax_lineage",
           "sample_accession", "submission_accession", 
           "study_accession", "experiment_accession",
           "library_selection", "fastq_bytes", "description",
           "read_count", "base_count", "library_strategy", "library_layout",
           "tag", "first_created"))

# Renaming columns to allow for comparison with metadata
names(contamination_df) <- str_replace(names(contamination_df),
                                       pattern = "^X",
                                       replacement = "")

names(contamination_df) <- str_replace(names(contamination_df),
                                       pattern = "_metaphlan$",
                                       replacement = "")

# Three samples are in the table twice due to aliases (cf. metadata_df)
duplicates <- data.frame(contamination_df$clade_name,
                         contamination_df$ERR594079, contamination_df$"0788",
                         contamination_df$ERR2204715, contamination_df$"7518",
                         contamination_df$ERR1512641, contamination_df$"394-2_staphylococcus-aureus")
( duplicates <- duplicates[rowSums(duplicates[2:5] != 0) > 0, ] )

# Filtering out the three duplicate samples
contamination_df <- contamination_df |>
  select(! c("0788", "7518", "394-2_staphylococcus-aureus")) |>
  # Renaming a strain with a super long name with its ERR ID instead
  rename(ERR1549309 = `Isolat27_IonXpress_012_R_2012_11_16_09_40_45_user_SN1-88_Auto_user_SN1-88_89`) |>
  identity()


# NCBI taxonomy DB download and querying
## Taxonomy retrieval ----
taxizedb::db_download_ncbi()

# Checking NCBI taxonomic IDs
taxize_NCBI_ID <-
  as.numeric(taxizedb::name2taxid(metadata_df$scientific_name))
table(taxize_NCBI_ID == metadata_df$tax_id)

fam_phyl <-
  function(x) {
    taxonomy <-
      taxizedb::classification(x)[[1]]
    taxonomy[taxonomy$rank %in% c("genus", "family", "phylum"), "name"]
  }

fam_phyl_columns <-
  data.frame(t(sapply(
  metadata_df$tax_id,
  fam_phyl)))

metadata_df$phylum <- fam_phyl_columns[, 1]
metadata_df$family <- fam_phyl_columns[, 2]
metadata_df$genus <- fam_phyl_columns[, 3]

# Special case for resolved name of length 1 word (i.e. only Genus)
# metadata_df$species[is.na(metadata_df$species)] <-
#   metadata_df$resolved_name[is.na(metadata_df$species)]

### Fetching family ----
# The step below can take a minute
unique_metadata_family <- unique(metadata_df$family)

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
                                   between(sra_abundance, 0, 5) ~ "Mismatch",
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
violin_df$SRA_of_interest <-
  violin_df$Name %in% c("SRR1522630",
                        "ERR212931")
violin_df <- violin_df[order(violin_df$family), ]

violin_plot <- ggplot(data = violin_df) +
  aes(x = "",
      y = sra_abundance) +
  ggbeeswarm::geom_beeswarm(aes(x = family,
                                y = sra_abundance,
                                fill = family),
                            shape = 22,
                            colour = "#666666",
                            method = "square",
                            cex = 1.66,
                            side = 1L,
                            size = 1.7) +
  scale_fill_manual(values = c("#ffb3f2", "#cea30d", "#44bebe")) +
  xlab("") +
  ylab("Observed relative abundance (%)\nof the expected (SRA) family") +
  ylim(c(0, 100)) +
  scale_y_continuous(breaks = c(0, 5, 80, 95, 100),
                     minor_breaks = seq(0, 100, 10)) +
  theme_minimal() +
  theme(axis.text = element_text(size = 8,
                                 family = "Arial")) +
  annotate("rect", fill = "#64ab05",
           xmin = 0.6, xmax = 0.8,
           ymin = 95, ymax = 100) +
  annotate("rect", fill = "#f3d819",
           xmin = 0.6, xmax = 0.8,
           ymin = 80, ymax = 95) +
  annotate("rect", fill = "#c86c44",
           xmin = 0.6, xmax = 0.8,
           ymin = 5, ymax = 80) +
  annotate("rect", fill = "#9c2a17",
           xmin = 0.6, xmax = 0.8,
           ymin = 0, ymax = 5) +
  geom_segment(x = 1.75, xend = 1.95,
               y = 94.7, yend = 99.2,
               arrow = arrow(length = unit(.2,"cm"),
                             type = "closed"),
               colour = "black",
               linewidth = 1.3) +
  annotate("text",
            x = 1.5, y = 92,
           label = "SRR1522630",
           colour = "black",
           size = 3.25,
           fontface = 2,
           family = "Arial") +
  geom_segment(x = 1.75, xend = 1.95,
               y = 10, yend = 5.5,
               arrow = arrow(length = unit(.2,"cm"),
                             type = "closed"),
               colour = "black",
               linewidth = 1.3) +
  annotate("text",
           x = 1.5, y = 12.5,
           label = "ERR212931",
           colour = "black",
           size = 3.25,
           fontface = 2,
           family = "Arial")

violin_plot <- violin_plot +
  coord_cartesian(xlim = c(1.2,3.2)) +
  theme(text = element_text(family = "Arial"),
        legend.text = element_text(family = "Arial"),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.275))

# No legend option
violin_plot <- violin_plot + theme(legend.position = "none")


png(filename = paste0(main_dir, "Percentage_matching_abundance.png"),
    width = 2600, height = 2600, res = 600,
    pointsize = 2)
violin_plot
dev.off()

## metaphlan UpSet ----

phylum_UpSet_metadata <- data.frame(
  family_and_phylum[, c("set","phylum")]
)


upset_plot <- ComplexUpset::upset(
  ### Basic input & parameters ----
  data = threshold_family,
  intersect = family_and_phylum$set,
  height_ratio = 0.8,
  # name refers to the main x label (below the intersection grid) 
  name = "Observed family intersections (bottom), and their sample counts (top)",
  # Removing empty intersections
  min_degree = 1,
  
  sort_sets = FALSE,
  sort_intersections_by = c('cardinality', 'degree', 'ratio'),
  
  stripes = upset_stripes(mapping = aes(color=phylum),
                          colors = c("Proteobacteria"="#eeeeeeee",
                                     "Firmicutes"="#ccccccee",
                                     "Actinobacteria"="#aaaaaaee",
                                     "Bacteroidota"="#888888ee"),
                          data = phylum_UpSet_metadata),
  
  ### Intersections ----
  base_annotations = list(
    "Intersection size" = intersection_size(
      text_mapping = aes(
        label = !!get_size_mode("exclusive_intersection"),
        # Shenanigan to place label on top of bar, there must exist a better way.
        y = ifelse(!!get_size_mode("exclusive_intersection"),
                   !!get_size_mode("exclusive_intersection") + .3,
                   0)),
      text_colors = c(on_background = "black",
                      on_bar = "black"),
      text = list(size = 4),
      mapping = aes(fill = SRA_Agreement)) +
      scale_fill_manual(values = c("95-100%"="#64ab05",
                                   "80-95%"="#f3d819",
                                   "5-80%"="#c86c44",
                                   "Mismatch"="#9c2a17")) +
      ylab("") +
      theme(axis.title =
              element_text(size = 14,
                           face = "bold",
                           hjust = 1,
                           family = "Arial")) +
      annotate("text",
               x = 6.35,
               y = 22,
               size = 5,
               label = "SRR1522630",
               colour = "black",
               fontface = 2,
               family = "Arial") +
      geom_segment(x = 4.72, xend = 4.72,
                 y = 21, yend = 14,
                 colour = "black",
                 arrow = arrow(length = unit(.2,"cm"),
                               type = "closed"),
                 linewidth = 1) +
      annotate("text",
               x = 10.2,
               y = 12,
               size = 5,
               label = "ERR212931",
               colour = "black",
               fontface = 2,
               family = "Arial") +
      geom_segment(x = 8.72, xend = 8.72,
                   y = 11, yend = 4,
                   colour = "black",
                   arrow = arrow(length = unit(.2,"cm"),
                                 type = "closed"),
                   linewidth = 1)
),
  
  ### Graphical adjustments ----
  set_sizes=FALSE,
  themes = upset_modify_themes(to_update = list(
    "intersections_matrix" = theme(text = element_text(size = 15)))
  )
)
upset_plot <- upset_plot + theme(text = element_text(family = "Arial"))


png(filename = paste0(main_dir, "UpSet_plot.png"),
    width = 6800, height = 4000, res = 600)
upset_plot
dev.off()

upset_plot <- ComplexUpset::upset(
  ### Basic input & parameters ----
  data = threshold_family,
  intersect = family_and_phylum$set,
  height_ratio = 0.8,
  # name refers to the main x label (below the intersection grid) 
  name = "Observed family intersections (bottom), and their sample counts (top)",
  # Removing empty intersections
  min_degree = 1,
  
  sort_sets = FALSE,
  sort_intersections_by = c('cardinality', 'degree', 'ratio'),
  
  stripes = upset_stripes(mapping = aes(color=phylum),
                          colors = c("Proteobacteria"="#eeeeeeee",
                                     "Firmicutes"="#ccccccee",
                                     "Actinobacteria"="#aaaaaaee",
                                     "Bacteroidota"="#888888ee"),
                          data = phylum_UpSet_metadata),
  
  ### Intersections ----
  base_annotations = list(
    "Intersection size" = intersection_size(
      text_mapping = aes(
        label = !!get_size_mode("exclusive_intersection"),
        # Shenanigan to place label on top of bar, there must exist a better way.
        y = ifelse(!!get_size_mode("exclusive_intersection"),
                   !!get_size_mode("exclusive_intersection") + .3,
                   0)),
      text_colors = c(on_background = "black",
                      on_bar = "black"),
      text = list(size = 4),
      mapping = aes(fill = family)) +
      scale_fill_manual(values = c("Enterobacteriaceae"="#ffb3f2",
                                   "Staphylococcaceae"="#cea30d",
                                   "Streptococcaceae"="#44bebe")) +
      ylab("") +
      theme(axis.title =
              element_text(size = 14,
                           face = "bold",
                           hjust = 1,
                           family = "Arial")) +
      annotate("text",
               x = 6.35,
               y = 22,
               size = 5,
               label = "SRR1522630",
               colour = "black",
               fontface = 2,
               family = "Arial") +
      geom_segment(x = 4.72, xend = 4.72,
                   y = 21, yend = 8,
                   colour = "black",
                   arrow = arrow(length = unit(.2,"cm"),
                                 type = "closed"),
                   linewidth = 1) +
      annotate("text",
               x = 10.2,
               y = 12,
               size = 5,
               label = "ERR212931",
               colour = "black",
               fontface = 2,
               family = "Arial") +
      geom_segment(x = 8.72, xend = 8.72,
                   y = 11, yend = 4,
                   colour = "black",
                   arrow = arrow(length = unit(.2,"cm"),
                                 type = "closed"),
                   linewidth = 1)
  ),
  
  ### Graphical adjustments ----
  set_sizes=FALSE,
  themes = upset_modify_themes(to_update = list(
    "intersections_matrix" = theme(text = element_text(size = 15)))
  )
)
upset_plot <- upset_plot + theme(text = element_text(family = "Arial"))

upset_plot

png(filename = paste0(main_dir, "UpSet_plot_alt.png"),
    width = 6800, height = 4000, res = 600)
upset_plot
dev.off()

# Saving tables
write_tsv(threshold_family,
          paste0(main_dir,
                 "Threshold_family.tsv"))

write_tsv(violin_df,
          paste0(main_dir,
                 "Abundance_phylum_family.tsv"))

write_tsv(metaphlan_sra_ref,
          paste0(main_dir,
                 "Distribution_table.tsv"))

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
# View(Staph_examples_df[order(Staph_examples_df$clade_name), ])

# sessionInfo() ----
writeLines(capture.output(sessionInfo(),
                          paste0(main_dir,
                                 "parsing_RsessionInfo.txt")))

