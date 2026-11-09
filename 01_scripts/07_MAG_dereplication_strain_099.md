# 7. Enterobacteriaceae strain-level (99% ANI) dereplication

This stage dereplicates all high-quality Enterobacteriaceae MAGs at a
stricter 99% ANI threshold to resolve *Escherichia coli* and related taxa
at the strain level, compares read-based (Sylph) vs. MAG-based
Enterobacteriaceae species composition, and computes strain-level recovery
prevalence (Fig. 4A–C, Supplementary Fig. 3).

Software version: dRep v3.5.0, same clustering options as the 95%
species-level run. Requires `enterobacteriaceae_HQ`, `qc_final_meta`,
and `meta` from the previous two folders.

## Step 1 — collect Enterobacteriaceae HQ MAGs and run dRep at 99% ANI

```bash
set -euo pipefail
cd 06_MAGS/drep
mkdir -p selected_HQ_enterobacteriaceae

while read -r f; do
  if [ -f "../qc_mags/results/seqkit/${f}" ]; then
    cp "../qc_mags/results/seqkit/${f}" selected_HQ_enterobacteriaceae/
  else
    echo "Missing: ${f}"
  fi
done < enterobacteriaceae_HQ_MAGs.txt

find "$PWD/selected_HQ_enterobacteriaceae" -type f -name "*.fasta" | sort > selected_HQ_enterobacteriaceae_paths.txt
wc -l selected_HQ_enterobacteriaceae_paths.txt

module load drep/3.5.0--pyhdfd78af_0

dRep dereplicate drep_HQ_enterobacteriaceae_099 \
  -g selected_HQ_enterobacteriaceae_paths.txt \
  -p 32 \
  -sa 0.99 \
  --genomeInfo genomeInfo_HQ.csv \
  --S_algorithm fastANI \
  --multiround_primary_clustering \
  --primary_chunksize 5000 \
  --greedy_secondary_clustering \
  --run_tertiary_clustering \
  --skip_plots
```

## Step 2 — build MAG-to-99%-representative membership and attach metadata/taxonomy

```r
library(tidyverse)
library(stringr)

drep_dir <- "06_MAGS/drep/drep_HQ_enterobacteriaceae_099/data_tables"
cdb <- read_csv(file.path(drep_dir, "Cdb.csv"), show_col_types = FALSE)
wdb <- read_csv(file.path(drep_dir, "Wdb.csv"), show_col_types = FALSE)
sdb <- read_csv(file.path(drep_dir, "Sdb.csv"), show_col_types = FALSE)

cdb_clean <- cdb %>%
  transmute(MAG_file = genome, MAG_ID = str_remove(genome, "\\.fasta$"),
            drep_primary_cluster = primary_cluster, drep_secondary_cluster = secondary_cluster,
            drep_original_secondary_cluster = original_secondary_cluster, greedy_representative)

wdb_clean <- wdb %>%
  transmute(Representative_file = genome, Representative_MAG_ID = str_remove(genome, "\\.fasta$"),
            drep_secondary_cluster = cluster, representative_score = score)

sdb_clean <- sdb %>%
  transmute(MAG_ID = str_remove(genome, "\\.fasta$"), drep_score = score)

entero_099_membership <- cdb_clean %>%
  left_join(wdb_clean, by = "drep_secondary_cluster") %>%
  left_join(sdb_clean, by = "MAG_ID")

extract_sample_id_from_mag <- function(x) {
  x %>% str_remove("^cleaned_long_") %>% str_replace("_bin_.*$", "")
}

meta <- read_csv("databases/metadata_long.csv", show_col_types = FALSE)

entero_099_membership_meta <- entero_099_membership %>%
  mutate(Sample_ID = extract_sample_id_from_mag(MAG_ID)) %>%
  left_join(
    meta %>% select(Sample_ID, Participant_ID, country, timepoint, treatment_arm,
                     TR_pre, HKW_pre, ASC_pre, ST_pre, Coinfection_pre,
                     TR_cured, HKW_cured, ASC_cured, AE_any, AE_GID) %>% distinct(),
    by = "Sample_ID"
  )

# Taxonomy, from the Enterobacteriaceae HQ MAG table built in Step 8 of the previous folder
entero_tax <- enterobacteriaceae_HQ %>%
  mutate(
    MAG_ID = postqc_genome_name,
    gtdb_family  = str_extract(gtdbtk_classification, "f__[^;]+") %>% str_remove("^f__"),
    gtdb_genus   = str_extract(gtdbtk_classification, "g__[^;]+") %>% str_remove("^g__"),
    gtdb_species = str_extract(gtdbtk_classification, "s__[^;]+") %>% str_remove("^s__")
  ) %>%
  mutate(
    gtdb_genus = if_else(is.na(gtdb_genus) | gtdb_genus == "", "Unclassified_genus", gtdb_genus),
    gtdb_species = if_else(is.na(gtdb_species) | gtdb_species == "", paste0(gtdb_genus, "_sp"), gtdb_species),
    gtdb_species_clean = gtdb_species %>% str_replace_all(" ", "_") %>%
      str_replace_all("[^A-Za-z0-9_]", "_") %>% str_replace_all("_+", "_") %>% str_remove("_$")
  ) %>%
  select(MAG_ID, gtdbtk_classification, gtdb_family, gtdb_genus, gtdb_species, gtdb_species_clean,
         checkm2_postqc_Completeness, checkm2_postqc_Contamination)

entero_099_membership_meta_tax <- entero_099_membership_meta %>%
  left_join(entero_tax, by = "MAG_ID")
```

## Step 3 — strain naming (per species, by cluster size) and the analysis table

```r
rep_species_table <- entero_099_membership_meta_tax %>%
  filter(MAG_ID == Representative_MAG_ID) %>%
  select(Representative_MAG_ID, gtdb_species, gtdb_species_clean, gtdb_genus, gtdb_family, representative_score) %>%
  distinct() %>%
  arrange(gtdb_species_clean, Representative_MAG_ID) %>%
  group_by(gtdb_species_clean) %>%
  mutate(strain_number = row_number(), strain_cluster_name = paste0(gtdb_species_clean, "_strain_", strain_number)) %>%
  ungroup()

entero_099_analysis_table <- entero_099_membership_meta_tax %>%
  left_join(
    rep_species_table %>% select(Representative_MAG_ID, representative_species = gtdb_species,
                                   representative_species_clean = gtdb_species_clean, representative_genus = gtdb_genus,
                                   representative_family = gtdb_family, strain_number, strain_cluster_name),
    by = "Representative_MAG_ID"
  ) %>%
  mutate(
    is_representative = MAG_ID == Representative_MAG_ID,
    country = factor(country),
    timepoint = factor(timepoint, levels = c("pre", "post")),
    treatment_arm = factor(treatment_arm, levels = c("ALB", "FDCx1", "FDCx3"))
  ) %>%
  arrange(strain_cluster_name, MAG_ID)

write_csv(entero_099_analysis_table, "06_MAGS/drep/entero_099_analysis_table_MAGs_metadata_taxonomy_strain_names.csv")
```

This yields **50 strain-level genomic clusters**, of which **31 belong to
*E. coli*** (matching the manuscript Results).

## Step 4 — Fig. 4B: strain cluster bar charts by country (all species, and *E. coli* only)

```r
library(ggsci)

strain_country_summary <- entero_099_analysis_table %>%
  count(strain_cluster_name, representative_species_clean, country, name = "n_MAGs") %>%
  group_by(strain_cluster_name) %>%
  mutate(total_MAGs = sum(n_MAGs)) %>%
  ungroup() %>%
  arrange(desc(total_MAGs), strain_cluster_name)

write_csv(strain_country_summary, "06_MAGS/drep/entero_099_strain_cluster_country_summary_named.csv")

p_entero_099_country <- strain_country_summary %>%
  mutate(
    strain_cluster_name = str_replace_all(strain_cluster_name, '_', ' '),
    strain_cluster_name = str_replace_all(strain_cluster_name, 'strain ', 's'),
    strain_cluster_name = fct_reorder(strain_cluster_name, total_MAGs)
  ) %>%
  ggplot(aes(x = strain_cluster_name, y = n_MAGs, fill = country)) +
  geom_col() +
  coord_flip() +
  scale_fill_bmj() +
  theme_bw() +
  labs(x = "Enterobacteriaceae 99% ANI strain cluster", y = "Number of HQ MAGs assigned to strain cluster", fill = "Country") +
  theme(legend.position = 'top')

ggsave("Figures/entero_099_strain_clusters_by_country_named.png", p_entero_099_country, width = 6, height = 7.5)
ggsave("Figures/entero_099_strain_clusters_by_country_named.pdf", p_entero_099_country, width = 6, height = 7.5)

p_entero_099_country_Ecoli <- strain_country_summary %>%
  filter(representative_species_clean == 'Escherichia_coli') %>%
  mutate(
    strain_cluster_name = str_replace_all(strain_cluster_name, '_', ' '),
    strain_cluster_name = str_replace_all(strain_cluster_name, 'strain ', 's'),
    strain_cluster_name = fct_reorder(strain_cluster_name, total_MAGs)
  ) %>%
  ggplot(aes(x = strain_cluster_name, y = n_MAGs, fill = country)) +
  geom_col() +
  coord_flip() +
  scale_fill_bmj() +
  theme_bw() +
  labs(x = "E. coli 99% ANI strain cluster", y = "Number of HQ MAGs assigned to strain cluster", fill = "Country") +
  theme(legend.position = 'top')

ggsave("Figures/entero_099_coli_strain_clusters_by_country_named.png", p_entero_099_country_Ecoli, width = 4, height = 5)
ggsave("Figures/entero_099_coli_strain_clusters_by_country_named.pdf", p_entero_099_country_Ecoli, width = 4, height = 5)
```

## Step 5 — strain cluster prevalence by country/arm/timepoint, and pre/post change

```r
meta_denominator <- meta %>%
  filter(!is.na(Sample_ID), !is.na(country), !is.na(treatment_arm), !is.na(timepoint)) %>%
  distinct(Sample_ID, country, treatment_arm, timepoint) %>%
  count(country, treatment_arm, timepoint, name = "n_samples")

strain_sample_presence <- entero_099_analysis_table %>%
  filter(!is.na(Sample_ID), !is.na(country), !is.na(treatment_arm), !is.na(timepoint), !is.na(strain_cluster_name)) %>%
  distinct(strain_cluster_name, representative_species_clean, country, treatment_arm, timepoint, Sample_ID) %>%
  count(strain_cluster_name, representative_species_clean, country, treatment_arm, timepoint, name = "n_positive_samples")

strain_prevalence_summary <- strain_sample_presence %>%
  left_join(meta_denominator, by = c("country", "treatment_arm", "timepoint")) %>%
  mutate(prevalence = n_positive_samples / n_samples) %>%
  group_by(strain_cluster_name) %>%
  mutate(max_prevalence = max(prevalence, na.rm = TRUE), total_positive_samples = sum(n_positive_samples, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(max_prevalence), strain_cluster_name)

write_csv(strain_prevalence_summary, "06_MAGS/drep/entero_099_strain_cluster_prevalence_by_group.csv")

strain_prepost_prevalence <- strain_prevalence_summary %>%
  select(strain_cluster_name, representative_species_clean, country, treatment_arm, timepoint, n_positive_samples, n_samples, prevalence) %>%
  pivot_wider(names_from = timepoint, values_from = c(n_positive_samples, n_samples, prevalence), values_fill = 0) %>%
  mutate(
    prevalence_pre = replace_na(prevalence_pre, 0),
    prevalence_post = replace_na(prevalence_post, 0),
    prevalence_delta_post_minus_pre = prevalence_post - prevalence_pre,
    prevalence_ratio_post_pre = (prevalence_post + 0.005) / (prevalence_pre + 0.005),
    log2_prevalence_ratio_post_pre = log2(prevalence_ratio_post_pre)
  ) %>%
  arrange(log2_prevalence_ratio_post_pre)

write_csv(strain_prepost_prevalence, "06_MAGS/drep/entero_099_strain_cluster_prepost_prevalence_change.csv")
```

## Step 6 — Fig. 4C: pre/post prevalence of the two dominant strains (*E. coli* s14, s29)

```r
library(hrbrthemes)

dominant_strains <- c("Escherichia coli s14", "Escherichia coli s29")
strain_prevalence_summary$strain_cluster_name <- str_replace(strain_prevalence_summary$strain_cluster_name, "Escherichia_coli_strain_29", "Escherichia coli s29")
strain_prevalence_summary$strain_cluster_name <- str_replace(strain_prevalence_summary$strain_cluster_name, "Escherichia_coli_strain_14", "Escherichia coli s14")

p_dominant_entero_strains <- strain_prevalence_summary %>%
  filter(strain_cluster_name %in% dominant_strains) %>%
  mutate(
    timepoint = factor(timepoint, levels = c("pre", "post"), labels = c("Pre-treatment", "Post-treatment")),
    strain_cluster_name = factor(strain_cluster_name, levels = dominant_strains)
  ) %>%
  ggplot(aes(x = timepoint, y = prevalence, fill = timepoint)) +
  geom_col(color = "black", position = 'dodge') +
  facet_grid(country + treatment_arm ~ strain_cluster_name) +
  scale_y_percent() + ylim(0, 0.3) +
  scale_fill_futurama() +
  theme_bw() +
  labs(x = NULL, y = "Prevalence") +
  theme(legend.position = "none", strip.text.x = element_text(face = "italic"), axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Figures/entero_099_dominant_strains_prepost.png", p_dominant_entero_strains, width = 4, height = 6)
ggsave("Figures/entero_099_dominant_strains_prepost.pdf", p_dominant_entero_strains, width = 8, height = 10)
```

## Step 7 — Fig. 4A: read-based (Sylph) vs. MAG-based Enterobacteriaceae species composition

Requires the per-sample Sylph/GlobDB clade-abundance files from module 2.

```r
library(scales)

read_sylph <- function(file) {
  dat <- read.delim(file, skip = 1, stringsAsFactors = FALSE)
  sample_name <- basename(file) %>% str_remove("^tax_") %>% str_remove("_qc.*$")
  dat %>% select(clade_name, relative_abundance) %>% mutate(Sample_ID = sample_name)
}

sylph_files <- list.files("03_TAXONOMY/tax_files_GlobDB/", pattern = "\\.sylphmpa$", full.names = TRUE)
sylph_df <- map_dfr(sylph_files, read_sylph)

sylph_tax_long <- sylph_df %>%
  separate(clade_name, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"),
           sep = "\\|", fill = "right", remove = FALSE) %>%
  mutate(
    across(c(Kingdom, Phylum, Class, Order, Family, Genus, Species, Strain), ~ str_remove(.x, "^[a-z]__")),
    across(c(Kingdom, Phylum, Class, Order, Family, Genus, Species, Strain), ~ na_if(.x, ""))
  )

sylph_entero_species_long <- sylph_tax_long %>%
  filter(Family == "Enterobacteriaceae", !is.na(Species), is.na(Strain)) %>%
  left_join(meta %>% distinct(Sample_ID, .keep_all = TRUE), by = "Sample_ID") %>%
  mutate(
    Species_clean = str_replace_all(Species, "_", " "),
    timepoint = factor(timepoint, levels = c("pre", "post"), labels = c("Pre-treatment", "Post-treatment")),
    treatment_arm = factor(treatment_arm, levels = c("ALB", "FDCx1", "FDCx3"))
  )

n_total_samples <- meta %>% distinct(Sample_ID) %>% nrow()

# Keep species detected at >=0.01% relative abundance in >=5% of samples; everything
# else is pooled into "Other Enterobacteriaceae"
sylph_entero_species_keep <- sylph_entero_species_long %>%
  filter(relative_abundance >= 0.01) %>%
  distinct(Sample_ID, Species_clean) %>%
  count(Species_clean, name = "n_positive_samples") %>%
  mutate(prevalence = n_positive_samples / n_total_samples) %>%
  filter(prevalence >= 0.05) %>%
  pull(Species_clean)

sample_meta <- meta %>%
  distinct(Sample_ID, .keep_all = TRUE) %>%
  select(Sample_ID, country, treatment_arm, timepoint) %>%
  mutate(
    timepoint = factor(timepoint, levels = c("pre", "post"), labels = c("Pre-treatment", "Post-treatment")),
    treatment_arm = factor(treatment_arm, levels = c("ALB", "FDCx1", "FDCx3"))
  )

entero_total_by_sample <- sylph_entero_species_long %>%
  group_by(Sample_ID) %>% summarise(total_entero_abundance = sum(relative_abundance, na.rm = TRUE), .groups = "drop")

kept_species_by_sample <- sylph_entero_species_long %>%
  filter(Species_clean %in% sylph_entero_species_keep) %>%
  group_by(Sample_ID, Species_clean) %>%
  summarise(relative_abundance = sum(relative_abundance, na.rm = TRUE), .groups = "drop")

kept_species_complete <- sample_meta %>%
  select(Sample_ID) %>%
  crossing(Species_clean = sylph_entero_species_keep) %>%
  left_join(kept_species_by_sample, by = c("Sample_ID", "Species_clean")) %>%
  mutate(relative_abundance = replace_na(relative_abundance, 0))

other_entero_by_sample <- sample_meta %>%
  select(Sample_ID) %>%
  left_join(entero_total_by_sample, by = "Sample_ID") %>%
  left_join(kept_species_complete %>% group_by(Sample_ID) %>% summarise(kept_entero_abundance = sum(relative_abundance, na.rm = TRUE), .groups = "drop"), by = "Sample_ID") %>%
  mutate(
    total_entero_abundance = replace_na(total_entero_abundance, 0),
    kept_entero_abundance = replace_na(kept_entero_abundance, 0),
    relative_abundance = pmax(total_entero_abundance - kept_entero_abundance, 0),
    Species_clean = "Other Enterobacteriaceae"
  ) %>%
  select(Sample_ID, Species_clean, relative_abundance)

sylph_entero_species_long_comp <- bind_rows(kept_species_complete, other_entero_by_sample) %>%
  left_join(sample_meta, by = "Sample_ID")

# Shared Enterobacteriaceae species colour palette (used for both panels of Fig. 4A)
entero_species_palette <- c(
  "Escherichia coli" = "#0072B5", "Escherichia coli F" = "#56B4E9", "Escherichia SPIREOTU 00968429" = "#4E79A7",
  "Klebsiella pneumoniae" = "#E69F00", "Klebsiella variicola" = "#F0E442", "Klebsiella quasipneumoniae" = "#D55E00",
  "Klebsiella quasivariicola" = "#CC79A7", "Klebsiella michiganensis" = "#B07AA1", "Klebsiella FMDMAG 09069" = "#F28E2B",
  "Enterobacter bugandensis" = "#009E73", "Enterobacter cloacae L" = "#59A14F", "Enterobacter cloacae" = "#39A14F",
  "Enterobacter cloacae M" = "#8CD17D", "Enterobacter hormaechei C" = "#00A087", "Enterobacter huaxiensis" = "springgreen2",
  "Enterobacter roggenkampii" = "#76B7B2", "Enterobacter dissolvens" = "#86BCB6", "Enterobacter kobei" = "#499894",
  "Cronobacter sakazakii" = "#7B68EE", "Pseudocitrobacter faecalis" = "#7B88EF", "Franconibacter helveticus" = "#7B38EF",
  "Other Enterobacteriaceae" = "#BDBDBD"
)

sylph_entero_species_group <- sylph_entero_species_long_comp %>%
  group_by(Species_clean, country, treatment_arm, timepoint) %>%
  summarise(mean_abundance = mean(relative_abundance, na.rm = TRUE), .groups = "drop")

p_sylph_entero_species_bar <- ggplot(sylph_entero_species_group, aes(x = timepoint, y = mean_abundance, fill = Species_clean)) +
  geom_col(color = "black", linewidth = 0.1) +
  scale_fill_manual(values = entero_species_palette, drop = FALSE, name = "Enterobacteriaceae species") +
  facet_grid(country ~ treatment_arm) +
  theme_bw() +
  labs(x = NULL, y = "Mean relative abundance (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "left", legend.title = element_blank())

# MAG-based recovery prevalence (summed across all species detected, labelled "MAG recovery summed prevalence")
mag_species_presence <- entero_099_analysis_table %>%
  filter(!is.na(Sample_ID), !is.na(country), !is.na(treatment_arm), !is.na(timepoint), !is.na(representative_species_clean)) %>%
  distinct(Sample_ID, country, treatment_arm, timepoint, representative_species_clean) %>%
  count(country, treatment_arm, timepoint, representative_species_clean, name = "n_positive_samples")

mag_denominator <- meta %>%
  filter(!is.na(Sample_ID), !is.na(country), !is.na(treatment_arm), !is.na(timepoint)) %>%
  distinct(Sample_ID, country, treatment_arm, timepoint) %>%
  count(country, treatment_arm, timepoint, name = "n_samples")

mag_species_prev <- mag_species_presence %>%
  left_join(mag_denominator, by = c("country", "treatment_arm", "timepoint")) %>%
  mutate(prevalence = n_positive_samples / n_samples, species_plot = str_replace_all(representative_species_clean, "_", " "))

top_mag_entero_species <- mag_species_prev %>%
  group_by(species_plot) %>%
  summarise(max_prevalence = max(prevalence, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_prevalence)) %>%
  slice_head(n = 10) %>%
  pull(species_plot)

mag_species_prev_plot <- mag_species_prev %>%
  mutate(species_plot = if_else(species_plot %in% top_mag_entero_species, species_plot, "Other Enterobacteriaceae")) %>%
  group_by(country, treatment_arm, timepoint, species_plot) %>%
  summarise(prevalence = sum(prevalence, na.rm = TRUE), n_positive_samples = sum(n_positive_samples, na.rm = TRUE),
            n_samples = max(n_samples, na.rm = TRUE), .groups = "drop") %>%
  mutate(timepoint = factor(timepoint, levels = c("pre", "post"), labels = c("Pre-treatment", "Post-treatment")))

p_mag_entero_species_prev_bar <- ggplot(mag_species_prev_plot, aes(x = timepoint, y = prevalence, fill = species_plot)) +
  geom_col(color = "black", linewidth = 0.1) +
  facet_grid(country ~ treatment_arm) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  scale_fill_manual(values = entero_species_palette, drop = FALSE, name = "Enterobacteriaceae species") +
  labs(x = NULL, y = "MAG recovery summed prevalence", fill = "Enterobacteriaceae species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right", legend.title = element_blank())

entero_comp_sylph_mags <- ggarrange(p_sylph_entero_species_bar, p_mag_entero_species_prev_bar)

ggsave("Figures/entero_099_sylph_mag_species_bar_timepoint.png", entero_comp_sylph_mags, width = 11, height = 4)
ggsave("Figures/entero_099_sylph_mag_species_bar_timepoint.pdf", entero_comp_sylph_mags, width = 11, height = 4)
```

## Step 8 — number of strain clusters per species (Supplementary Fig. 3A)

```r
strains_per_species <- entero_099_analysis_table %>%
  distinct(representative_species_clean, strain_cluster_name) %>%
  count(representative_species_clean, name = "n_strain_clusters") %>%
  arrange(desc(n_strain_clusters))
```

## Step 9 — strain-level phylogeny with country-prevalence heatmap (Supplementary Fig. 3B)

Tree built upstream with GTDB-Tk identify/align (bac120) + FastTree, on the
50 dereplicated Enterobacteriaceae strain-level representative genomes.

```r
library(ape)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)

rep_country_pattern <- entero_099_analysis_table %>%
  distinct(Representative_MAG_ID, strain_cluster_name, representative_species_clean, country) %>%
  filter(!is.na(country)) %>%
  mutate(present = TRUE) %>%
  pivot_wider(names_from = country, values_from = present, values_fill = FALSE) %>%
  mutate(country_pattern = case_when(
    Ethiopia == TRUE & Mozambique == TRUE ~ "Shared",
    Ethiopia == TRUE & Mozambique == FALSE ~ "Ethiopia-only",
    Ethiopia == FALSE & Mozambique == TRUE ~ "Mozambique-only",
    TRUE ~ "Unknown"
  ))

country_denominator <- meta %>%
  filter(!is.na(Sample_ID), !is.na(country)) %>%
  distinct(Sample_ID, country) %>%
  count(country, name = "n_samples_country")

strain_country_presence <- entero_099_analysis_table %>%
  filter(!is.na(Sample_ID), !is.na(country), !is.na(strain_cluster_name)) %>%
  distinct(strain_cluster_name, representative_species_clean, Representative_MAG_ID, country, Sample_ID) %>%
  count(strain_cluster_name, representative_species_clean, Representative_MAG_ID, country, name = "n_positive_samples_country") %>%
  left_join(country_denominator, by = "country") %>%
  mutate(country_prevalence = n_positive_samples_country / n_samples_country)

strain_country_prev_wide <- strain_country_presence %>%
  select(Representative_MAG_ID, strain_cluster_name, representative_species_clean, country, country_prevalence, n_positive_samples_country) %>%
  pivot_wider(names_from = country, values_from = c(country_prevalence, n_positive_samples_country), values_fill = 0)

tree_meta_entero_099 <- rep_species_table %>%
  left_join(rep_country_pattern %>% select(Representative_MAG_ID, country_pattern), by = "Representative_MAG_ID") %>%
  left_join(strain_country_prev_wide, by = c("Representative_MAG_ID", "strain_cluster_name", "gtdb_species_clean" = "representative_species_clean")) %>%
  mutate(
    country_pattern = replace_na(country_pattern, "Unknown"),
    country_prevalence_Ethiopia = replace_na(country_prevalence_Ethiopia, 0),
    country_prevalence_Mozambique = replace_na(country_prevalence_Mozambique, 0),
    species_label = str_replace_all(gtdb_species_clean, "_", " "),
    tree_label = Representative_MAG_ID,
    strain_label = str_replace_all(strain_cluster_name, "_", " ") %>% str_replace(" strain ", " s")
  )

write_csv(tree_meta_entero_099, "06_MAGS/drep/entero_099_tree_metadata_country_prevalence.csv")

tree <- read.tree("06_MAGS/drep/entero_099_gtdbtk_bac120_fasttree.nwk")
tree$tip.label <- tree$tip.label %>% basename() %>% str_remove("\\.fasta$") %>% str_remove("\\.fa$") %>% str_remove("\\.fna$")

tree_meta <- tree_meta_entero_099

heat_data <- tree_meta %>%
  transmute(id = tree_label, Ethiopia = country_prevalence_Ethiopia, Mozambique = country_prevalence_Mozambique) %>%
  pivot_longer(cols = c(Ethiopia, Mozambique), names_to = "country", values_to = "prevalence") %>%
  mutate(country = factor(country, levels = c("Ethiopia", "Mozambique")))

p_tree_heat <- ggtree(tree, layout = "circular", branch.length = "none", size = 0.35) %<+% tree_meta +
  geom_tippoint(aes(color = species_label), size = 2) +
  scale_color_manual(values = entero_species_palette, drop = FALSE, name = "Species") +
  ggnewscale::new_scale_fill() +
  geom_fruit(data = heat_data, geom = geom_tile, mapping = aes(y = id, x = country, fill = prevalence),
             color = 'grey90', offset = 0.15, pwidth = 0.25,
             axis.params = list(axis = "x", text.angle = 50, text.size = 2.5, hjust = 1)) +
  scale_fill_gradient(low = "white", high = "darkgoldenrod2", name = "MAG recovery\nprevalence", labels = scales::percent) +
  theme(legend.position = "right")

ggsave("Figures/entero_099_tree.png", p_tree_heat, width = 9, height = 6)
ggsave("Figures/entero_099_tree.pdf", p_tree_heat, width = 9, height = 6)
```

