# 6. Species-level (95% ANI) MAG dereplication

This stage dereplicates the 22,032 high-quality MAGs into species-level
representative genomes, classifies their taxonomic novelty, and builds the
bacterial and archaeal phylogenies with country/prevalence/cluster-size
annotation rings (Fig. 4A right panel, Extended Data Fig. 8–9).

Software version: dRep v3.5.0, using fastANI as the secondary-clustering
ANI algorithm, multi-round primary clustering, greedy secondary clustering,
and tertiary clustering (Methods). Requires `qc_final_meta`, `hq_mags`, and
`genomeInfo_HQ.csv` from the previous folder.

## Step 1 — collect HQ MAG fastas and run dRep at 95% ANI

```bash
set -euo pipefail
cd 06_MAGS/drep

while read -r f; do
  if [ -f "../qc_mags/results/seqkit/${f}" ]; then
    cp "../qc_mags/results/seqkit/${f}" selected_HQ_mags/
  else
    echo "Missing: ${f}"
  fi
done < ../qc_mags/selected_HQ_mags.txt

find "$PWD/selected_HQ_mags" -type f -name "*.fasta" | sort > selected_HQ_mags_paths.txt
wc -l selected_HQ_mags_paths.txt

module load drep/3.5.0--pyhdfd78af_0

dRep dereplicate drep_HQ_095 \
  -g selected_HQ_mags_paths.txt \
  -p 32 \
  -sa 0.95 \
  --genomeInfo genomeInfo_HQ.csv \
  --S_algorithm fastANI \
  --multiround_primary_clustering \
  --primary_chunksize 5000 \
  --greedy_secondary_clustering \
  --run_tertiary_clustering \
  --skip_plots

echo "Representative genome count:"
ls drep_HQ_095/dereplicated_genomes/*.fasta | wc -l
```

dRep writes `Cdb.csv` (every genome's cluster assignment), `Wdb.csv` (the
winning representative per cluster), `Sdb.csv` (per-genome dRep score),
`Bdb.csv` (genome path/length), and `Widb.csv` (cluster membership
counts/distances) — all five tables are used below.

## Step 2 — build the full MAG-to-representative membership table (R)

```r
library(tidyverse)
library(stringr)

drep_095_dir <- "06_MAGS/drep/drep_HQ_095/data_tables"

cdb_095 <- read_csv(file.path(drep_095_dir, "Cdb.csv"), show_col_types = FALSE)
wdb_095 <- read_csv(file.path(drep_095_dir, "Wdb.csv"), show_col_types = FALSE)
sdb_095 <- read_csv(file.path(drep_095_dir, "Sdb.csv"), show_col_types = FALSE)
bdb_095 <- read_csv(file.path(drep_095_dir, "Bdb.csv"), show_col_types = FALSE)
widb_095 <- read_csv(file.path(drep_095_dir, "Widb.csv"), show_col_types = FALSE)

cdb_095_clean <- cdb_095 %>%
  transmute(
    MAG_file = genome,
    MAG_ID = str_remove(genome, "\\.fasta$"),
    drep_primary_cluster_095 = primary_cluster,
    drep_secondary_cluster_095 = secondary_cluster,
    drep_original_secondary_cluster_095 = original_secondary_cluster,
    drep_greedy_representative_095 = greedy_representative
  )

wdb_095_clean <- wdb_095 %>%
  transmute(
    Representative_file_095 = genome,
    Representative_MAG_ID_095 = str_remove(genome, "\\.fasta$"),
    drep_secondary_cluster_095 = cluster,
    representative_score_095 = score
  )

sdb_095_clean <- sdb_095 %>%
  transmute(MAG_ID = str_remove(genome, "\\.fasta$"), drep_score_095 = score)

bdb_095_clean <- bdb_095 %>%
  transmute(MAG_ID = str_remove(genome, "\\.fasta$"), genome_path = location, genome_length = length)

widb_095_clean <- widb_095 %>%
  transmute(
    Representative_MAG_ID_095 = str_remove(genome, "\\.fasta$"),
    drep_secondary_cluster_095 = cluster,
    cluster_members_095 = cluster_members,
    closest_cluster_member_095 = closest_cluster_member,
    furthest_cluster_member_095 = furthest_cluster_member
  )

hq_095_membership <- cdb_095_clean %>%
  left_join(wdb_095_clean, by = "drep_secondary_cluster_095") %>%
  left_join(sdb_095_clean, by = "MAG_ID") %>%
  left_join(bdb_095_clean, by = "MAG_ID") %>%
  left_join(widb_095_clean, by = c("Representative_MAG_ID_095", "drep_secondary_cluster_095"))

extract_sample_id_from_mag <- function(x) {
  x %>% str_remove("^cleaned_long_") %>% str_replace("_bin_.*$", "")
}

meta <- read_csv("databases/metadata_long.csv", show_col_types = FALSE)

hq_mag_taxonomy <- hq_mags %>%
  mutate(
    MAG_ID = postqc_genome_name,
    gtdb_domain  = str_extract(gtdbtk_classification, "d__[^;]+") %>% str_remove("^d__"),
    gtdb_phylum  = str_extract(gtdbtk_classification, "p__[^;]+") %>% str_remove("^p__"),
    gtdb_class   = str_extract(gtdbtk_classification, "c__[^;]+") %>% str_remove("^c__"),
    gtdb_order   = str_extract(gtdbtk_classification, "o__[^;]+") %>% str_remove("^o__"),
    gtdb_family  = str_extract(gtdbtk_classification, "f__[^;]+") %>% str_remove("^f__"),
    gtdb_genus   = str_extract(gtdbtk_classification, "g__[^;]+") %>% str_remove("^g__"),
    gtdb_species = str_extract(gtdbtk_classification, "s__[^;]+") %>% str_remove("^s__")
  ) %>%
  select(MAG_ID, gtdbtk_classification, gtdb_domain, gtdb_phylum, gtdb_class, gtdb_order,
         gtdb_family, gtdb_genus, gtdb_species, checkm2_postqc_Completeness, checkm2_postqc_Contamination) %>%
  distinct(MAG_ID, .keep_all = TRUE)

hq_095_membership_meta_tax <- hq_095_membership %>%
  mutate(Sample_ID = extract_sample_id_from_mag(MAG_ID)) %>%
  left_join(meta, by = "Sample_ID") %>%
  left_join(hq_mag_taxonomy, by = "MAG_ID")

hq_095_membership_meta_tax %>%
  summarise(
    n_MAGs = n(),
    n_representatives = n_distinct(Representative_MAG_ID_095),
    n_missing_country = sum(is.na(country)),
    n_missing_taxonomy = sum(is.na(gtdbtk_classification))
  )
```

## Step 3 — representative-level taxonomy and taxonomic-assignment depth

```r
rep_095_taxonomy <- hq_095_membership_meta_tax %>%
  filter(MAG_ID == Representative_MAG_ID_095) %>%
  transmute(
    Representative_MAG_ID_095, drep_secondary_cluster_095, representative_score_095, cluster_members_095,
    representative_gtdbtk_classification = gtdbtk_classification,
    representative_phylum = gtdb_phylum, representative_class = gtdb_class,
    representative_order = gtdb_order, representative_family = gtdb_family,
    representative_genus = gtdb_genus, representative_species = gtdb_species,
    representative_species_clean = gtdb_species %>%
      str_replace_all(" ", "_") %>% str_replace_all("[^A-Za-z0-9_]", "_") %>%
      str_replace_all("_+", "_") %>% str_remove("_$"),
    representative_completeness = checkm2_postqc_Completeness,
    representative_contamination = checkm2_postqc_Contamination
  ) %>%
  distinct(Representative_MAG_ID_095, .keep_all = TRUE) %>%
  mutate(
    taxonomic_assignment_depth = case_when(
      !is.na(representative_species) & representative_species != "" ~ "Species assigned",
      !is.na(representative_genus) & representative_genus != "" ~ "Genus assigned, species missing",
      !is.na(representative_family) & representative_family != "" ~ "Family assigned, genus missing",
      !is.na(representative_order) & representative_order != "" ~ "Order assigned, family missing",
      !is.na(representative_class) & representative_class != "" ~ "Class assigned, order missing",
      !is.na(representative_phylum) & representative_phylum != "" ~ "Phylum assigned, class missing",
      TRUE ~ "Unclassified"
    )
  )
```

## Step 4 — country pattern and MAG recovery prevalence per representative

```r
sample_denominator <- meta %>%
  filter(!is.na(Sample_ID), !is.na(country), !is.na(timepoint), !is.na(treatment_arm)) %>%
  distinct(Sample_ID, country, timepoint, treatment_arm) %>%
  count(country, timepoint, treatment_arm, name = "n_samples")

country_denominator <- meta %>%
  filter(!is.na(Sample_ID), !is.na(country)) %>%
  distinct(Sample_ID, country) %>%
  count(country, name = "n_samples_country")

rep_095_country_presence <- hq_095_membership_meta_tax %>%
  filter(!is.na(Sample_ID), !is.na(country), !is.na(Representative_MAG_ID_095)) %>%
  distinct(Representative_MAG_ID_095, drep_secondary_cluster_095, Sample_ID, country)

rep_095_country_prev_wide <- rep_095_country_presence %>%
  count(Representative_MAG_ID_095, drep_secondary_cluster_095, country, name = "n_positive_samples_country") %>%
  left_join(country_denominator, by = "country") %>%
  mutate(prevalence_country = n_positive_samples_country / n_samples_country) %>%
  select(Representative_MAG_ID_095, drep_secondary_cluster_095, country, prevalence_country, n_positive_samples_country) %>%
  pivot_wider(names_from = country, values_from = c(prevalence_country, n_positive_samples_country), values_fill = 0)

rep_095_country_pattern <- rep_095_country_presence %>%
  distinct(Representative_MAG_ID_095, drep_secondary_cluster_095, country) %>%
  mutate(present = TRUE) %>%
  pivot_wider(names_from = country, values_from = present, values_fill = FALSE) %>%
  mutate(country_pattern = case_when(
    Ethiopia == TRUE & Mozambique == TRUE ~ "Shared",
    Ethiopia == TRUE & Mozambique == FALSE ~ "Ethiopia-only",
    Ethiopia == FALSE & Mozambique == TRUE ~ "Mozambique-only",
    TRUE ~ "Unknown"
  ))

rep_095_tree_meta_arch <- rep_095_taxonomy %>%
  left_join(rep_095_country_pattern, by = c("Representative_MAG_ID_095", "drep_secondary_cluster_095")) %>%
  left_join(rep_095_country_prev_wide, by = c("Representative_MAG_ID_095", "drep_secondary_cluster_095")) %>%
  mutate(
    country_pattern = replace_na(country_pattern, "Unknown"),
    prevalence_country_Ethiopia = replace_na(prevalence_country_Ethiopia, 0),
    prevalence_country_Mozambique = replace_na(prevalence_country_Mozambique, 0),
    n_positive_samples_country_Ethiopia = replace_na(n_positive_samples_country_Ethiopia, 0),
    n_positive_samples_country_Mozambique = replace_na(n_positive_samples_country_Mozambique, 0)
  )
```

## Step 5 — manual archaeal taxonomy patch

GTDB-Tk's bacterial (bac120) marker set does not classify archaea; the six
archaeal representative genomes were separately identified and classified
via GTDB-Tk's ANI/ar53 workflow, and are patched in here by MAG ID.

```r
archaeal_taxonomy_manual <- tibble::tribble(
  ~Representative_MAG_ID_095, ~manual_domain, ~manual_phylum, ~manual_class, ~manual_order, ~manual_family, ~manual_genus, ~manual_species, ~archaea_status,

  "cleaned_long_2031469.9_bin_71_orgin_filtered_kept_contigs",
  "Archaea", "Methanobacteriota", "Methanobacteria", "Methanobacteriales",
  "Methanobacteriaceae", "Methanocatella", "Methanocatella smithii", "GTDB-Tk ANI classified",

  "cleaned_long_2078735.6_bin_63_orgin_filtered_kept_contigs",
  "Archaea", "Methanobacteriota", "Methanobacteria", "Methanobacteriales",
  "Methanobacteriaceae", "Methanosphaera", "Methanosphaera cuniculi", "GTDB-Tk ANI classified",

  "cleaned_long_ET_0059_V0_bin_53_orgin_filtered_kept_contigs",
  "Archaea", "Methanobacteriota", "Methanobacteria", "Methanobacteriales",
  "Methanobacteriaceae", "Methanosphaera", "Methanosphaera stadtmanae", "GTDB-Tk ANI classified",

  "cleaned_long_ET_0061_V0_bin_66_orgin_filtered_kept_contigs",
  "Archaea", "Methanobacteriota", "Methanobacteria", "Methanobacteriales",
  "Methanobacteriaceae", "Methanocatella", "Methanocatella sp030539875", "GTDB-Tk ANI classified",

  "cleaned_long_2078710.3_bin_32_permissive_filtered_kept_contigs",
  "Archaea", "Unclassified Archaea", NA_character_, NA_character_,
  NA_character_, NA_character_, NA_character_, "Archaeal marker-positive, GTDB unresolved",

  "cleaned_long_ET_0477_V0_bin_42_strict_filtered_kept_contigs",
  "Archaea", "Unclassified Archaea", NA_character_, NA_character_,
  NA_character_, NA_character_, NA_character_, "Archaeal marker-positive, GTDB unresolved"
)

rep_095_tree_meta_arch <- rep_095_tree_meta_arch %>%
  left_join(archaeal_taxonomy_manual, by = "Representative_MAG_ID_095") %>%
  mutate(
    representative_domain = if_else(!is.na(manual_domain), manual_domain, "Bacteria"),
    representative_phylum = if_else(!is.na(manual_phylum), manual_phylum, representative_phylum),
    representative_class  = if_else(!is.na(manual_class), manual_class, representative_class),
    representative_order  = if_else(!is.na(manual_order), manual_order, representative_order),
    representative_family = if_else(!is.na(manual_family), manual_family, representative_family),
    representative_genus  = if_else(!is.na(manual_genus), manual_genus, representative_genus),
    representative_species = if_else(!is.na(manual_species), manual_species, representative_species),
    representative_species_clean = if_else(
      !is.na(manual_species),
      manual_species %>% str_replace_all(" ", "_") %>% str_replace_all("[^A-Za-z0-9_]", "_") %>%
        str_replace_all("_+", "_") %>% str_remove("_$"),
      representative_species_clean
    )
  ) %>%
  select(-manual_domain, -manual_phylum, -manual_class, -manual_order, -manual_family, -manual_genus, -manual_species)

# Updated, domain-aware taxonomic-assignment depth and novelty classification
rep_095_tree_meta_arch <- rep_095_tree_meta_arch %>%
  mutate(
    taxonomic_assignment_depth_updated = case_when(
      representative_domain == "Archaea" & archaea_status == "GTDB-Tk ANI classified" ~ "Archaeal species assigned",
      representative_domain == "Archaea" & archaea_status == "Archaeal marker-positive, GTDB unresolved" ~ "Archaeal unresolved",
      representative_domain == "Bacteria" & !is.na(representative_species) & representative_species != "" ~ "Bacterial species assigned",
      representative_domain == "Bacteria" & !is.na(representative_genus) & representative_genus != "" ~ "Bacterial genus assigned, species missing",
      representative_domain == "Bacteria" & !is.na(representative_family) & representative_family != "" ~ "Bacterial family assigned, genus missing",
      representative_domain == "Bacteria" & !is.na(representative_order) & representative_order != "" ~ "Bacterial order assigned, family missing",
      representative_domain == "Bacteria" & !is.na(representative_class) & representative_class != "" ~ "Bacterial class assigned, order missing",
      representative_domain == "Bacteria" & !is.na(representative_phylum) & representative_phylum != "" ~ "Bacterial phylum assigned, class missing",
      TRUE ~ "Unclassified"
    ),
    novelty_category_updated = case_when(
      taxonomic_assignment_depth_updated == "Bacterial species assigned" ~ "GTDB species-assigned bacteria",
      taxonomic_assignment_depth_updated == "Bacterial genus assigned, species missing" ~ "Putative novel bacterial species",
      taxonomic_assignment_depth_updated == "Bacterial family assigned, genus missing" ~ "Putative novel bacterial genus",
      taxonomic_assignment_depth_updated %in% c("Bacterial order assigned, family missing", "Bacterial class assigned, order missing", "Bacterial phylum assigned, class missing") ~ "Deeply under-classified bacteria",
      taxonomic_assignment_depth_updated == "Archaeal species assigned" ~ "GTDB species-assigned Archaea",
      taxonomic_assignment_depth_updated == "Archaeal unresolved" ~ "Unclassified Archaea",
      TRUE ~ "Unclassified"
    ),
    novelty_category_updated = factor(novelty_category_updated, levels = c(
      "GTDB species-assigned bacteria", "Putative novel bacterial species", "Putative novel bacterial genus",
      "Deeply under-classified bacteria", "GTDB species-assigned Archaea", "Unclassified Archaea", "Unclassified"
    ))
  )
```

## Step 6 — country-pattern / cluster-size / phylum summary bars (part of Extended Data Fig. 8)

```r
library(ggsci)
library(ggpubr)

country_pattern_cols <- c("Shared" = "#009E73", "Ethiopia-only" = "#0072B5", "Mozambique-only" = "#E69F00", "Unknown" = "grey70")

p_country_pattern_bar <- rep_095_tree_meta_arch %>%
  count(country_pattern) %>%
  mutate(country_pattern = factor(country_pattern, levels = c("Shared", "Ethiopia-only", "Mozambique-only", "Unknown"))) %>%
  ggplot(aes(x = country_pattern, y = n, fill = country_pattern)) +
  geom_col(color = "black") +
  geom_text(aes(label = n), vjust = -1, size = 4) +
  scale_fill_manual(values = country_pattern_cols) +
  ylim(0, 600) +
  theme_bw() +
  labs(x = NULL, y = "Number of representative MAGs") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

p_country_pattern_box <- rep_095_tree_meta_arch %>%
  mutate(country_pattern = factor(country_pattern, levels = c("Shared", "Ethiopia-only", "Mozambique-only", "Unknown"))) %>%
  ggplot(aes(x = country_pattern, y = cluster_members_095, fill = country_pattern)) +
  geom_boxplot() +
  scale_y_log10() +
  scale_fill_manual(values = country_pattern_cols) +
  theme_bw() +
  labs(x = NULL, y = "Cluster size (n) - Log scale") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

rep_phylum_country <- rep_095_tree_meta_arch %>%
  drop_na(representative_phylum) %>%
  mutate(country_pattern = factor(country_pattern, levels = c("Shared", "Ethiopia-only", "Mozambique-only", "Unknown"))) %>%
  count(representative_phylum, country_pattern, name = "n") %>%
  group_by(representative_phylum) %>%
  mutate(total_phylum = sum(n)) %>%
  ungroup()

p_rep_phylum_country_stacked <- rep_phylum_country %>%
  mutate(representative_phylum = fct_reorder(representative_phylum, total_phylum)) %>%
  ggplot(aes(x = representative_phylum, y = n, fill = country_pattern)) +
  geom_col(color = "black", linewidth = 0.15) +
  coord_flip() +
  scale_fill_manual(values = country_pattern_cols) +
  theme_bw() +
  labs(x = "GTDB phylum", y = "Number of representative MAGs") +
  theme(legend.position = "none")

down <- ggarrange(p_country_pattern_bar, p_country_pattern_box, p_rep_phylum_country_stacked, widths = c(1, 1, 2), nrow = 1)
```

## Step 7 — bacterial phylogeny with phylum/country/cluster-size rings (Fig. 4A, Extended Data Fig. 8)

Tree built upstream with GTDB-Tk identify/align + FastTree (bac120 marker
set, Whelan–Goldman model).

```r
library(ape)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(pals)

rep_tree <- read.tree("06_MAGS/drep/drep_HQ_095/rep_095_gtdbtk_bac120_fasttree.nwk")
rep_tree$tip.label <- rep_tree$tip.label %>% basename() %>%
  str_remove("\\.fasta$") %>% str_remove("\\.fa$") %>% str_remove("\\.fna$")

tree_meta <- rep_095_tree_meta_arch %>%
  mutate(
    tree_label = Representative_MAG_ID_095,
    representative_phylum = replace_na(representative_phylum, "Unclassified"),
    country_pattern = factor(replace_na(country_pattern, "Unknown"), levels = c("Shared", "Ethiopia-only", "Mozambique-only", "Unknown")),
    prevalence_country_Mozambique = replace_na(prevalence_country_Mozambique, 0),
    prevalence_country_Ethiopia = replace_na(prevalence_country_Ethiopia, 0)
  ) %>%
  select(tree_label, Representative_MAG_ID_095, representative_phylum, country_pattern,
         prevalence_country_Mozambique, prevalence_country_Ethiopia,
         putative_novelty_status = novelty_category_updated, cluster_members_095)

phylum_cols <- pals::glasbey(14)

phylum_ring <- tree_meta %>% transmute(id = tree_label, ring = "Phylum", representative_phylum)
country_ring <- tree_meta %>% transmute(id = tree_label, ring = "Country pattern", country_pattern)
cluster_ring <- tree_meta %>% transmute(id = tree_label, ring = "Cluster size (n)", cluster_members_095)

p_rep095_tree <- ggtree(rep_tree, layout = "circular", size = 0.12) + theme_void()

p_rep095_tree_phy <- p_rep095_tree +
  geom_fruit(data = phylum_ring, geom = geom_tile, mapping = aes(y = id, x = ring, fill = representative_phylum),
             offset = 0.02, pwidth = 0.06, color = NA) +
  scale_fill_manual(values = phylum_cols, name = "Phylum", guide = guide_legend(ncol = 2, byrow = TRUE))

p_rep095_tree_phy_cont <- p_rep095_tree_phy +
  ggnewscale::new_scale_fill() +
  geom_fruit(data = country_ring, geom = geom_tile, mapping = aes(y = id, x = ring, fill = country_pattern),
             offset = 0.028, pwidth = 0.06, color = NA) +
  scale_fill_manual(values = country_pattern_cols, name = "Country pattern")

p_rep095_tree_size <- p_rep095_tree_phy_cont +
  ggnewscale::new_scale_fill() +
  geom_fruit(data = cluster_ring, geom = geom_tile, mapping = aes(y = id, x = ring, fill = cluster_members_095),
             offset = 0.030, pwidth = 0.06, color = NA) +
  scale_fill_gradientn(colors = c("snow1", "goldenrod1", "#FF4D6D"), name = "MAG cluster size")

rep_tree_plot_arranged <- ggarrange(p_rep095_tree_size, plot.new(), down, ncol = 1, heights = c(1.3, 0.05, 0.6))

ggsave("Figures/rep_095_tree_clusters.png", rep_tree_plot_arranged, width = 10, height = 9)
ggsave("Figures/rep_095_tree_clusters.pdf", rep_tree_plot_arranged, width = 10, height = 9)
```

## Step 8 — archaeal phylogeny with cluster-size/country-prevalence rings (Extended Data Fig. 8)

Only the six archaeal representatives (identified via GTDB-Tk's ANI
workflow) are plotted, on a separate ar53-marker-set tree.

```r
archaea_meta_plot <- rep_095_tree_meta_arch %>%
  filter(representative_domain == "Archaea") %>%
  mutate(
    tree_label = Representative_MAG_ID_095,
    species_label = if_else(!is.na(representative_species) & representative_species != "", representative_species, "Unclassified Archaea"),
    label_face = if_else(species_label == "Unclassified Archaea", "plain", "italic"),
    cluster_size_plot = replace_na(cluster_members_095, 0),
    moz_prev_plot = replace_na(prevalence_country_Mozambique, 0),
    eth_prev_plot = replace_na(prevalence_country_Ethiopia, 0)
  ) %>%
  select(tree_label, Representative_MAG_ID_095, species_label, label_face, cluster_size_plot, moz_prev_plot, eth_prev_plot)

archaea_tree <- read.tree("06_MAGS/drep/drep_HQ_095/rep_095_gtdbtk_ar53_fasttree.nwk")
archaea_tree$tip.label <- archaea_tree$tip.label %>% basename() %>%
  str_remove("\\.fasta$") %>% str_remove("\\.fa$") %>% str_remove("\\.fna$")

archaea_label_meta <- archaea_meta_plot %>% select(tree_label, Representative_MAG_ID_095, species_label, label_face)
cluster_ring_archaea <- archaea_meta_plot %>% transmute(id = tree_label, ring = "Cluster size", cluster_size_value = cluster_size_plot)
moz_ring_archaea <- archaea_meta_plot %>% transmute(id = tree_label, ring = "Mozambique", moz_prev_value = moz_prev_plot)
eth_ring_archaea <- archaea_meta_plot %>% transmute(id = tree_label, ring = "Ethiopia", eth_prev_value = eth_prev_plot)

p_tmp <- ggtree(archaea_tree, layout = "rectangular", branch.length = "none")
tree_xmax <- max(p_tmp$data$x, na.rm = TRUE)

p_archaea_tree_base <- p_tmp %<+% archaea_label_meta +
  geom_tree(linewidth = 1) +
  geom_tiplab(aes(label = species_label, fontface = label_face), size = 7, align = TRUE, linesize = 0.1, offset = 0.05) +
  xlim(0, tree_xmax + 20) +
  theme_void() +
  theme(plot.margin = margin(5, 40, 5, 5))

p_archaea_tree_cluster <- p_archaea_tree_base +
  geom_fruit(data = cluster_ring_archaea, geom = geom_tile, mapping = aes(y = id, x = ring, fill = cluster_size_value),
             offset = 2, pwidth = 0.08, width = 1, color = "grey80") +
  scale_fill_gradientn(colors = c("snow1", "goldenrod1", "#FF4D6D"), name = "Cluster size\n(n MAGs)")

p_archaea_tree_moz <- p_archaea_tree_cluster +
  ggnewscale::new_scale_fill() +
  geom_fruit(data = moz_ring_archaea, geom = geom_tile, mapping = aes(y = id, x = ring, fill = moz_prev_value),
             offset = 0.13, pwidth = 0.08, width = 1, color = "grey80") +
  scale_fill_gradient(low = "snow1", high = "#E69F00", labels = scales::percent, name = "Mozambique\nprevalence")

p_archaea_tree_final <- p_archaea_tree_moz +
  ggnewscale::new_scale_fill() +
  geom_fruit(data = eth_ring_archaea, geom = geom_tile, mapping = aes(y = id, x = ring, fill = eth_prev_value),
             offset = 0.13, pwidth = 0.08, width = 1, color = "grey80") +
  scale_fill_gradient(low = "snow1", high = "#0072B5", labels = scales::percent, name = "Ethiopia\nprevalence") +
  theme(legend.position = "right", plot.margin = margin(5, 20, 5, 5))

ggsave("Figures/rep_095_archaea_ar53_tree_cluster_prevalence.png", p_archaea_tree_final, width = 7.5, height = 5, dpi = 300)
ggsave("Figures/rep_095_archaea_ar53_tree_cluster_prevalence.pdf", p_archaea_tree_final, width = 7.5, height = 5)
```

## Step 9 — under-classified ("putative novel") bacteria summary (Extended Data Fig. 9)

```r
n_samples_total <- country_denominator %>% summarise(n = sum(n_samples_country, na.rm = TRUE)) %>% pull(n)

rep_095_underclassified <- rep_095_tree_meta_arch %>%
  filter(novelty_category_updated != "GTDB species-assigned bacteria") %>%
  mutate(
    representative_phylum_plot = if_else(representative_domain == "Archaea", representative_phylum,
                                          if_else(is.na(representative_phylum) | representative_phylum == "", "Unclassified bacteria", representative_phylum)),
    total_positive_samples = n_positive_samples_country_Ethiopia + n_positive_samples_country_Mozambique,
    overall_prevalence = total_positive_samples / n_samples_total,
    representative_order_plot = if_else(is.na(representative_order) | representative_order == "", "Unclassified order", representative_order),
    representative_order_plot = str_wrap(representative_order_plot, width = 18)
  )

order_n_labels <- rep_095_underclassified %>%
  count(representative_phylum_plot, representative_order_plot, name = "n_new_rep_MAGs") %>%
  mutate(representative_order_plot_n = paste0(representative_order_plot, "\n(n=", n_new_rep_MAGs, ")"))

rep_095_underclassified <- rep_095_underclassified %>%
  left_join(order_n_labels, by = c("representative_phylum_plot", "representative_order_plot"))

# Phylum x novelty-category stacked bar chart
p_underclassified_phylum_country <- rep_095_underclassified %>%
  drop_na(representative_gtdbtk_classification) %>%
  mutate(representative_phylum_plot = fct_lump_n(representative_phylum_plot, n = 12, other_level = "Other phyla")) %>%
  count(representative_phylum_plot, novelty_category_updated, country_pattern) %>%
  group_by(representative_phylum_plot) %>%
  mutate(total_phylum = sum(n)) %>%
  ungroup() %>%
  mutate(representative_phylum_plot = fct_reorder(representative_phylum_plot, total_phylum)) %>%
  ggplot(aes(x = representative_phylum_plot, y = n, fill = country_pattern)) +
  geom_col(color = "black", linewidth = 0.15) +
  facet_grid(~ novelty_category_updated, scales = "free_y") +
  coord_flip() +
  scale_fill_manual(values = country_pattern_cols) +
  theme_bw() +
  labs(x = "GTDB phylum", y = "Number of under-classified representative MAGs", fill = "") +
  theme(strip.background = element_rect(fill = "grey90", color = "black"), strip.text = element_text(size = 8), legend.position = "bottom")

# Prevalence x cluster-size jittered dot plots, split Bacillota vs. other phyla for layout
p_underclassified_prevalence_size_baci <- rep_095_underclassified %>%
  filter(representative_phylum_plot == 'Bacillota') %>%
  ggplot(aes(x = representative_order_plot_n, y = overall_prevalence * 100, size = cluster_members_095, fill = country_pattern)) +
  geom_jitter(shape = 21, color = "black", alpha = 0.8, width = 0.18, height = 0) +
  facet_wrap(~ representative_phylum_plot, scales = "free_x") +
  scale_fill_manual(values = country_pattern_cols, guide = "none") +
  theme_bw() +
  labs(x = NULL, y = "MAG prevalence (%)") + ylim(0, 12.5) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8), strip.text.x = element_text(size = 8))

p_underclassified_prevalence_size_nobaci <- rep_095_underclassified %>%
  filter(representative_phylum_plot != 'Bacillota') %>%
  ggplot(aes(x = representative_order_plot_n, y = overall_prevalence * 100, size = cluster_members_095, fill = country_pattern)) +
  geom_jitter(shape = 21, color = "black", alpha = 0.8, width = 0.18, height = 0) +
  facet_grid(~ representative_phylum_plot, scales = "free_x", space = "free_x") +
  scale_size_continuous(name = "Cluster size\n(n MAGs)", range = c(1.5, 8)) +
  scale_fill_manual(values = country_pattern_cols, guide = "none") +
  ylim(0, 12.5) +
  theme_bw() +
  labs(x = NULL, y = "MAG prevalence (%)") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8), strip.text.x = element_text(size = 8))

p_underclassified_combined_up <- ggarrange(p_underclassified_prevalence_size_baci, p_underclassified_prevalence_size_nobaci, nrow = 2, heights = c(1, 1.4))
p_underclassified_combined_all <- ggarrange(p_underclassified_phylum_country, p_underclassified_combined_up, nrow = 2, heights = c(1, 2))

ggsave("Figures/rep_095_novel_bact.png", p_underclassified_combined_all, width = 7, height = 8, dpi = 300)
ggsave("Figures/rep_095_novel_bact.pdf", p_underclassified_combined_all, width = 7, height = 8)
```
