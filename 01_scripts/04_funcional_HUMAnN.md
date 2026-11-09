# 4. Functional (metabolic pathway) profiling and analysis

This stage profiles host-depleted reads for metabolic pathway content
using HUMAnN, builds a community-level phyloseq object from the
unstratified pathway table, applies the same baseline/longitudinal
differential-abundance framework used for taxonomy and ARGs, and uses the
species-stratified pathway table to identify which taxa contribute to
specific pathways of interest.

Software versions: HUMAnN v3.8; MetaPhlAn v4.0.6 with the
mpa_vOct22_CHOCOPhlAnSGB_202212 marker database; full ChocoPhlAn
nucleotide and UniRef90 protein databases.

## Part A — generating functional profiles (bash)

HUMAnN expects a single concatenated input per sample for paired-end data.

```bash
set -euo pipefail

UNMAPPED_DIR="02_QC/unmapped"
OUT_DIR="05_HUMAnN"
SAMPLE_LIST="sample_list.txt"

mkdir -p "${OUT_DIR}/concat" "${OUT_DIR}/per_sample"

module load humann/3.8
module load metaphlan/4.0.6

export METAPHLAN_DB_DIR="/path/to/metaphlan_db"
export CHOCOPHLAN_DB_DIR="/path/to/chocophlan"
export UNIREF90_DB_DIR="/path/to/uniref90"

while read -r sample_id; do
  R1="${UNMAPPED_DIR}/${sample_id}_qc_dedup_unmap_1.fastq.gz"
  R2="${UNMAPPED_DIR}/${sample_id}_qc_dedup_unmap_2.fastq.gz"
  CONCAT="${OUT_DIR}/concat/${sample_id}_concat.fastq.gz"

  cat "${R1}" "${R2}" > "${CONCAT}"

  humann \
    --input "${CONCAT}" \
    --output "${OUT_DIR}/per_sample/${sample_id}" \
    --metaphlan-options "--index mpa_vOct22_CHOCOPhlAnSGB_202212 --bowtie2db ${METAPHLAN_DB_DIR}" \
    --nucleotide-database "${CHOCOPHLAN_DB_DIR}" \
    --protein-database "${UNIREF90_DB_DIR}" \
    --threads 8

  rm -f "${CONCAT}"
done < "${SAMPLE_LIST}"
```

Merge per-sample pathway abundance tables and normalise to CPM (the gold
standard normalisation for HUMAnN):

```bash
humann_join_tables \
  --input "${OUT_DIR}/per_sample" \
  --output "${OUT_DIR}/pathabundance_unstratified_raw.tsv" \
  --file_name pathabundance \
  --search-subdirectories

humann_renorm_table \
  --input "${OUT_DIR}/pathabundance_unstratified_raw.tsv" \
  --output "${OUT_DIR}/pathabundance_cpm_stratified.tsv" \
  --units cpm

humann_split_stratified_table \
  --input "${OUT_DIR}/pathabundance_cpm_stratified.tsv" \
  --output "${OUT_DIR}"

mv "${OUT_DIR}/pathabundance_cpm_stratified_unstratified.tsv" "${OUT_DIR}/pathabundance_cpm_unstratified.tsv"
```

## Part B — building the phyloseq object (R, unstratified / community-level)

The UNMAPPED and UNINTEGRATED summary rows are excluded, since they are
HUMAnN QC values, not biological pathways.

```r
library(tidyverse)
library(phyloseq)
source("R/functions/analysis_functions.R")

metadata_long <- read_csv("databases/metadata_long.csv", show_col_types = FALSE)
reads_per_sample <- read_tsv("02_QC/read_counts.tsv", show_col_types = FALSE) %>%
  transmute(sample = sample, Reads = bowtie_hostrm_reads)

path_abund <- read_tsv("05_HUMAnN/pathabundance_cpm_unstratified.tsv", show_col_types = FALSE)
colnames(path_abund)[1] <- "pathway"
path_abund <- path_abund %>%
  filter(!pathway %in% c("UNMAPPED", "UNINTEGRATED"))
colnames(path_abund) <- colnames(path_abund) %>%
  str_replace("_Abundance$", "") %>%
  str_replace("_CPM$", "")

path_mat <- path_abund %>% column_to_rownames("pathway") %>% as.matrix()
common_samples <- intersect(colnames(path_mat), metadata_long$Sample_ID)
path_mat <- path_mat[, common_samples, drop = FALSE]

meta_sub <- metadata_long %>%
  filter(Sample_ID %in% common_samples) %>%
  distinct(Sample_ID, .keep_all = TRUE) %>%
  arrange(match(Sample_ID, colnames(path_mat))) %>%
  column_to_rownames("Sample_ID")

path_tax <- tibble(pathway = rownames(path_mat), Feature = rownames(path_mat)) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

physeq_path <- phyloseq(
  otu_table(path_mat, taxa_are_rows = TRUE),
  tax_table(path_tax),
  sample_data(meta_sub)
)

MIN_PREV_PATH <- 0.05
physeq_path_f <- filter_by_prevalence(physeq_path, MIN_PREV_PATH)

physeq_path_pre <- subset_samples(physeq_path_f, timepoint == "pre")
physeq_path_pre <- prune_taxa(taxa_sums(physeq_path_pre) > 0, physeq_path_pre)
```

## Part C — baseline diversity and differential abundance

```r
num_vars_all <- c("age", "BMI", "TR_epg_pre", "HKW_epg_pre", "ASC_epg_pre")
cat_vars_all <- c("sex", "country", "Coinfection_pre", "TR_pre", "HKW_pre", "ASC_pre", "ST_pre")

alpha_path_overall <- run_alpha_tests(physeq_path_pre, num_vars_all, cat_vars_all, "Overall")
beta_path_overall  <- run_beta_tests(physeq_path_pre, num_vars_all, cat_vars_all, "Overall")

physeq_path_pre_eth <- subset_samples(physeq_path_pre, country == "Ethiopia")
physeq_path_pre_moz <- subset_samples(physeq_path_pre, country == "Mozambique")

da_path_country <- run_da_baseline(
  phy_list = list(pathway = physeq_path_pre),
  reads_df = reads_per_sample,
  main_var = "country",
  tax_level_field = "Feature",
  prefix = "baseline_country",
  out_root = "MaasLin3/HUMAnN/baseline"
)

da_path_coinf_eth <- run_da_baseline(
  phy_list = list(pathway = physeq_path_pre_eth),
  reads_df = reads_per_sample,
  main_var = "Coinfection_pre",
  tax_level_field = "Feature",
  prefix = "eth_coinfection",
  out_root = "MaasLin3/HUMAnN/baseline"
)

da_path_coinf_moz <- run_da_baseline(
  phy_list = list(pathway = physeq_path_pre_moz),
  reads_df = reads_per_sample,
  main_var = "ASC_epg_pre",
  tax_level_field = "Feature",
  prefix = "moz_coinfection",
  out_root = "MaasLin3/HUMAnN/baseline"
)
```

## Part D — longitudinal (pre/post treatment) pathway differential abundance

```r
build_path_arm_subset <- function(country_val, arm_val) {
  ps <- subset_samples(physeq_path_f, country == country_val & treatment_arm == arm_val)
  prune_taxa(taxa_sums(ps) > 0, ps)
}

arms <- expand_grid(country = c("Ethiopia", "Mozambique"), arm = c("ALB", "FDCx1", "FDCx3"))

da_path_treatment <- purrr::pmap(arms, function(country, arm) {
  ps <- build_path_arm_subset(country, arm)
  run_da_longitudinal(
    phy_list = list(pathway = ps),
    reads_df = reads_per_sample,
    main_var = "timepoint",
    tax_level_field = "Feature",
    prefix = paste0(tolower(substr(country, 1, 3)), "_", tolower(arm)),
    out_root = "MaasLin3/HUMAnN/treatment"
  )
})
names(da_path_treatment) <- paste0(arms$country, "_", arms$arm)
```

## Part E — pathway-category classification (Fig. 3C)

Significant pathways are grouped into broad functional categories (amino
acid metabolism, redox metabolism, Gram-negative cell wall, etc.) based on
their MetaCyc annotation text.

```r
classify_pathway_theme <- function(pathway) {
  p <- str_to_lower(pathway)
  case_when(
    str_detect(p, "peptidoglycan|muramoyl|cell wall|lipopolysaccharide|o-antigen|colanic acid|kdo|polymyxin resistance") ~ "Gram-negative cell wall",
    str_detect(p, "fatty acid|palmitate|oleate|phospholipid|beta-oxidation|fatty acids biosynthesis|fatty acid elongation") ~ "Lipid metabolism",
    str_detect(p, "glycolysis|tca cycle|glyoxylate|fermentation|nitrate reduction|anaerobic energy metabolism") ~ "Energy metabolism",
    str_detect(p, "arginine|aspartate|lysine|histidine|methionine|isoleucine|valine|tryptophan|threonine|branched chain amino acid|glutamate|glutamine") ~ "Amino acid metabolism",
    str_detect(p, "biotin|coenzyme a|phosphopantothenate|heme|menaquin|nad salvage|nad de novo|thiamine|folate|queuosine|s-adenosyl-l-methionine") ~ "Redox metabolism",
    str_detect(p, "purine|pyrimidine|nucleotide|nucleoside|ump|guanosine|adenine|t.?rna charging|ribonucleotide") ~ "Translation-related metabolism",
    str_detect(p, "catechol|protocatechuate|salicylate|aromatic|phosphonate") ~ "Aromatic degradation",
    str_detect(p, "galactarate|galacturonate|glucarate|glucuron|fructuronate|xylose|mannan|chitin|ascorbate") ~ "Host-diet substrate metabolism",
    TRUE ~ "Other"
  )
}
```

## Part F — species-stratified pathway contributors (Fig. 3D, Extended Data Fig. 7)

For a pathway of interest, this returns the mean species-stratified
contribution (CPM), before vs. after treatment, for the top contributing
taxa — used to confirm that a given pathway's decrease is predominantly
attributable to *E. coli*.

```r
path_strat <- read_tsv("05_HUMAnN/pathabundance_cpm_stratified.tsv", show_col_types = FALSE)
colnames(path_strat)[1] <- "pathway"
colnames(path_strat) <- colnames(path_strat) %>%
  str_replace("_Abundance$", "") %>%
  str_replace("_CPM$", "")

get_pathway_contributors <- function(pathway_id, metadata, top_n = 5) {
  path_strat %>%
    filter(str_detect(pathway, fixed(pathway_id))) %>%
    filter(str_detect(pathway, "\\|")) %>%
    mutate(taxon = str_replace(pathway, "^.*\\|", "")) %>%
    pivot_longer(cols = -c(pathway, taxon), names_to = "Sample_ID", values_to = "CPM") %>%
    left_join(metadata, by = "Sample_ID") %>%
    group_by(taxon, timepoint) %>%
    summarise(mean_CPM = mean(CPM, na.rm = TRUE), .groups = "drop") %>%
    group_by(taxon) %>%
    mutate(max_CPM = max(mean_CPM, na.rm = TRUE)) %>%
    ungroup() %>%
    slice_max(max_CPM, n = top_n * 2) # top_n taxa x 2 timepoints
}
```

## Part G — Fig. 3C: number of significant pathways per functional theme

```r
eth_fdcx3_path <- da_path_treatment[["Ethiopia_FDCx3"]]

ancom_sig_path <- extract_ancom_results(eth_fdcx3_path$ancom) %>%
  filter(diff_robust_timepointpost == TRUE) %>%
  transmute(feature = taxon)

maaslin_sig_path <- read_maaslin_results("eth_fdcx3", "pathway", "MaasLin3/HUMAnN/treatment") %>%
  filter(!is.na(coef)) %>%
  mutate(significance = ifelse(qval_individual <= 0.01, "Significant", "Not Significant"))

consensus_sig_path <- ancom_sig_path %>%
  inner_join(maaslin_sig_path, by = "feature") %>%
  mutate(theme = classify_pathway_theme(feature))

p_theme_counts <- consensus_sig_path %>%
  count(theme) %>%
  mutate(theme = fct_reorder(theme, n)) %>%
  ggplot(aes(x = theme, y = n)) +
  geom_col(fill = "#C4B7A6") +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = "Number of significant pathways")

ggsave("Figures/pathway_theme_counts_FDCx3.png", p_theme_counts, width = 6, height = 4)
```

## Part H — Fig. 3D: species-stratified contribution to selected treatment-decreased pathways

```r
pathways_of_interest <- c("D-galactarate Deg I", "FA elongation", "lipid IVA BS", "SP of D-glucarate and D-galactarate Deg")

contributor_plots <- purrr::map(pathways_of_interest, function(p) {
  df <- get_pathway_contributors(p, meta_sub, top_n = 5) %>%
    mutate(timepoint = factor(timepoint, levels = c("pre", "post"), labels = c("Pre-treatment", "Post-treatment")))

  ggplot(df, aes(x = mean_CPM, y = fct_reorder(taxon, max_CPM), fill = timepoint)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pre-treatment" = "orange", "Post-treatment" = "darkred")) +
    labs(title = p, x = "CPM", y = NULL) +
    theme_bw() +
    theme(legend.position = "bottom")
})

p_pathway_contributors <- patchwork::wrap_plots(contributor_plots, ncol = 2) + patchwork::plot_layout(guides = "collect")
ggsave("Figures/pathway_species_contributors_FDCx3.png", p_pathway_contributors, width = 10, height = 8)
```

## Part I — Fig. 5B: relaxed (MaAsLin3-only) cross-country / cross-arm pathway intersection

```r
sig_sets_path <- purrr::map(names(da_path_treatment), function(nm) {
  parts <- str_split(nm, "_", simplify = TRUE)
  prefix <- paste0(tolower(substr(parts[1], 1, 3)), "_", tolower(parts[2]))
  res <- read_maaslin_results(prefix, "pathway", "MaasLin3/HUMAnN/treatment")
  unique(res$feature)
})
names(sig_sets_path) <- names(da_path_treatment)

library(ComplexUpset)

upset_df_path <- tibble(feature = unique(unlist(sig_sets_path))) %>%
  mutate(
    `ETH-FDCx1` = feature %in% sig_sets_path[["Ethiopia_FDCx1"]],
    `ETH-FDCx3` = feature %in% sig_sets_path[["Ethiopia_FDCx3"]],
    `MOZ-FDCx1` = feature %in% sig_sets_path[["Mozambique_FDCx1"]],
    `MOZ-FDCx3` = feature %in% sig_sets_path[["Mozambique_FDCx3"]]
  )

p_upset_pathways <- ComplexUpset::upset(
  upset_df_path,
  intersect = c("ETH-FDCx1", "ETH-FDCx3", "MOZ-FDCx1", "MOZ-FDCx3"),
  name = "",
  base_annotations = list("Intersection size" = intersection_size(counts = TRUE))
)

ggsave("Figures/upset_pathways_crossarm.png", p_upset_pathways, width = 8, height = 5)
```

## Part J — Fig. 5E: the 28-pathway core significant in all four country x arm combinations

```r
core_pathways <- Reduce(intersect, sig_sets_path)
length(core_pathways) # expect 28, per manuscript Results

core_pathway_coefs <- purrr::map_dfr(names(da_path_treatment), function(nm) {
  parts <- str_split(nm, "_", simplify = TRUE)
  prefix <- paste0(tolower(substr(parts[1], 1, 3)), "_", tolower(parts[2]))
  read_maaslin_results(prefix, "pathway", "MaasLin3/HUMAnN/treatment") %>%
    filter(feature %in% core_pathways) %>%
    transmute(feature, coef, comparison = nm)
})

p_core_pathway_heatmap <- core_pathway_coefs %>%
  mutate(
    theme = classify_pathway_theme(feature),
    feature_short = str_wrap(str_replace(feature, "^.*?:\\s*", ""), width = 40)
  ) %>%
  ggplot(aes(x = comparison, y = feature_short, fill = coef)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", coef)), size = 2.8) +
  scale_fill_gradient(low = "white", high = "#355F5B", name = "MaAsLin3\ncoefficient") +
  facet_grid(theme ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text.y = element_text(angle = 0))

ggsave("Figures/core_28_pathways_heatmap.png", p_core_pathway_heatmap, width = 9, height = 10)
```

