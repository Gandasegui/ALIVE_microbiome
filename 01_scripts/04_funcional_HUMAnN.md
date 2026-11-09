# 4. Read-based functional (metabolic pathway) profiling and analysis

This stage profiles host-depleted reads for metabolic pathway content
using HUMAnN, builds a community-level phyloseq object from the
unstratified pathway table, applies the same baseline/longitudinal
differential-abundance framework used for taxonomy and ARGs, and uses the
species-stratified pathway table to identify which taxa contribute to
specific pathways of interest.

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
  main_var = "Coinfection_pre",
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
