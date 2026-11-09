# 3. Resistome (ARG) profiling and analysis

This stage identifies antimicrobial resistance genes from host-depleted
reads, builds RPKM-normalised and count-scale phyloseq objects at the
individual-gene level, and applies the same baseline/longitudinal
differential-abundance framework used for taxonomy.

Software/database versions: RGI bwt module v6.0.8; CARD data v4.0.1; CARD
prevalence/resistance/variant data v4.0.2.

## Part A — generating ARG profiles (bash)

Reads are aligned directly against CARD reference sequences using RGI in
"bwt" (reads) mode.

```bash
set -euo pipefail

UNMAPPED_DIR="02_QC/unmapped"
OUT_DIR="04_ARG/rgi_bwt"
SAMPLE_LIST="sample_list.txt"
CARD_DB="/path/to/card_v4.0.1"   # localized CARD reference database

mkdir -p "${OUT_DIR}"

module load rgi/6.0.8
rgi load --card_json "${CARD_DB}/card.json" --local
rgi card_annotation -i "${CARD_DB}/card.json" > card_annotation.log 2>&1

while read -r sample_id; do
  R1="${UNMAPPED_DIR}/${sample_id}_qc_dedup_unmap_1.fastq.gz"
  R2="${UNMAPPED_DIR}/${sample_id}_qc_dedup_unmap_2.fastq.gz"

  rgi bwt \
    --read_one "${R1}" \
    --read_two "${R2}" \
    --aligner bowtie2 \
    --output_file "${OUT_DIR}/${sample_id}" \
    --local \
    --threads 8
done < "${SAMPLE_LIST}"
```

`rgi bwt` writes `<sample_id>.gene_mapping_data.txt` per sample, containing,
per ARO term: "Completely Mapped Reads", "Mapped Reads with Flanking
Sequence" (summed as "All Mapped Reads"), AMR Gene Family, Drug Class,
Resistance Mechanism, and Reference Length.

## Part B — parsing RGI output and normalising to RPKM (R)

Gene mappings supported by fewer than 25 reads are dropped (Methods), then
normalised to reads per kilobase of gene length per million host-depleted
reads.

```r
library(tidyverse)
library(phyloseq)
source("R/functions/analysis_functions.R")

metadata_long <- read_csv("databases/metadata_long.csv", show_col_types = FALSE)
reads_per_sample <- read_tsv("02_QC/read_counts.tsv", show_col_types = FALSE) %>%
  transmute(sample = sample, Reads = bowtie_hostrm_reads)

rgi_anno_from_file <- function(path) {
  x <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
  sample_id <- basename(path) %>% str_remove("\\.gene_mapping_data\\.txt$")

  if ("All Mapped Reads" %in% names(x)) {
    count_vec <- x[["All Mapped Reads"]]
  } else {
    count_vec <- as.numeric(x[["Completely Mapped Reads"]]) +
      as.numeric(x[["Mapped Reads with Flanking Sequence"]])
  }

  tibble(
    sample          = sample_id,
    gene_name       = x[["ARO Term"]],
    ARO             = paste0("ARO:", x[["ARO Accession"]]),
    AMR_gene_family = x[["AMR Gene Family"]],
    drug_class      = x[["Drug Class"]],
    resistance_mech = x[["Resistance Mechanism"]],
    ref_len_bp      = as.numeric(x[["Reference Length"]]),
    count           = as.numeric(count_vec)
  ) %>%
    filter(!is.na(gene_name), gene_name != "", !is.na(ARO), ARO != "ARO:", !is.na(count), count > 0)
}

rgi_files <- list.files("04_ARG/rgi_bwt", pattern = "\\.gene_mapping_data\\.txt$", full.names = TRUE)
rgi_long <- map_dfr(rgi_files, rgi_anno_from_file) %>%
  group_by(sample, gene_name, ARO, AMR_gene_family, drug_class, resistance_mech, ref_len_bp) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  filter(count >= 25) %>%
  mutate(gene_id = paste0("RGI|", ARO, "|", gene_name))

rgi_norm <- rgi_long %>%
  left_join(reads_per_sample, by = "sample") %>%
  mutate(
    rpkm = if_else(
      !is.na(ref_len_bp) & ref_len_bp > 0 & !is.na(Reads) & Reads > 0,
      (count * 1e9) / (ref_len_bp * Reads),
      0
    )
  )
```

## Part C — building gene-level phyloseq objects (counts and RPKM)

```r
build_matrix <- function(df, value_col) {
  df %>%
    select(gene_id, sample, !!sym(value_col)) %>%
    pivot_wider(names_from = sample, values_from = !!sym(value_col), values_fill = 0) %>%
    column_to_rownames("gene_id") %>%
    as.matrix()
}

gene_counts_mat <- build_matrix(rgi_norm, "count")
gene_rpkm_mat   <- build_matrix(rgi_norm, "rpkm")

gene_tax <- rgi_norm %>%
  distinct(gene_id, gene_name, ARO, AMR_gene_family, drug_class, resistance_mech, ref_len_bp) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()
gene_tax <- cbind(Gene = rownames(gene_tax), gene_tax)

common_samples <- intersect(colnames(gene_counts_mat), metadata_long$Sample_ID)
meta_ps <- metadata_long %>%
  filter(Sample_ID %in% common_samples) %>%
  distinct(Sample_ID, .keep_all = TRUE) %>%
  column_to_rownames("Sample_ID")

physeq_ARG_counts <- phyloseq(
  otu_table(gene_counts_mat[, common_samples], taxa_are_rows = TRUE),
  tax_table(gene_tax),
  sample_data(meta_ps)
)
physeq_ARG_rpkm <- phyloseq(
  otu_table(gene_rpkm_mat[, common_samples], taxa_are_rows = TRUE),
  tax_table(gene_tax),
  sample_data(meta_ps)
)

physeq_ARG_counts_pre <- subset_samples(physeq_ARG_counts, timepoint == "pre")
physeq_ARG_counts_pre <- prune_taxa(taxa_sums(physeq_ARG_counts_pre) > 0, physeq_ARG_counts_pre)
physeq_ARG_rpkm_pre <- subset_samples(physeq_ARG_rpkm, timepoint == "pre")
physeq_ARG_rpkm_pre <- prune_taxa(taxa_sums(physeq_ARG_rpkm_pre) > 0, physeq_ARG_rpkm_pre)
```

## Part D — baseline diversity and differential abundance

ARGs are required to occur in at least 5% of samples (Methods; no separate
relative-abundance threshold is applied, since RPKM is already a
normalised scale). Feeds Supplementary Fig. 4 and Extended Data Fig. 3.

```r
MIN_PREV_ARG <- 0.05

num_vars_all <- c("age", "BMI", "TR_epg_pre", "HKW_epg_pre", "ASC_epg_pre")
cat_vars_all <- c("sex", "country", "Coinfection_pre", "TR_pre", "HKW_pre", "ASC_pre", "ST_pre")

alpha_ARG_overall <- run_alpha_tests(physeq_ARG_rpkm_pre, num_vars_all, cat_vars_all, "Overall")
beta_ARG_overall  <- run_beta_tests(physeq_ARG_rpkm_pre, num_vars_all, cat_vars_all, "Overall")

physeq_ARG_rpkm_pre_eth <- subset_samples(physeq_ARG_rpkm_pre, country == "Ethiopia")
physeq_ARG_rpkm_pre_moz <- subset_samples(physeq_ARG_rpkm_pre, country == "Mozambique")

physeq_ARG_counts_pre_f <- filter_by_prevalence(physeq_ARG_counts_pre, MIN_PREV_ARG)
physeq_ARG_rpkm_pre_f   <- filter_by_prevalence(physeq_ARG_rpkm_pre, MIN_PREV_ARG)
common_genes <- intersect(taxa_names(physeq_ARG_counts_pre_f), taxa_names(physeq_ARG_rpkm_pre_f))
physeq_ARG_counts_pre_f <- prune_taxa(common_genes, physeq_ARG_counts_pre_f)
physeq_ARG_rpkm_pre_f   <- prune_taxa(common_genes, physeq_ARG_rpkm_pre_f)

da_ARG_country <- run_da_baseline(
  phy_list = list(gene = physeq_ARG_rpkm_pre_f),
  reads_df = reads_per_sample,
  main_var = "country",
  tax_level_field = "Gene",
  prefix = "country",
  out_root = "MaasLin3/ARG/baseline"
)
physeq_ARG_eth_rpkm_pre_f <- subset_samples(physeq_ARG_rpkm_pre_f, country == "Ethiopia")
da_ARG_coinf_eth <- run_da_baseline(
  phy_list = list(pathway = physeq_ARG_eth_rpkm_pre_f),
  reads_df = reads_per_sample,
  main_var = "Coinfection_pre",
  tax_level_field = "Feature",
  prefix = "eth_coinfection",
  out_root = "MaasLin3/HUMAnN/baseline"
)
physeq_ARG_moz_rpkm_pre_f <- subset_samples(physeq_ARG_rpkm_pre_f, country == "Mozambique")
da_path_coinf_moz <- run_da_baseline(
  phy_list = list(pathway = physeq_ARG_moz_rpkm_pre_f),
  reads_df = reads_per_sample,
  main_var = "ASC_epg_pre",
  tax_level_field = "Feature",
  prefix = "moz_coinfection",
  out_root = "MaasLin3/HUMAnN/baseline"
)
```

## Part E — longitudinal (pre/post treatment) ARG differential abundance

Run per country × treatment-arm combination (Fig. 3B, Extended Data
Fig. 6). Feeds the 15 ARGs identified as significantly decreased after
FDCx3 treatment (Fig. 3B) and the tetracycline-resistance genes found to
increase concordantly across FDC regimens (Fig. 5C).

```r
build_ARG_arm_subset <- function(country_val, arm_val) {
  ps_counts <- subset_samples(physeq_ARG_counts, country == country_val & treatment_arm == arm_val)
  ps_counts <- prune_taxa(taxa_sums(ps_counts) > 0, ps_counts)
  ps_rpkm <- subset_samples(physeq_ARG_rpkm, country == country_val & treatment_arm == arm_val)
  ps_rpkm <- prune_taxa(taxa_sums(ps_rpkm) > 0, ps_rpkm)
  ps_rpkm_f <- filter_by_prevalence(ps_rpkm, MIN_PREV_ARG)
  list(counts = prune_taxa(taxa_names(ps_rpkm_f), ps_counts), rpkm = ps_rpkm_f)
}

arms <- expand_grid(country = c("Ethiopia", "Mozambique"), arm = c("ALB", "FDCx1", "FDCx3"))

da_ARG_treatment <- purrr::pmap(arms, function(country, arm) {
  subs <- build_ARG_arm_subset(country, arm)
  run_da_longitudinal(
    phy_list = list(gene = subs$rpkm),
    reads_df = reads_per_sample,
    main_var = "timepoint",
    tax_level_field = "Gene",
    prefix = paste0(tolower(substr(country, 1, 3)), "_", tolower(arm)),
    out_root = "MaasLin3/ARG/treatment"
  )
})
names(da_ARG_treatment) <- paste0(arms$country, "_", arms$arm)
```

## Part F — Fig. 3B: coefficient plot of the 15 consensus ARGs decreased after FDCx3 (Ethiopia), with drug-class/resistance-mechanism annotation

```r
library(ggnewscale)
library(cowplot)

eth_fdcx3 <- da_ARG_treatment[["Ethiopia_FDCx3"]]

ancom_sig_arg <- extract_ancom_results(eth_fdcx3$ancom) %>%
  filter(diff_robust_timepointpost == TRUE) %>%
  transmute(feature = taxon)

maaslin_sig_arg <- read_maaslin_results("eth_fdcx3", "gene", "MaasLin3/ARG/treatment") %>%
  filter(!is.na(coef)) %>%
  mutate(significance = ifelse(qval_individual <= 0.01, "Significant", "Not Significant"))

df_plot_arg <- ancom_sig_arg %>%
  inner_join(maaslin_sig_arg, by = "feature") %>%
  left_join(
    rgi_norm %>% distinct(gene_id, gene_name, AMR_gene_family, drug_class, resistance_mech),
    by = c("feature" = "gene_id")
  )

drug_abbrev_map <- c(
  "aminoglycoside antibiotic" = "AMG", "aminocoumarin antibiotic" = "AMU", "carbapenem" = "IPM",
  "cephalosporin" = "CEF", "diaminopyrimidine antibiotic" = "TMP", "fluoroquinolone antibiotic" = "FLO",
  "glycylcycline" = "TET", "lincosamide antibiotic" = "CLI", "macrolide antibiotic" = "MAC",
  "nitrofuran antibiotic" = "NIT", "peptide antibiotic" = "PEP", "rifamycin antibiotic" = "RIF",
  "sulfonamide antibiotic" = "SLF", "tetracycline antibiotic" = "TET", "penicillin beta-lactam" = "PBL"
)

abbrev_drug_class <- function(x) {
  x %>% str_split(";") %>% lapply(str_trim) %>%
    lapply(function(classes) { m <- drug_abbrev_map[classes]; m[is.na(m)] <- "OTHERS"; m }) %>%
    sapply(function(v) paste(unique(v), collapse = " / "))
}

df_plot_arg <- df_plot_arg %>%
  mutate(drug_class_abbrev = abbrev_drug_class(drug_class)) %>%
  arrange(desc(coef)) %>%
  mutate(y = row_number())

drug_combo_pal <- setNames(grDevices::hcl.colors(n_distinct(df_plot_arg$drug_class_abbrev), palette = "Roma"),
                            unique(df_plot_arg$drug_class_abbrev))

p_annot_arg <- ggplot(df_plot_arg, aes(y = y)) +
  geom_text(aes(x = 1.095, label = gene_name), hjust = 1, size = 3) +
  geom_tile(aes(x = 0.9, fill = drug_class_abbrev), color = "black", width = 0.05, height = 0.8) +
  scale_fill_manual(name = "Drug Class", values = drug_combo_pal, na.value = "grey90") +
  new_scale_fill() +
  geom_tile(aes(x = 0.95, fill = resistance_mech), color = "black", width = 0.05, height = 0.8) +
  scale_fill_discrete(name = "Resistance mechanism") +
  scale_y_continuous(breaks = df_plot_arg$y, labels = NULL) +
  xlim(0.85, 1.1) +
  theme_void() +
  theme(legend.position = "left")

p_coef_arg <- ggplot(df_plot_arg, aes(x = coef, y = y)) +
  geom_point(aes(color = significance, shape = model), size = 3) +
  geom_errorbarh(aes(xmin = coef - stderr, xmax = coef + stderr), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_y_continuous(breaks = df_plot_arg$y, labels = NULL) +
  labs(x = "MaAsLin3 coef. (>0, decrease after treatment)", y = NULL) +
  theme_bw() +
  theme(axis.text.y = element_blank(), legend.position = "right")

p_ARG_DA <- cowplot::plot_grid(p_annot_arg, p_coef_arg, ncol = 2, rel_widths = c(1, 1.5))

ggsave("Figures/ARG_DA_treatment_FDCx3.png", p_ARG_DA, width = 10, height = 6)
```

This yields the **15 consensus ARGs** (confirmed by both MaAsLin3 and
ANCOM-BC2) found to decrease after FDCx3 in Ethiopia, matching the
manuscript Results and Fig. 3B — these are the same 15 genes passed to
module 8's Fig. 4D analysis to confirm which Enterobacteriaceae species
carry them at the MAG level.

