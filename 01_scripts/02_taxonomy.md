# 2. Taxonomic profiling and analysis

This stage classifies host-depleted, QC'd reads taxonomically, builds a
phyloseq object, and tests for associations between microbiome composition
and country, coinfection status, and treatment (pre/post), using the
MaAsLin3 + ANCOM-BC2 consensus framework described in Methods.

Software/database versions: Sylph v0.8.1, sylph-tax v1.2.0, GlobDB release
226, GTDB release 226 (used for benchmarking only, see Supplementary
Fig. 1C–E).

## Part A — generating taxonomic profiles (bash)

Reads are classified using Sylph against both GlobDB and GTDB. GlobDB was
selected as the primary reference database after comparing classification
rates and the number of genomes detected per sample between the two
(Supplementary Fig. 1C–E; Methods).

```bash
set -euo pipefail

UNMAPPED_DIR="02_QC/unmapped"
OUT_DIR="03_TAXONOMY"
SAMPLE_LIST="sample_list.txt"

mkdir -p "${OUT_DIR}/tax_files_GlobDB" "${OUT_DIR}/tax_files_GTDB" "${OUT_DIR}/sketches"

module load sylph/0.8.1--ha6fb395_0
module load sylph-tax/1.2.0

GLOBDB_SKETCH="/path/to/globdb_r226.syldb"     # pre-built GlobDB release 226 genome sketch
GTDB_SKETCH="/path/to/gtdb_r226.syldb"          # pre-built GTDB release 226 genome sketch
GLOBDB_TAX="/path/to/globdb_r226_taxonomy.tsv"
GTDB_TAX="/path/to/gtdb_r226_taxonomy.tsv"

while read -r sample_id; do
  R1="${UNMAPPED_DIR}/${sample_id}_qc_dedup_unmap_1.fastq.gz"
  R2="${UNMAPPED_DIR}/${sample_id}_qc_dedup_unmap_2.fastq.gz"

  # Sketch the sample's reads
  sylph sketch -1 "${R1}" -2 "${R2}" -d "${OUT_DIR}/sketches" -o "${sample_id}"

  # Profile against GlobDB (primary reference) — minimum read-sequence
  # identity of 99.5%, with unknown-sequence correction (-u), per Methods
  sylph profile "${GLOBDB_SKETCH}" "${OUT_DIR}/sketches/${sample_id}.paired.sylsp" \
    --minimum-ani 99.5 -u -t 8 \
    > "${OUT_DIR}/tax_files_GlobDB/tax_${sample_id}_qc_dedup_unmap.sylphmpa"

  sylph-tax taxprof "${OUT_DIR}/tax_files_GlobDB/tax_${sample_id}_qc_dedup_unmap.sylphmpa" \
    -t "${GLOBDB_TAX}" -o "${OUT_DIR}/tax_files_GlobDB/tax_${sample_id}_qc_dedup_unmap"

  # Profile against GTDB (benchmark comparison only)
  sylph profile "${GTDB_SKETCH}" "${OUT_DIR}/sketches/${sample_id}.paired.sylsp" \
    --minimum-ani 99.5 -u -t 8 \
    > "${OUT_DIR}/tax_files_GTDB/tax_${sample_id}_qc_dedup_unmap.sylphmpa"

  sylph-tax taxprof "${OUT_DIR}/tax_files_GTDB/tax_${sample_id}_qc_dedup_unmap.sylphmpa" \
    -t "${GTDB_TAX}" -o "${OUT_DIR}/tax_files_GTDB/tax_${sample_id}_qc_dedup_unmap"
done < "${SAMPLE_LIST}"
```

Database comparison summary (classified %, unclassified %, genomes
detected per sample), feeding Supplementary Fig. 1C–E:

```bash
OUT_CMP="${OUT_DIR}/sylph_db_comparison.tsv"
echo -e "Sample\tClassified_pct\tUnclassified_pct\tGenomes_detected\tDatabase" > "${OUT_CMP}"

for db in GlobDB GTDB; do
  for f in "${OUT_DIR}/tax_files_${db}"/*.sylphmpa; do
    sample=$(basename "$f" _qc_dedup_unmap.sylphmpa)
    classified_pct=$(awk 'NR>1{sum+=$NF} END{print sum}' "$f")
    unclassified_pct=$(echo "100 - ${classified_pct}" | bc -l)
    genomes_detected=$(($(wc -l < "$f") - 1))
    echo -e "${sample}\t${classified_pct}\t${unclassified_pct}\t${genomes_detected}\t${db}" >> "${OUT_CMP}"
  done
done
```

## Part B — building the phyloseq object (R)

Sylph's clade-abundance output is parsed into an OTU table and a taxonomy
table, which are combined with the sample metadata (one row per
`Sample_ID`, with `country`, `timepoint`, `treatment_arm`,
`Participant_ID`, and all baseline covariates described in Methods and
Table 1) into a single phyloseq object.

```r
library(tidyverse)
library(phyloseq)
source("R/functions/analysis_functions.R")

metadata_long <- read_csv("databases/metadata_long.csv", show_col_types = FALSE)

read_sylph <- function(file) {
  dat <- read.delim(file, skip = 1, stringsAsFactors = FALSE)
  sample_name <- basename(file) %>%
    str_remove("^tax_") %>%
    str_remove("_qc.*$")
  dat %>%
    select(clade_name, relative_abundance) %>%
    mutate(Sample = sample_name)
}

clean_tax <- function(x) {
  out <- gsub("^[a-z]__", "", x)
  out[out == "" | is.na(out)] <- NA
  out
}

files <- list.files("03_TAXONOMY/tax_files_GlobDB/", pattern = ".sylphmpa$", full.names = TRUE)
sylph_df <- bind_rows(lapply(files, read_sylph))

otu_df <- sylph_df %>%
  pivot_wider(names_from = clade_name, values_from = relative_abundance, values_fill = 0)

sample_ids <- otu_df$Sample
otu_matrix <- as.matrix(otu_df[, -1])
rownames(otu_matrix) <- sample_ids

taxa_names <- colnames(otu_matrix)
tax_mat <- str_split(taxa_names, "\\|", simplify = TRUE)
tax_mat <- apply(tax_mat, 2, clean_tax)
rownames(tax_mat) <- taxa_names
colnames(tax_mat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:ncol(tax_mat)]

physeq <- phyloseq(
  otu_table(otu_matrix, taxa_are_rows = FALSE),
  tax_table(tax_mat)
)
metadata <- metadata_long %>% column_to_rownames("Sample_ID")
physeq <- merge_phyloseq(physeq, sample_data(metadata))

physeq_pre  <- physeq %>% subset_samples(timepoint == "pre")
physeq_pre  <- prune_taxa(taxa_sums(physeq_pre) > 0, physeq_pre)
physeq_post <- physeq %>% subset_samples(timepoint == "post")
physeq_post <- prune_taxa(taxa_sums(physeq_post) > 0, physeq_post)

# Host-depleted read counts, needed for ANCOM-BC2 pseudo-count rescaling
reads_per_sample <- read_tsv("02_QC/read_counts.tsv", show_col_types = FALSE) %>%
  transmute(sample = sample, Reads = bowtie_hostrm_reads)
```

## Part C — baseline diversity testing

Overall and country-stratified Shannon-diversity and Bray–Curtis
PERMANOVA testing against age, sex, BMI, and STH infection-related
variables (Fig. 1C–D).

```r
num_vars_all <- c("age", "BMI", "TR_epg_pre", "HKW_epg_pre", "ASC_epg_pre")
cat_vars_all <- c("sex", "country", "Coinfection_pre", "TR_pre", "HKW_pre", "ASC_pre", "ST_pre")

alpha_overall <- run_alpha_tests(physeq_pre, num_vars_all, cat_vars_all, "Overall")
beta_overall  <- run_beta_tests(physeq_pre, num_vars_all, cat_vars_all, "Overall")

physeq_pre_eth <- subset_samples(physeq_pre, country == "Ethiopia")
physeq_pre_eth <- prune_taxa(taxa_sums(physeq_pre_eth) > 0, physeq_pre_eth)
physeq_pre_moz <- subset_samples(physeq_pre, country == "Mozambique")
physeq_pre_moz <- prune_taxa(taxa_sums(physeq_pre_moz) > 0, physeq_pre_moz)

alpha_eth <- run_alpha_tests(physeq_pre_eth, c("age", "BMI", "HKW_epg_pre"),
                              c("sex", "Coinfection_pre", "HKW_pre", "ST_pre"), "Ethiopia")
beta_eth  <- run_beta_tests(physeq_pre_eth, c("age", "BMI", "HKW_epg_pre"),
                             c("sex", "Coinfection_pre", "HKW_pre", "ST_pre"), "Ethiopia")

alpha_moz <- run_alpha_tests(physeq_pre_moz, c("age", "BMI", "TR_epg_pre", "ASC_epg_pre"),
                              c("sex", "Coinfection_pre", "TR_pre", "ASC_pre"), "Mozambique")
beta_moz  <- run_beta_tests(physeq_pre_moz, c("age", "BMI", "TR_epg_pre", "ASC_epg_pre"),
                             c("sex", "Coinfection_pre", "TR_pre", "ASC_pre"), "Mozambique")

baseline_diversity_summary <- bind_rows(alpha_overall, alpha_eth, alpha_moz, beta_overall, beta_eth, beta_moz)
write_csv(baseline_diversity_summary, "03_TAXONOMY/taxonomy_diversity_results_summary.csv")
```

## Part D — baseline differential abundance (country; within-country coinfection)

Taxa are first collapsed to phylum/family/genus/species and filtered to a
minimum relative abundance of 0.5% and minimum prevalence of 5% (Methods),
then tested with the paired MaAsLin3 + ANCOM-BC2 framework (Fig. 1E,
Extended Data Fig. 1–2).

```r
MIN_DETECT <- 0.5 / 100
MIN_PREV   <- 5 / 100

make_phy_list <- function(ps, detect = MIN_DETECT, prev = MIN_PREV) {
  list(
    phylum  = create_compositional_phyloseq(ps, "Phylum",  detect, prev),
    family  = create_compositional_phyloseq(ps, "Family",  detect, prev),
    genus   = create_compositional_phyloseq(ps, "Genus",   detect, prev),
    species = create_compositional_phyloseq(ps, "Species", detect, prev)
  )
}

# Country (all participants)
da_country <- run_da_baseline(
  phy_list = make_phy_list(physeq_pre),
  reads_df = reads_per_sample,
  main_var = "country",
  prefix = "country",
  out_root = "MaasLin3/baseline"
)

# Coinfection within Ethiopia (hookworm-S. stercoralis)
da_coinf_eth <- run_da_baseline(
  phy_list = make_phy_list(physeq_pre_eth),
  reads_df = reads_per_sample,
  main_var = "Coinfection_pre",
  prefix = "eth_coinfection",
  out_root = "MaasLin3/baseline"
)

# Coinfection within Mozambique (T. trichiura-A. lumbricoides)
da_coinf_moz <- run_da_baseline(
  phy_list = make_phy_list(physeq_pre_moz),
  reads_df = reads_per_sample,
  main_var = "Coinfection_pre",
  prefix = "moz_coinfection",
  out_root = "MaasLin3/baseline"
)
```

## Part E — longitudinal (pre/post treatment) differential abundance

Run separately for each country × treatment-arm combination, with
`Participant_ID` as a random effect to account for the paired study design
(Fig. 3A, Extended Data Fig. 6).

```r
build_arm_subset <- function(country_val, arm_val) {
  ps <- subset_samples(physeq, country == country_val & treatment_arm == arm_val)
  prune_taxa(taxa_sums(ps) > 0, ps)
}

arms <- expand_grid(country = c("Ethiopia", "Mozambique"), arm = c("ALB", "FDCx1", "FDCx3"))

da_treatment <- purrr::pmap(arms, function(country, arm) {
  ps <- build_arm_subset(country, arm)
  run_da_longitudinal(
    phy_list = make_phy_list(ps),
    reads_df = reads_per_sample,
    main_var = "timepoint",
    prefix = paste0(tolower(substr(country, 1, 3)), "_", tolower(arm)),
    out_root = "MaasLin3/treatment"
  )
})
names(da_treatment) <- paste0(arms$country, "_", arms$arm)
```

## Part F — Fig. 1C: PERMANOVA barplots (overall, Ethiopia, Mozambique)

```r
library(plyr)

beta_df <- bind_rows(
  beta_overall %>% mutate(subgroup = "Overall"),
  beta_eth %>% mutate(subgroup = "Ethiopia"),
  beta_moz %>% mutate(subgroup = "Mozambique")
)

beta_df$variable <- plyr::revalue(as.factor(beta_df$variable), c(
  "age" = "Age", "ASC_epg_pre" = "Ascaris EPG", "ASC_pre" = "Ascaris presence", "BMI" = "BMI",
  "Coinfection_pre" = "Coinfection", "country" = "Country", "HKW_epg_pre" = "Hookworm EPG",
  "HKW_pre" = "Hookworm presence", "sex" = "Sex", "ST_pre" = "Strongyloides presence",
  "TR_epg_pre" = "Trichuris EPG", "TR_pre" = "Trichuris presence"
))
beta_df$subgroup <- factor(beta_df$subgroup, levels = c("Overall", "Ethiopia", "Mozambique"))

p_beta_sig <- ggplot(beta_df, aes(x = variable, y = statistic, fill = p_adj < 0.05)) +
  geom_col(position = "dodge") +
  facet_wrap(~subgroup, scales = "free_y") +
  coord_flip() +
  theme_bw() +
  labs(x = "Variable", y = "PERMANOVA statistic", fill = "Significantly associated") +
  theme(legend.position = "bottom")

ggsave("Figures/beta_div_sig.png", p_beta_sig, width = 8, height = 3)
```

## Part G — Fig. 1D: PCoA by country, and by coinfection status within Ethiopia

```r
library(ggExtra)
library(ggpubr)

dist_bc <- phyloseq::distance(physeq_pre, method = "bray")
ord_bc <- ordinate(physeq_pre, method = "PCoA", distance = dist_bc)
eigen_vec <- rownames_to_column(data.frame(ord_bc$vectors), var = "Sample_ID")

alpha_div <- microbiome::alpha(physeq_pre, index = "diversity_shannon")
alpha_df <- data.frame(Sample_ID = rownames(alpha_div), alpha_div) %>%
  left_join(as.data.frame(sample_data(physeq_pre)) %>% rownames_to_column("Sample_ID"), by = "Sample_ID")
beta_alpha_div_df <- left_join(alpha_df, eigen_vec, by = "Sample_ID")

p_beta_country <- ggplot(beta_alpha_div_df, aes(x = Axis.1, y = Axis.2, col = country)) +
  geom_point(size = 2, alpha = 0.8, stroke = 0.3) +
  scale_color_bmj(name = "") +
  xlab("PC1") + ylab("PC2") +
  theme_bw() +
  theme(legend.position = "bottom")
p_beta_country <- ggMarginal(p_beta_country, type = "boxplot", groupColour = TRUE, groupFill = TRUE, margins = "both", size = 6)

# Within-Ethiopia PCoA, highlighting coinfection pattern
dist_bc_eth <- phyloseq::distance(physeq_pre_eth, method = "bray")
ord_bc_eth <- ordinate(physeq_pre_eth, method = "PCoA", distance = dist_bc_eth)
eigen_vec_eth <- rownames_to_column(data.frame(ord_bc_eth$vectors), var = "Sample_ID")
beta_alpha_div_df_eth <- left_join(alpha_df, eigen_vec_eth, by = "Sample_ID") %>%
  mutate(newvar = case_when(
    country == "Ethiopia" & ST_pre == "No" & Coinfection_pre == "0" ~ "Hookworm only",
    country == "Ethiopia" & ST_pre == "Yes" & Coinfection_pre == "0" ~ "Strongyloides only",
    country == "Ethiopia" & ST_pre == "Yes" & Coinfection_pre == "1" ~ "Coinfection"
  ))

p_beta_coinfect_eth <- beta_alpha_div_df_eth %>%
  filter(country == "Ethiopia") %>%
  ggplot(aes(x = Axis.1, y = Axis.2, col = newvar)) +
  geom_point(size = 3, alpha = 0.8, stroke = 0.3) +
  scale_color_manual(values = c("Hookworm only" = "blue3", "Strongyloides only" = "purple", "Coinfection" = "pink"),
                      name = NULL, labels = c("Coinf.", "Hookworm", "S. stercoralis")) +
  xlab("PC1") + ylab("PC2") +
  theme_bw() +
  theme(legend.position = "right")
p_beta_coinfect_eth <- ggMarginal(p_beta_coinfect_eth, type = "boxplot", groupColour = TRUE, groupFill = TRUE, margins = "both", size = 6)

p_beta_PCoA <- ggarrange(p_beta_country, plot.new(), p_beta_coinfect_eth, nrow = 1, widths = c(1, 0.10, 1))
ggsave("Figures/beta_div_PCoA.png", p_beta_PCoA, width = 6.5, height = 3.5)
```

## Part H — Fig. 1E: MaAsLin3 coefficient plot for the consensus (MaAsLin3 + ANCOM-BC2) baseline taxa

```r
# ancom_sig: significant ANCOM-BC2 features across ranks (from da_country$ancom)
# significant_results: MaAsLin3 significant_results.tsv rows across ranks (from read_maaslin_results)
ancom_sig <- extract_ancom_results(da_country$ancom) %>%
  filter(diff_robust_countryMozambique == TRUE) %>%
  transmute(taxon = str_remove(feature, "_A$"), taxa_level = level)

significant_results <- read_maaslin_results("country", c("phylum", "family", "genus", "species"), "MaasLin3/baseline") %>%
  filter(!is.na(coef)) %>%
  mutate(significance = ifelse(qval_individual <= 0.01, "Significant", "Not Significant"),
         feature = str_replace(feature, "_A|_R", ""))

shared_sig <- ancom_sig %>%
  inner_join(significant_results, by = c("taxon" = "feature"))

p_DA_baseline <- ggplot(shared_sig, aes(x = coef, y = taxon)) +
  geom_hline(aes(yintercept = taxon), color = "lightgray", linetype = "dotted") +
  geom_point(aes(color = significance, shape = model), size = 3) +
  geom_errorbarh(aes(xmin = coef - stderr, xmax = coef + stderr), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "MaAsLin3 coef. (> 0, higher in Moz)", y = "Microbial Feature") +
  facet_wrap(~taxa_level, scales = "free_y", ncol = 2) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10, face = "italic"), legend.position = "right",
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())

ggsave("Figures/DA_baseline.png", p_DA_baseline, width = 8, height = 5)
```

This yields **24 consensus (MaAsLin3 + ANCOM-BC2-confirmed) taxa** across
phylum/family/genus/species, matching the manuscript Results and Fig. 1E.

## Part I — relaxed cross-country / cross-arm intersection analysis

Because requiring concordant significance by both MaAsLin3 and ANCOM-BC2
is conservative, features significant by MaAsLin3 alone are additionally
intersected across countries (within each treatment arm) and across
FDCx1/FDCx3 (within each country), feeding the UpSet plots in Fig. 5A.

```r
sig_sets <- purrr::map(names(da_treatment), function(nm) {
  parts <- str_split(nm, "_", simplify = TRUE)
  prefix <- paste0(tolower(substr(parts[1], 1, 3)), "_", tolower(parts[2]))
  res <- read_maaslin_results(prefix, c("phylum", "family", "genus", "species"), "MaasLin3/treatment")
  unique(res$feature)
})
names(sig_sets) <- names(da_treatment)

# `sig_sets` is a named list of character vectors (significant taxa per
# country x arm combination) — feed directly into ComplexUpset for Fig. 5A,
# or into intersect()/setdiff() for pairwise comparisons.
```

Fig. 5A itself (UpSet plot across the four country x arm combinations):

```r
library(ComplexUpset)

upset_df <- tibble(feature = unique(unlist(sig_sets))) %>%
  mutate(
    `ETH-FDCx1` = feature %in% sig_sets[["Ethiopia_FDCx1"]],
    `ETH-FDCx3` = feature %in% sig_sets[["Ethiopia_FDCx3"]],
    `MOZ-FDCx1` = feature %in% sig_sets[["Mozambique_FDCx1"]],
    `MOZ-FDCx3` = feature %in% sig_sets[["Mozambique_FDCx3"]]
  )

p_upset_taxonomy <- ComplexUpset::upset(
  upset_df,
  intersect = c("ETH-FDCx1", "ETH-FDCx3", "MOZ-FDCx1", "MOZ-FDCx3"),
  name = "",
  base_annotations = list("Intersection size" = intersection_size(counts = TRUE))
)

ggsave("Figures/upset_taxonomy_crossarm.png", p_upset_taxonomy, width = 8, height = 5)
```

