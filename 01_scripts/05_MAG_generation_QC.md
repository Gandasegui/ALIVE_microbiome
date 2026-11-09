# 5. Metagenome-assembled genome (MAG) generation and quality control

This stage assembles, bins, refines, taxonomically classifies, and
quality-filters MAGs from every sample's host-depleted reads, and produces
the Supplementary QC figure summarising assembly/binning/decontamination
performance across the cohort.

Two purpose-built Nextflow pipelines are used in sequence: **mags_maker**
(assembly, binning, bin refinement) and **qc_mags** (taxonomic
classification, quality assessment, decontamination). Software versions:
mags_maker v1.1.1 (metaWRAP v1.3.2-based); qc_mags v1.0.0; GTDB-Tk (GTDB
release 226); QUAST v5.3.0; CheckM2 v1.0.2; GUNC v1.0.6; MDMCleaner v0.8.7;
seqkit v2.10.0.

## Step 1 — mags_maker: assembly, binning, and bin refinement

```bash
set -euo pipefail

UNMAPPED_DIR="02_QC/unmapped"
OUT_DIR="06_MAGS"
SAMPLE_LIST="sample_list.txt"

mkdir -p "${OUT_DIR}/mags_maker"

MANIFEST="${OUT_DIR}/manifest_mags_maker.csv"
echo "ID,R1,R2" > "${MANIFEST}"
while read -r sample_id; do
  echo "${sample_id},${UNMAPPED_DIR}/${sample_id}_qc_dedup_unmap_1.fastq.gz,${UNMAPPED_DIR}/${sample_id}_qc_dedup_unmap_2.fastq.gz" >> "${MANIFEST}"
done < "${SAMPLE_LIST}"

module load mags_maker/v1.1.1
mags_maker \
  --manifest "${MANIFEST}" \
  --qc skip \
  --outdir "${OUT_DIR}/mags_maker"
```

## Step 2 — qc_mags: taxonomy, decontamination, and re-assessment

A manifest is built from mags_maker's best-bin output per sample. qc_mags
produces, per sample, a `<Sample_ID>_final_report.tsv` (one row per MAG,
with both pre- and post-decontamination CheckM2/GUNC/QUAST/GTDB-Tk fields)
and a `<Sample_ID>_filter_summary.tsv`.

```bash
QC_MANIFEST="${OUT_DIR}/manifest_qc_mags.csv"
echo "ID,mags_dir" > "${QC_MANIFEST}"
find "${OUT_DIR}/mags_maker" -mindepth 2 -maxdepth 2 -type d -name "*_best_bins" | \
  awk -F'/' '{sample=$(NF-1); print sample","$0}' >> "${QC_MANIFEST}"

module load qc_mags/v1.0.0
qc_mags --manifest "${QC_MANIFEST}" --fasta_ext fasta
```

## Step 3 — load and merge every sample's QC report with participant metadata (R)

```r
library(tidyverse)
library(ggsci)
library(ggthemes)
library(ggpubr)

report_dir    <- "06_MAGS/qc_mags/final_reports"
filter_dir    <- "06_MAGS/qc_mags/filter_summaries"
metadata_file <- "databases/metadata_long.csv"

read_final_report <- function(file) {
  sample_id <- basename(file) %>% str_remove("_final_report.tsv$")
  read_tsv(file, show_col_types = FALSE) %>% mutate(Sample_ID = sample_id)
}
read_filter_summary <- function(file) {
  sample_id <- basename(file) %>% str_remove("_filter_summary.tsv$")
  read_tsv(file, show_col_types = FALSE) %>% mutate(Sample_ID = sample_id)
}

qc_final <- list.files(report_dir, pattern = "_final_report.tsv$", full.names = TRUE) %>%
  map_dfr(read_final_report)
qc_filter <- list.files(filter_dir, pattern = "_filter_summary.tsv$", full.names = TRUE) %>%
  map_dfr(read_filter_summary)

colnames(qc_final) <- colnames(qc_final) %>%
  str_remove_all("# ") %>%
  str_replace_all(" ", "_")

meta <- read_csv(metadata_file, show_col_types = FALSE)

qc_final_meta <- qc_final %>%
  left_join(meta, by = "Sample_ID") %>%
  mutate(
    timepoint = recode(timepoint, "pre" = "Pre-treatment", "post" = "Post-treatment"),
    timepoint = factor(timepoint, levels = c("Pre-treatment", "Post-treatment"))
  )
```

## Step 4 — classify MAG quality

High-quality: completeness ≥90%, contamination ≤5%, GUNC pass. Medium-
quality: completeness ≥50%, contamination ≤10%, GUNC pass. Everything
else is low-quality/fail. All thresholds are assessed on the
**post-decontamination** (post-MDMCleaner/seqkit) fields.

```r
qc_final_meta <- qc_final_meta %>%
  mutate(
    mag_quality = case_when(
      checkm2_postqc_Completeness >= 90 & checkm2_postqc_Contamination <= 5 &
        gunc_postqc_pass.GUNC == TRUE ~ "High-quality",
      checkm2_postqc_Completeness >= 50 & checkm2_postqc_Contamination <= 10 &
        gunc_postqc_pass.GUNC == TRUE ~ "Medium-quality",
      TRUE ~ "Low-quality / fail"
    ),
    pass_high_medium = mag_quality %in% c("High-quality", "Medium-quality"),
    pass_high = mag_quality == "High-quality"
  )

qc_final_meta$mag_quality <- factor(
  qc_final_meta$mag_quality,
  levels = c("High-quality", "Medium-quality", "Low-quality / fail")
)
```

## Step 5 — per-sample and per-treatment-group QC summary tables

```r
sample_summary <- qc_final_meta %>%
  group_by(Sample_ID, country, timepoint, treatment_arm) %>%
  summarise(
    n_mags_total = n(),
    n_high_quality = sum(mag_quality == "High-quality", na.rm = TRUE),
    n_medium_quality = sum(mag_quality == "Medium-quality", na.rm = TRUE),
    n_high_medium = sum(pass_high_medium, na.rm = TRUE),
    median_completeness = median(checkm2_postqc_Completeness, na.rm = TRUE),
    median_contamination = median(checkm2_postqc_Contamination, na.rm = TRUE),
    prop_gunc_pass = mean(gunc_postqc_pass.GUNC == TRUE, na.rm = TRUE),
    median_n50 = median(quast_postqc_N50, na.rm = TRUE),
    median_genome_size = median(checkm2_postqc_Genome_Size, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(sample_summary, "06_MAGS/qc_mags/MAG_QC_sample_summary.csv")

ttx_summary <- sample_summary %>%
  group_by(country, treatment_arm, timepoint) %>%
  summarise(
    n_samples = n(),
    total_mags = sum(n_mags_total, na.rm = TRUE),
    total_high_quality = sum(n_high_quality, na.rm = TRUE),
    total_high_medium = sum(n_high_medium, na.rm = TRUE),
    median_mags_per_sample = median(n_mags_total, na.rm = TRUE),
    median_high_quality_per_sample = median(n_high_quality, na.rm = TRUE),
    median_high_medium_per_sample = median(n_high_medium, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(ttx_summary, "06_MAGS/qc_mags/MAG_QC_treatment_summary.csv")
```

## Step 6 — the QC summary figure: MAGs per sample, completeness/contamination, GUNC pass rate, quality composition

```r
# Panel: number of MAGs generated per sample
mags_per_sample <- qc_final_meta %>%
  group_by(Sample_ID, country, treatment_arm, timepoint) %>%
  summarise(n_mags_total = n(), .groups = "drop")

p_mags_per_sample <- ggplot(mags_per_sample, aes(x = n_mags_total, color = timepoint, fill = timepoint)) +
  geom_density(alpha = 0.25) +
  facet_grid(country ~ treatment_arm, scales = "free_y") +
  scale_color_futurama(name = "") +
  scale_fill_futurama(name = "") +
  theme_bw() +
  labs(x = "Number of MAGs generated per sample", y = "Density") +
  theme(legend.position = 'top')

# Panel: CheckM2 completeness and contamination, pre- vs. post-decontamination
p_comp <- ggplot(qc_final_meta, aes(x = checkm2_preqc_Completeness, y = checkm2_postqc_Completeness, color = timepoint)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_grid(country ~ treatment_arm) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
  scale_color_futurama() +
  theme_bw() +
  labs(x = "Pre-QC completeness (%)", y = "Post-QC completeness (%)")

p_cont <- ggplot(qc_final_meta, aes(x = checkm2_preqc_Contamination, y = checkm2_postqc_Contamination, color = timepoint)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_grid(country ~ treatment_arm) +
  coord_cartesian(xlim = c(0, 25), ylim = c(0, 25)) +
  scale_color_futurama() +
  theme_bw() +
  labs(x = "Pre-QC contamination (%)", y = "Post-QC contamination (%)")

p_comp_cont <- ggarrange(p_comp, p_cont, ncol = 2, common.legend = TRUE, legend = "none")
top <- ggarrange(p_mags_per_sample, p_comp_cont, common.legend = TRUE, widths = c(1, 2))

# Panel: GUNC pass proportion, pre- vs. post-decontamination
gunc_prop <- qc_final_meta %>%
  group_by(country, treatment_arm) %>%
  summarise(
    pre_pass = sum(gunc_preqc_pass.GUNC == TRUE, na.rm = TRUE),
    pre_total = sum(!is.na(gunc_preqc_pass.GUNC)),
    post_pass = sum(gunc_postqc_pass.GUNC == TRUE, na.rm = TRUE),
    post_total = sum(!is.na(gunc_postqc_pass.GUNC)),
    .groups = "drop"
  ) %>%
  mutate(pre_prop = pre_pass / pre_total, post_prop = post_pass / post_total) %>%
  pivot_longer(cols = c(pre_prop, post_prop), names_to = "qc_stage", values_to = "prop_gunc_pass") %>%
  mutate(
    qc_stage = recode(qc_stage, pre_prop = "Pre-QC", post_prop = "Post-QC"),
    qc_stage = factor(qc_stage, levels = c("Pre-QC", "Post-QC")),
    n_pass = ifelse(qc_stage == "Pre-QC", pre_pass, post_pass),
    n_total = ifelse(qc_stage == "Pre-QC", pre_total, post_total),
    label = paste0(n_pass, "/", n_total, "\n", scales::percent(prop_gunc_pass, accuracy = 0.1))
  )

p_gunc_per <- ggplot(gunc_prop, aes(x = qc_stage, y = prop_gunc_pass, fill = qc_stage)) +
  geom_col(width = 0.7, color = "black") +
  geom_text(aes(label = label), vjust = -0.25, size = 3) +
  facet_grid(country ~ treatment_arm) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.08)) +
  scale_fill_tron(name = "QC stage") +
  theme_bw() +
  labs(x = NULL, y = "Proportion of MAGs passing GUNC") +
  theme(legend.position = "none")

gunc_total <- qc_final_meta %>%
  summarise(Pre_QC = sum(gunc_preqc_pass.GUNC == TRUE, na.rm = TRUE),
            Post_QC = sum(gunc_postqc_pass.GUNC == TRUE, na.rm = TRUE)) %>%
  pivot_longer(everything(), names_to = "qc_stage", values_to = "n_mags") %>%
  mutate(qc_stage = recode(qc_stage, Pre_QC = "Pre-QC", Post_QC = "Post-QC"),
         qc_stage = factor(qc_stage, levels = c("Pre-QC", "Post-QC")))

p_gunc_total <- ggplot(gunc_total, aes(x = qc_stage, y = n_mags, fill = qc_stage)) +
  geom_col(width = 0.7, color = "black") +
  geom_text(aes(label = n_mags), vjust = -0.3, size = 4) +
  scale_fill_tron() +
  coord_cartesian(ylim = c(0, 60000)) +
  theme_bw() +
  labs(x = NULL, y = "Total MAGs passing GUNC") +
  theme(legend.position = "none")

p_gunc <- ggarrange(p_gunc_total, p_gunc_per, widths = c(1, 2))

# Panel: MAG quality composition (donut)
quality_pie_df <- qc_final_meta %>%
  count(mag_quality) %>%
  mutate(
    mag_quality = factor(mag_quality, levels = c("High-quality", "Medium-quality", "Low-quality / fail"))
  ) %>%
  arrange(desc(mag_quality)) %>%
  mutate(
    prop = n / sum(n), ymax = cumsum(prop), ymin = lag(ymax, default = 0), ymid = (ymax + ymin) / 2,
    short_label = paste0(n, "\n", scales::percent(prop, accuracy = 0.1))
  )

p_quality_pie <- ggplot(quality_pie_df, aes(x = 1, y = prop, fill = mag_quality)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y", clip = "off") +
  geom_text(aes(y = ymid, label = short_label), x = 1, size = 4, color = "black", lineheight = 0.9, fontface = "bold") +
  scale_fill_ptol(name = "MAG quality") +
  theme_void() + theme(legend.position = 'bottom')

down <- ggarrange(p_gunc, p_quality_pie)
mags_qc <- ggarrange(top, down, ncol = 1)

ggsave('Figures/MAGs_qc.png', mags_qc, width = 16, height = 10)
ggsave('Figures/MAGs_qc.pdf', mags_qc, width = 16, height = 10)
```

## Step 7 — GTDB-Tk taxonomic assignment depth

```r
qc_final_meta <- qc_final_meta %>%
  mutate(
    tax_level = case_when(
      str_detect(gtdbtk_classification, "s__[^;]+") ~ "Species",
      str_detect(gtdbtk_classification, "g__[^;]+") ~ "Genus",
      str_detect(gtdbtk_classification, "f__[^;]+") ~ "Family",
      TRUE ~ "Higher / unclassified"
    )
  )

tax_summary <- qc_final_meta %>%
  group_by(tax_level) %>%
  summarise(n_mags = n(), n_unique_classifications = n_distinct(gtdbtk_classification), .groups = "drop")
write_csv(tax_summary, "06_MAGS/qc_mags/MAG_QC_taxonomic_assignment_summary.csv")
```

## Step 8 — select MAGs for downstream dereplication

Two selection tiers are written out: a looser set (≥50% completeness,
≤5% contamination, GUNC pass) used for an initial, more permissive
dereplication sensitivity check, and the final **high-quality (HQ)** set
(≥90% completeness, ≤5% contamination, GUNC pass) used throughout the
manuscript's primary analyses (22,032 MAGs). A third selection restricts
the HQ set to Enterobacteriaceae only, feeding the 99% ANI strain-level
dereplication two folders on.

```r
# Looser set (sensitivity check)
selected_mags <- qc_final_meta %>%
  filter(checkm2_postqc_Completeness >= 50, checkm2_postqc_Contamination <= 5, gunc_postqc_pass.GUNC == TRUE) %>%
  mutate(fasta_file = paste0(postqc_genome_name, ".fasta")) %>%
  select(Sample_ID, country, timepoint, treatment_arm, postqc_genome_name, fasta_file,
         checkm2_postqc_Completeness, checkm2_postqc_Contamination, gunc_postqc_pass.GUNC,
         gtdbtk_classification, tax_level)
write_tsv(selected_mags, "06_MAGS/qc_mags/selected_mags_50comp_5contam_guncpass.tsv")
write_lines(selected_mags$fasta_file, "06_MAGS/qc_mags/selected_mags_50comp_5contam_guncpass.txt")

# High-quality (HQ) set — the 22,032 MAGs used throughout the manuscript
hq_mags <- qc_final_meta %>%
  filter(checkm2_postqc_Completeness >= 90, checkm2_postqc_Contamination <= 5, gunc_postqc_pass.GUNC == TRUE) %>%
  mutate(fasta_file = paste0(postqc_genome_name, ".fasta")) %>%
  select(Sample_ID, country, timepoint, treatment_arm, postqc_genome_name, fasta_file,
         checkm2_postqc_Completeness, checkm2_postqc_Contamination, gunc_postqc_pass.GUNC,
         gtdbtk_classification, tax_level)
write_tsv(hq_mags, "06_MAGS/qc_mags/selected_HQ_mags.tsv")
write_lines(hq_mags$fasta_file, "06_MAGS/qc_mags/selected_HQ_mags.txt")

# genomeInfo.csv, required by dRep's --genomeInfo flag
genomeInfo_HQ <- hq_mags %>%
  transmute(genome = fasta_file, completeness = checkm2_postqc_Completeness, contamination = checkm2_postqc_Contamination)
write_csv(genomeInfo_HQ, "06_MAGS/drep/genomeInfo_HQ.csv")

# HQ MAGs assigned to family Enterobacteriaceae (feeds the 99% ANI dereplication)
enterobacteriaceae_HQ <- qc_final_meta %>%
  filter(checkm2_postqc_Completeness >= 90, checkm2_postqc_Contamination <= 5, gunc_postqc_pass.GUNC == TRUE,
         str_detect(gtdbtk_classification, "f__Enterobacteriaceae")) %>%
  mutate(fasta_file = paste0(postqc_genome_name, ".fasta")) %>%
  select(Sample_ID, country, timepoint, treatment_arm, postqc_genome_name, fasta_file,
         gtdbtk_classification, checkm2_postqc_Completeness, checkm2_postqc_Contamination)
write_tsv(enterobacteriaceae_HQ, "06_MAGS/drep/enterobacteriaceae_HQ_MAGs.tsv")
write_lines(enterobacteriaceae_HQ$fasta_file, "06_MAGS/drep/enterobacteriaceae_HQ_MAGs.txt")
```

This retains **22,032 high-quality MAGs**, matching the manuscript
Results, and feeds directly into the species-level (95% ANI) and
Enterobacteriaceae strain-level (99% ANI) dereplication steps in the
following folders. `qc_final_meta`, `hq_mags`, and `enterobacteriaceae_HQ`
are required inputs for both of those folders' R code — keep this R
session's objects (or re-run Steps 3–8) before proceeding.

