################################################################################
# ALIVE microbiome study — shared analysis functions
#
# These functions are used identically across all three read-based data types
# (taxonomy, ARG/resistome, HUMAnN functional pathways). They were originally
# written independently for each data type and then consolidated here into one
# generalised set, per the manuscript's description of applying the same
# analytical framework (prevalence/abundance filtering, MaAsLin3 + ANCOM-BC2
# consensus differential abundance, alpha/beta diversity testing) across all
# three profiles.
#
################################################################################

library(phyloseq)
library(microbiome)
library(vegan)
library(tidyverse)
library(maaslin3)
library(ANCOMBC)

# ------------------------------------------------------------------------------
# 1. Prevalence/abundance filtering and compositional transform
#    (originally used for taxonomy; generalised here for ARG/pathway tables by
#    working directly on the OTU/feature table rather than re-aggregating by
#    taxonomic rank, which only applies to taxonomy)
# ------------------------------------------------------------------------------

#' Collapse a phyloseq object to a given taxonomic rank, convert to
#' compositional (relative abundance) data, and remove rare features.
#' Used for TAXONOMY only (phylum/family/genus/species aggregation).
#'
#' @param physeq_object phyloseq object (taxa_are_rows = FALSE), unfiltered
#' @param tax_rank one of "Phylum", "Family", "Genus", "Species"
#' @param detect minimum relative-abundance threshold (e.g. 0.5/100)
#' @param prev minimum prevalence threshold (e.g. 5/100)
create_compositional_phyloseq <- function(physeq_object, tax_rank, detect, prev) {
  otu_df <- as.data.frame(t(otu_table(physeq_object))) %>%
    rownames_to_column("OTU")
  tax_df <- as.data.frame(tax_table(physeq_object)) %>%
    rownames_to_column("OTU")
  meta_df <- sample_data(physeq_object) %>%
    as("data.frame") %>%
    rownames_to_column("SampleID")

  physeq_tidy <- otu_df %>%
    pivot_longer(cols = -OTU, names_to = "SampleID", values_to = "Percentage") %>%
    dplyr::left_join(tax_df, by = "OTU") %>%
    dplyr::left_join(meta_df, by = "SampleID")

  physeq_aggregated_abund <- physeq_tidy %>%
    drop_na(!!sym(tax_rank)) %>%
    group_by(SampleID, !!sym(tax_rank)) %>%
    dplyr::summarise(Abundance = sum(Percentage, na.rm = TRUE), .groups = "drop")

  otu_table_new <- physeq_aggregated_abund %>%
    pivot_wider(names_from = SampleID, values_from = Abundance, values_fill = 0) %>%
    column_to_rownames(tax_rank) %>%
    as.matrix()

  tax_table_new <- physeq_tidy %>%
    drop_na(!!sym(tax_rank)) %>%
    group_by(!!sym(tax_rank)) %>%
    dplyr::summarise(
      Kingdom = first(na.omit(Kingdom)),
      Phylum  = first(na.omit(Phylum)),
      Class   = first(na.omit(Class)),
      Order   = first(na.omit(Order)),
      Family  = first(na.omit(Family)),
      Genus   = first(na.omit(Genus)),
      Species = first(na.omit(Species))
    ) %>%
    ungroup()
  rownames(tax_table_new) <- tax_table_new[[tax_rank]]
  tax_table_new <- as.matrix(tax_table_new)

  metadata <- sample_data(physeq_object)

  physeq_new <- phyloseq(
    otu_table(otu_table_new, taxa_are_rows = TRUE),
    tax_table(tax_table_new),
    sample_data(metadata)
  )

  physeq_comp <- microbiome::transform(physeq_new, "compositional")
  physeq_final <- aggregate_rare(physeq_comp, level = tax_rank, detection = detect, prevalence = prev)
  return(physeq_final)
}

#' Filter an already-built phyloseq object (ARG or pathway table — no
#' taxonomic rank aggregation needed) by prevalence (fraction of samples in
#' which the feature is detected, i.e. abundance > 0).
#' Used for ARG and HUMAnN pathway tables.
#'
#' @param physeq_obj phyloseq object
#' @param min_prev minimum fraction of samples a feature must be detected in
filter_by_prevalence <- function(physeq_obj, min_prev = 0.05) {
  mat <- as(otu_table(physeq_obj), "matrix")
  if (!taxa_are_rows(physeq_obj)) mat <- t(mat)
  prev <- rowSums(mat > 0) / ncol(mat)
  keep <- names(prev[prev >= min_prev])
  prune_taxa(keep, physeq_obj)
}

# ------------------------------------------------------------------------------
# 2. Converting relative-abundance / normalised feature tables to pseudo-counts
#    for ANCOM-BC2 (which expects count-like data)
# ------------------------------------------------------------------------------

#' Convert a relative-abundance (or already-normalised) phyloseq object to
#' sequencing-depth-scaled pseudo-counts, by multiplying each sample's
#' relative abundances by that sample's total host-depleted read count and
#' rounding to the nearest integer. This keeps the total "read depth" per
#' sample consistent for ANCOM-BC2's internal bias-correction model.
#'
#' @param physeq_rel phyloseq object with relative-abundance (or RPKM/CPM) data
#' @param reads_df data.frame with columns `sample` and `Reads` (host-depleted
#'   read count per sample)
convert_rel_to_counts <- function(physeq_rel, reads_df) {
  mat_rel <- as(otu_table(physeq_rel), "matrix")
  if (!taxa_are_rows(physeq_rel)) mat_rel <- t(mat_rel)

  reads_vec <- reads_df$Reads
  names(reads_vec) <- reads_df$sample

  common_samples <- intersect(colnames(mat_rel), names(reads_vec))
  mat_rel <- mat_rel[, common_samples]
  reads_vec <- reads_vec[common_samples]

  mat_counts <- sweep(mat_rel, 2, reads_vec, `*`)
  mat_counts <- round(mat_counts)

  otu_new <- otu_table(mat_counts, taxa_are_rows = TRUE)
  phyloseq(otu_new, tax_table(physeq_rel), sample_data(physeq_rel))
}

#' Convert binary metadata columns encoded as 0/1 (or already "Yes"/"No") into
#' consistent two-level factors ("No","Yes"), in place on a phyloseq object's
#' sample_data.
convert_binary_vars <- function(ps_obj, vars) {
  meta_df <- as.data.frame(sample_data(ps_obj))
  for (v in vars) {
    if (v %in% colnames(meta_df)) {
      meta_df[[v]] <- ifelse(meta_df[[v]] %in% c(1, "1", "Yes"), "Yes", "No")
      meta_df[[v]] <- factor(meta_df[[v]], levels = c("No", "Yes"))
    }
  }
  sample_data(ps_obj) <- meta_df
  ps_obj
}

# ------------------------------------------------------------------------------
# 3. Alpha and beta diversity testing
#    (one generalised pair of functions used for taxonomy, ARG, and HUMAnN
#    Shannon diversity and Bray-Curtis PERMANOVA testing)
# ------------------------------------------------------------------------------

#' Test associations between Shannon alpha diversity and a set of numeric
#' (Spearman correlation) and categorical (Wilcoxon / Kruskal-Wallis)
#' variables, within a phyloseq object. Returns a tidy results data.frame.
#'
#' @param ps phyloseq object
#' @param num_vars character vector of numeric/continuous metadata columns
#' @param cat_vars character vector of categorical metadata columns
#' @param subgroup_label label for this subgroup (e.g. "Overall", "Ethiopia")
run_alpha_tests <- function(ps, num_vars, cat_vars, subgroup_label = "Overall") {
  alpha_div <- microbiome::alpha(ps, index = c("observed", "dominance_dbp", "diversity_shannon"))
  meta <- as.data.frame(sample_data(ps)) %>% rownames_to_column("Sample_ID")
  alpha_df <- data.frame(Sample_ID = rownames(alpha_div), alpha_div) %>%
    left_join(meta, by = "Sample_ID")

  num_res <- purrr::map_dfr(num_vars, function(v) {
    x <- as.numeric(alpha_df[[v]]); y <- as.numeric(alpha_df$diversity_shannon)
    if (length(unique(na.omit(x))) < 2) return(tibble())
    res <- suppressWarnings(cor.test(y, x, method = "spearman"))
    tibble(variable = v, test = "Spearman", statistic = as.numeric(res$estimate), p_value = res$p.value)
  })

  cat_res <- purrr::map_dfr(cat_vars, function(v) {
    y <- as.numeric(alpha_df$diversity_shannon); x <- alpha_df[[v]]
    x <- droplevels(factor(x))
    if (nlevels(x) < 2) return(tibble())
    if (nlevels(x) == 2) {
      res <- wilcox.test(y ~ x); test_name <- "Wilcoxon"
    } else {
      res <- kruskal.test(y ~ x); test_name <- "Kruskal-Wallis"
    }
    tibble(variable = v, test = test_name, statistic = as.numeric(res$statistic), p_value = res$p.value)
  })

  bind_rows(num_res, cat_res) %>%
    mutate(subgroup = subgroup_label, p_adj = p.adjust(p_value, method = "bonferroni"))
}

#' Test associations between Bray-Curtis beta diversity and a set of
#' variables using PERMANOVA (999 permutations).
run_beta_tests <- function(ps, num_vars, cat_vars, subgroup_label = "Overall") {
  dist_bc <- phyloseq::distance(ps, method = "bray")
  metadata <- as.data.frame(sample_data(ps))

  out <- purrr::map_dfr(c(num_vars, cat_vars), function(var) {
    if (!(var %in% colnames(metadata))) return(tibble())
    if (length(unique(na.omit(metadata[[var]]))) < 2) return(tibble())
    formula <- as.formula(paste("dist_bc ~", var))
    res <- adonis2(formula, data = metadata, permutations = 999)
    tibble(variable = var, test = "PERMANOVA", statistic = res$F[1], p_value = res$`Pr(>F)`[1])
  })

  out %>% mutate(subgroup = subgroup_label, p_adj = p.adjust(p_value, method = "bonferroni"))
}

# ------------------------------------------------------------------------------
# 4. Differential abundance: MaAsLin3 + ANCOM-BC2, baseline (cross-sectional)
#    comparisons. One generalised function, applied identically to taxonomy
#    (per-rank), ARG, and HUMAnN pathway tables.
# ------------------------------------------------------------------------------

#' Run the paired MaAsLin3 + ANCOM-BC2 differential-abundance framework for a
#' single cross-sectional comparison (e.g. country, coinfection status).
#'
#' For TAXONOMY, pass a list of already-rank-collapsed, prevalence/abundance-
#' filtered phyloseq objects (one per rank) via `phy_list`.
#' For ARG / HUMAnN pathways, pass a single already-prevalence-filtered
#' phyloseq object via `phy_list = list(feature = phy)`.
#'
#' @param phy_list named list of phyloseq objects (relative abundance / CPM /
#'   RPKM scale), one entry per "level" to test (e.g. phylum/family/genus/
#'   species for taxonomy, or a single "feature" entry for ARG/HUMAnN)
#' @param reads_df data.frame with `sample` and `Reads` columns, used to
#'   rescale to pseudo-counts for ANCOM-BC2
#' @param main_var the metadata variable being tested
#' @param tax_level_field the tax_table column name ANCOM-BC2 should use to
#'   label features (e.g. "Phylum", "Gene", "Feature")
#' @param prefix output subdirectory name for MaAsLin3 results
#' @param out_root base directory for MaAsLin3 output (differs by data type)
run_da_baseline <- function(phy_list, reads_df, main_var, tax_level_field = NULL,
                             prefix = "output", out_root = "MaasLin3/baseline") {

  maaslin_results <- list()
  ancom_results <- list()

  for (level in names(phy_list)) {
    phy <- phy_list[[level]]

    otu_df <- as.data.frame(t(otu_table(phy)))
    meta_df <- sample_data(phy) %>% as("data.frame")

    output_dir <- file.path(out_root, paste0(prefix, "_", level))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    set.seed(1)
    fit <- maaslin3(
      input_data = otu_df,
      input_metadata = meta_df,
      output = output_dir,
      formula = paste0("~ ", main_var),
      normalization = "TSS",
      transform = "LOG",
      augment = TRUE,
      standardize = TRUE,
      max_significance = 0.01,
      median_comparison_abundance = TRUE,
      median_comparison_prevalence = FALSE,
      max_pngs = 10,
      cores = 1
    )
    maaslin_results[[level]] <- fit

    phy_counts <- convert_rel_to_counts(phy, reads_df)
    tax_field <- if (!is.null(tax_level_field)) tax_level_field else level
    ancom_results[[level]] <- ancombc2(
      data = phy_counts,
      assay_name = "counts",
      tax_level = tax_field,
      fix_formula = main_var,
      rand_formula = NULL,
      p_adj_method = "bonferroni",
      pseudo = 0,
      alpha = 0.05
    )
  }

  list(maaslin = maaslin_results, ancom = ancom_results)
}

# ------------------------------------------------------------------------------
# 5. Differential abundance: MaAsLin3 + ANCOM-BC2, longitudinal (pre/post
#    treatment) comparisons, with participant ID as a random effect.
# ------------------------------------------------------------------------------

#' Run the paired MaAsLin3 + ANCOM-BC2 differential-abundance framework for a
#' longitudinal (paired pre/post treatment) comparison, with participant
#' identifier included as a random effect.
#'
#' @param phy_list named list of phyloseq objects, as in run_da_baseline()
#' @param reads_df data.frame with `sample` and `Reads` columns
#' @param main_var the metadata variable being tested (typically "timepoint")
#' @param co_var optional additional fixed-effect covariate (e.g. "country",
#'   when pooling both countries in one model)
#' @param random_effect_var participant identifier column for the random
#'   effect (default "Participant_ID")
run_da_longitudinal <- function(phy_list, reads_df, main_var, co_var = NULL,
                                 tax_level_field = NULL, random_effect_var = "Participant_ID",
                                 prefix = "output", out_root = "MaasLin3/treatment") {

  rand_formula_str <- paste0("(1|", random_effect_var, ")")

  fixed_terms <- main_var
  if (!is.null(co_var) && nzchar(co_var)) fixed_terms <- paste0(main_var, " + ", co_var)

  maaslin_formula <- paste0("~ ", fixed_terms, " + ", rand_formula_str)
  ancom_fix_formula <- fixed_terms
  ancom_rand_formula <- paste0("(1 | ", random_effect_var, ")")

  maaslin_results <- list()
  ancom_results <- list()

  for (level in names(phy_list)) {
    phy <- phy_list[[level]]

    otu_df <- as.data.frame(t(otu_table(phy)))
    meta_df <- sample_data(phy) %>% as("data.frame")

    output_dir <- file.path(out_root, paste0(prefix, "_", level))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    set.seed(1)
    fit <- maaslin3(
      input_data = otu_df,
      input_metadata = meta_df,
      output = output_dir,
      formula = maaslin_formula,
      normalization = "TSS",
      transform = "LOG",
      augment = TRUE,
      standardize = TRUE,
      max_significance = 0.01,
      median_comparison_abundance = TRUE,
      median_comparison_prevalence = FALSE,
      max_pngs = 10,
      cores = 1
    )
    maaslin_results[[level]] <- fit

    phy_counts <- convert_rel_to_counts(phy, reads_df)
    tax_field <- if (!is.null(tax_level_field)) tax_level_field else level
    ancom_results[[level]] <- ancombc2(
      data = phy_counts,
      assay_name = "counts",
      tax_level = tax_field,
      fix_formula = ancom_fix_formula,
      rand_formula = ancom_rand_formula,
      p_adj_method = "bonferroni",
      pseudo = 0,
      alpha = 0.05
    )
  }

  list(maaslin = maaslin_results, ancom = ancom_results)
}

# ------------------------------------------------------------------------------
# 6. Extracting significant features from MaAsLin3 / ANCOM-BC2 results
# ------------------------------------------------------------------------------

#' Read MaAsLin3's significant_results.tsv for one or more "levels" (ranks /
#' a single feature table) and combine into one tidy data.frame.
read_maaslin_results <- function(prefix, levels, out_root = "MaasLin3/baseline") {
  res_list <- list()
  for (lvl in levels) {
    file_path <- file.path(out_root, paste0(prefix, "_", lvl), "significant_results.tsv")
    if (file.exists(file_path)) {
      res_list[[lvl]] <- read_tsv(file_path, show_col_types = FALSE) %>%
        mutate(level = lvl)
    }
  }
  bind_rows(res_list)
}

#' Extract the ANCOM-BC2 results table (res$res) for one or more levels,
#' tagging each with its level name.
extract_ancom_results <- function(ancom_res_list) {
  purrr::map_dfr(names(ancom_res_list), function(level) {
    as.data.frame(ancom_res_list[[level]]$res) %>%
      rownames_to_column("feature") %>%
      mutate(level = level)
  })
}
