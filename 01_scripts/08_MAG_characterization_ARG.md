# 8. Linking treatment-associated ARGs to their bacterial hosts

This final stage functionally annotates the high-quality MAGs for
antimicrobial resistance genes using multiple independent detection tools
(CARD/RGI, DeepARG, fARGene), harmonises their output with
**hAMRonization**, and combines confirmed ARG calls with GTDB-Tk taxonomy
to determine which bacterial species/family carries each treatment-
associated ARG identified in the read-based resistome analysis (module 3),
producing Fig. 4D and Fig. 5D.

## Step 1 — ARG annotation of MAGs (bash)

MAGs are functionally annotated using the merged_annotator Nextflow
workflow in metagenomic mode, with ARG annotation enabled (CARD/RGI,
DeepARG, fARGene all run; eggNOG-mapper skipped for runtime reasons).

```bash
set -euo pipefail

OUT_DIR="06_MAGS/annotator/all_HQ"
mkdir -p "${OUT_DIR}"

MANIFEST="${OUT_DIR}/manifest_HQ_mags.csv"
echo "ID,fasta" > "${MANIFEST}"
for f in 06_MAGS/drep/selected_HQ_mags/*.fasta; do
  id=$(basename "$f" .fasta)
  echo "${id},${f}" >> "${MANIFEST}"
done

module load nextflow/25.10.3
module load ISG/singularity/3.11.4

nextflow run merged_annotator/main.nf \
  --mode meta \
  --manifest "${MANIFEST}" \
  --outdir "${OUT_DIR}" \
  --func_annotate \
  --func_skip_eggnogmapper \
  --arg_annotate

# merged_annotator runs hAMRonization internally to harmonise the output of
# CARD/RGI, DeepARG, and fARGene into a single combined report.
```

## Step 2 — load and standardise the combined hAMRonization report (R)

```r
library(tidyverse)
library(janitor)
library(stringr)

ham_raw <- read_tsv("06_MAGS/annotator/all_HQ/hamronization_combined.tsv", show_col_types = FALSE)

ham <- ham_raw %>%
  mutate(
    mag_id = input_file_name %>%
      str_remove("\\.mapping\\.potential\\.ARG\\.deeparg$") %>%
      str_remove("\\.mapping\\.ARG\\.deeparg$") %>%
      str_remove("\\.txt\\.rgi$") %>%
      str_remove("_retrieved-genes-.*\\.out\\.fargene$") %>%
      str_remove("\\.txt\\.amrfinderplus$") %>%
      str_remove("\\.txt\\.abricate$"),
    tool = case_when(
      reference_database_name == "CARD" | analysis_software_name == "rgi" ~ "CARD_RGI",
      analysis_software_name == "deeparg" ~ "DeepARG",
      analysis_software_name == "fargene" ~ "fARGene",
      analysis_software_name == "amrfinderplus" ~ "AMRFinderPlus",
      analysis_software_name == "abricate" ~ "ABRICATE",
      TRUE ~ analysis_software_name
    ),
    gene_symbol_clean = na_if(str_squish(gene_symbol), ""),
    gene_name_clean = na_if(str_squish(gene_name), ""),
    drug_class_clean = na_if(str_squish(drug_class), ""),
    antimicrobial_agent_clean = na_if(str_squish(antimicrobial_agent), ""),
    resistance_mechanism_clean = na_if(str_squish(resistance_mechanism), ""),
    broad_class = coalesce(drug_class_clean, antimicrobial_agent_clean, resistance_mechanism_clean, gene_symbol_clean, gene_name_clean, "Unclassified"),
    arg_label = coalesce(gene_symbol_clean, gene_name_clean, antimicrobial_agent_clean, broad_class)
  )

write_tsv(ham, "06_MAGS/annotator/all_HQ/hamronization_all_hits_standardised.tsv")
```

## Step 3 — CARD/RGI hits, confirmed by cross-tool support

CARD/RGI is used as the primary ARG caller (it provides the detailed
drug-class/resistance-mechanism/ARO annotation used throughout the
manuscript). A CARD/RGI call is additionally flagged as supported if the
same MAG also carries a DeepARG and/or fARGene hit in the same broad
resistance class — this is the multi-tool confirmation referred to in
Methods.

```r
card_hits <- ham %>%
  filter(tool == "CARD_RGI") %>%
  transmute(
    mag_id, card_gene_symbol = gene_symbol_clean, card_gene_name = gene_name_clean,
    card_antimicrobial_agent = antimicrobial_agent_clean, card_drug_class = drug_class_clean,
    card_resistance_mechanism = resistance_mechanism_clean, card_aro = reference_accession,
    card_sequence_identity = sequence_identity, card_coverage_percentage = coverage_percentage,
    card_predicted_phenotype = predicted_phenotype, card_input_sequence_id = input_sequence_id
  )

deeparg_support <- ham %>%
  filter(tool == "DeepARG") %>%
  transmute(mag_id, deeparg_drug_class = drug_class_clean, deeparg_antimicrobial_agent = antimicrobial_agent_clean,
            deeparg_resistance_mechanism = resistance_mechanism_clean, deeparg_gene_symbol = gene_symbol_clean) %>%
  mutate(support_key = coalesce(deeparg_drug_class, deeparg_antimicrobial_agent, deeparg_resistance_mechanism, deeparg_gene_symbol)) %>%
  distinct()

fargene_support <- ham %>%
  filter(tool == "fARGene") %>%
  transmute(mag_id, fargene_gene_symbol = gene_symbol_clean, fargene_gene_name = gene_name_clean, fargene_input_file = input_file_name) %>%
  mutate(support_key = case_when(
    str_detect(str_to_lower(fargene_input_file), "tet_") ~ "tetracycline",
    TRUE ~ coalesce(fargene_gene_symbol, fargene_gene_name, "fargene_other")
  )) %>%
  distinct()

card_hits_supported <- card_hits %>%
  mutate(support_key = coalesce(card_drug_class, card_antimicrobial_agent, card_resistance_mechanism, card_gene_symbol)) %>%
  left_join(deeparg_support %>% group_by(mag_id) %>% summarise(has_deeparg = TRUE, .groups = "drop"), by = "mag_id") %>%
  left_join(fargene_support %>% group_by(mag_id) %>% summarise(has_fargene = TRUE, .groups = "drop"), by = "mag_id") %>%
  mutate(
    has_deeparg = replace_na(has_deeparg, FALSE),
    has_fargene = replace_na(has_fargene, FALSE),
    support_level = case_when(
      has_deeparg & has_fargene ~ "CARD + DeepARG + fARGene",
      has_deeparg ~ "CARD + DeepARG",
      has_fargene ~ "CARD + fARGene",
      TRUE ~ "CARD only"
    )
  ) %>%
  select(-support_key)

write_tsv(card_hits_supported, "06_MAGS/annotator/all_HQ/hamronization_card_supported_hits.tsv")
```

## Step 4 — MAG-level ARG burden: proportion of high-quality MAGs with CARD+DeepARG-confirmed ARGs

```r
# qc_final_meta is from module 5; sdb_095_clean (dRep Sdb.csv, 95% run) from module 6
mag_universe <- sdb_095_clean %>% distinct(mag_id = MAG_ID)

card_mags <- ham %>% filter(tool == "CARD_RGI") %>% distinct(mag_id)
deeparg_mags <- ham %>% filter(tool == "DeepARG") %>% distinct(mag_id)
highconf_mags <- inner_join(card_mags, deeparg_mags, by = "mag_id") %>% distinct(mag_id)

hq_mag_meta <- qc_final_meta %>%
  transmute(mag_id = postqc_genome_name, Sample_ID, country, timepoint, treatment_arm,
            participant_id = Participant_ID, gtdbtk_classification,
            completeness = checkm2_postqc_Completeness, contamination = checkm2_postqc_Contamination) %>%
  distinct()

mag_arg_hq_summary <- hq_mag_meta %>%
  mutate(has_card_deeparg = mag_id %in% highconf_mags$mag_id)

sample_burden <- mag_arg_hq_summary %>%
  group_by(Sample_ID, country, timepoint, treatment_arm) %>%
  summarise(n_hq_mags = n(), n_highconf_mags = sum(has_card_deeparg, na.rm = TRUE),
            prop_highconf_mags = n_highconf_mags / n_hq_mags, .groups = "drop")

p_burden <- ggplot(sample_burden, aes(x = timepoint, y = prop_highconf_mags, fill = timepoint)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
  geom_point(aes(color = timepoint), position = position_jitter(width = 0.12, height = 0), size = 1.5, alpha = 0.7) +
  facet_grid(country ~ treatment_arm) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_futurama() + scale_color_futurama() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Proportion of MAGs\nwith ARG presence")
```

## Step 5 — Fig. 4D: prevalence of the 15 FDCx3-associated ARGs among high-confidence Enterobacteriaceae MAGs

```r
library(scales)
library(viridis)
library(cowplot)

common_genes_eth_fdcx3 <- c(
  "RGI|ARO:3000216|acrB", "RGI|ARO:3000254|emrY", "RGI|ARO:3000491|acrD", "RGI|ARO:3000508|gadX",
  "RGI|ARO:3000793|mdtB", "RGI|ARO:3000794|mdtC", "RGI|ARO:3000833|evgS", "RGI|ARO:3001216|mdtH",
  "RGI|ARO:3001328|Escherichia coli mdfA", "RGI|ARO:3003577|ugd", "RGI|ARO:3003578|PmrF",
  "RGI|ARO:3003838|gadW", "RGI|ARO:3003950|msbA", "RGI|ARO:3003952|YojI", "RGI|ARO:3005044|OmpA"
)

target_tbl <- tibble(target_label = common_genes_eth_fdcx3) %>%
  mutate(aro = str_extract(target_label, "[0-9]+"), gene_label = str_remove(target_label, "^RGI\\|ARO:[0-9]+\\|"))

highconf_mag_ids <- mag_arg_hq_summary %>% filter(has_card_deeparg) %>% distinct(mag_id)

mag_tax <- hq_mag_meta %>%
  transmute(
    mag_id, Sample_ID, country, timepoint, treatment_arm, gtdbtk_classification,
    species = str_extract(gtdbtk_classification, "s__[^;]+") %>% str_remove("^s__"),
    genus = str_extract(gtdbtk_classification, "g__[^;]+") %>% str_remove("^g__"),
    family = str_extract(gtdbtk_classification, "f__[^;]+") %>% str_remove("^f__")
  ) %>%
  mutate(
    species = case_when(
      !is.na(species) & species != "" ~ species,
      !is.na(genus) & genus != "" ~ paste0(genus, "_sp"),
      !is.na(family) & family != "" ~ paste0(family, "_fam"),
      TRUE ~ "Unclassified"
    ),
    species = str_replace_all(species, " ", "_"),
    timepoint = factor(timepoint, levels = c("Pre-treatment", "Post-treatment")),
    treatment_arm = factor(treatment_arm, levels = c("ALB", "FDCx1", "FDCx3"))
  ) %>%
  distinct(mag_id, .keep_all = TRUE)

card_target_hits <- ham %>%
  filter(tool == "CARD_RGI") %>%
  semi_join(highconf_mag_ids, by = "mag_id") %>%
  mutate(reference_accession = as.character(reference_accession)) %>%
  inner_join(target_tbl, by = c("reference_accession" = "aro")) %>%
  distinct(mag_id, gene_label, .keep_all = TRUE) %>%
  left_join(mag_tax, by = "mag_id")

species_denominator <- mag_tax %>%
  filter(mag_id %in% highconf_mag_ids$mag_id) %>%
  distinct(mag_id, species, country, timepoint, treatment_arm) %>%
  count(species, country, timepoint, treatment_arm, name = "n_total_mags")

species_numerator <- card_target_hits %>%
  distinct(mag_id, gene_label, species, country, timepoint, treatment_arm) %>%
  count(gene_label, species, country, timepoint, treatment_arm, name = "n_positive_mags")

gene_species_prev <- species_numerator %>%
  left_join(species_denominator, by = c("species", "country", "timepoint", "treatment_arm")) %>%
  mutate(prevalence = n_positive_mags / n_total_mags, species_plot = str_replace_all(species, "_", " ")) %>%
  arrange(country, treatment_arm, timepoint, gene_label, desc(prevalence))

species_levels <- gene_species_prev %>%
  mutate(species_plot = fct_reorder(species_plot, prevalence, .fun = max, .desc = TRUE)) %>%
  pull(species_plot) %>% levels()

species_support <- mag_tax %>%
  filter(mag_id %in% highconf_mag_ids$mag_id) %>%
  mutate(species_plot = str_replace_all(species, "_", " ")) %>%
  distinct(mag_id, country, timepoint, treatment_arm, species_plot) %>%
  count(country, timepoint, treatment_arm, species_plot, name = "n_taxon_mags") %>%
  filter(species_plot %in% species_levels) %>%
  mutate(species_plot = factor(species_plot, levels = species_levels),
         timepoint = factor(timepoint, levels = c("Pre-treatment", "Post-treatment")))

# Heatmap: ARG prevalence by MAG species, pre/post FDCx3 (the regimen the 15 genes
# were originally identified as significantly decreased in)
p_gene_species_prev <- gene_species_prev %>%
  mutate(gene_label = factor(gene_label, levels = target_tbl$gene_label),
         species_plot = str_replace_all(species, "_", " "),
         species_plot = fct_reorder(species_plot, prevalence, .fun = max, .desc = TRUE)) %>%
  filter(treatment_arm == "FDCx3") %>%
  ggplot(aes(x = gene_label, y = species_plot, fill = prevalence)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_grid(country ~ treatment_arm + timepoint, scales = "free", space = "free_y") +
  scale_fill_gradientn(colors = c("grey92", viridis(100)), limits = c(0, 1),
                        breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("Absent", "25%", "50%", "75%", "100%"),
                        name = "Prevalence",
                        guide = guide_colorbar(title.position = "top", barheight = unit(4, "cm"), barwidth = unit(0.45, "cm"))) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey90", color = "black"), strip.text = element_text(size = 8),
        plot.margin = margin(8, 4, 8, 8), legend.position = "right") +
  labs(x = "", y = "MAG species")

# Left bar: number of high-confidence MAGs per species per facet (the support/denominator)
p_species_count_bar <- species_support %>%
  filter(treatment_arm == "FDCx3") %>%
  mutate(timepoint_short = factor(timepoint, levels = c("Pre-treatment", "Post-treatment"), labels = c("Pre", "Post"))) %>%
  ggplot(aes(x = n_taxon_mags, y = species_plot, fill = n_taxon_mags)) +
  geom_col(width = 0.75) +
  facet_grid(country ~ treatment_arm + timepoint_short, scales = "free_y", space = "free_y") +
  scale_fill_viridis_c(option = "C", guide = "none") +
  scale_x_continuous(breaks = c(0, 15, 30), expand = expansion(mult = c(0, 0.05))) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(),
        strip.background = element_rect(fill = "grey90", color = "black"), plot.margin = margin(8, 8, 8, 0), legend.position = "none") +
  labs(x = "Number of MAGs", y = NULL)

p_gene_species_prev_heatmap <- plot_grid(p_gene_species_prev, p_species_count_bar, ncol = 2, rel_widths = c(10, 2), align = "h", axis = "tb")

ggsave("Figures/MAGs_ARGs_entero.png", p_gene_species_prev_heatmap, width = 12, height = 4)
```

## Step 6 — Fig. 5D: prevalence of the four concordantly-increased tetracycline genes, across all bacterial families

Unlike Step 5, this is **not** restricted to Enterobacteriaceae, and each
CARD/RGI tet-gene call is confirmed only by DeepARG support at the
tetracycline-resistance class level (not fARGene, since fARGene's
tetracycline detection is organised by mechanism sub-type rather than
individual gene).

```r
library(tidytext)

target_tet_genes <- c("tet(W)", "tet(40)", "tet(O/32/O)", "tet(O/M/O)")

hq_mag_family <- qc_final_meta %>%
  filter(checkm2_postqc_Completeness >= 90, checkm2_postqc_Contamination <= 5, gunc_postqc_pass.GUNC == TRUE) %>%
  transmute(
    mag_id = postqc_genome_name, Sample_ID, country, timepoint, treatment_arm,
    family = str_extract(gtdbtk_classification, "f__[^;]+") %>% str_remove("^f__")
  ) %>%
  mutate(family = if_else(is.na(family) | family == "", "Unclassified family", family)) %>%
  distinct(mag_id, .keep_all = TRUE)

card_target_tet_hits <- ham %>%
  filter(tool == "CARD_RGI", gene_symbol_clean %in% target_tet_genes) %>%
  transmute(mag_id, gene = gene_symbol_clean, aro = as.character(reference_accession)) %>%
  distinct(mag_id, gene, .keep_all = TRUE) %>%
  semi_join(hq_mag_family, by = "mag_id")

deeparg_tet_mags <- ham %>%
  filter(tool == "DeepARG") %>%
  mutate(deeparg_search = str_c(coalesce(gene_symbol_clean, ""), coalesce(gene_name_clean, ""),
                                 coalesce(drug_class_clean, ""), coalesce(antimicrobial_agent_clean, ""),
                                 coalesce(resistance_mechanism_clean, ""), sep = " | ")) %>%
  filter(str_detect(deeparg_search, regex("tetracycline|\\btet(?:\\(|[_ -]?[a-z0-9])", ignore_case = TRUE))) %>%
  distinct(mag_id)

supported_tet_hits <- card_target_tet_hits %>%
  semi_join(deeparg_tet_mags, by = "mag_id") %>%
  inner_join(hq_mag_family, by = "mag_id") %>%
  distinct(mag_id, gene, family, .keep_all = TRUE)

family_denominator <- hq_mag_family %>% count(family, name = "n_HQ_MAGs")

family_gene_numerator <- supported_tet_hits %>%
  distinct(mag_id, family, gene) %>%
  count(family, gene, name = "n_gene_positive_MAGs")

families_with_target_gene <- supported_tet_hits %>% distinct(family) %>% pull(family)

family_gene_prevalence <- expand_grid(family = families_with_target_gene, gene = target_tet_genes) %>%
  left_join(family_gene_numerator, by = c("family", "gene")) %>%
  left_join(family_denominator, by = "family") %>%
  mutate(
    n_gene_positive_MAGs = replace_na(n_gene_positive_MAGs, 0L),
    prevalence = n_gene_positive_MAGs / n_HQ_MAGs,
    prevalence_plot = if_else(n_gene_positive_MAGs == 0, NA_real_, prevalence),
    tile_label = if_else(n_gene_positive_MAGs == 0, "", paste0(n_gene_positive_MAGs, "/", n_HQ_MAGs))
  )

# Fig. 5D heatmap: families carrying >=1 of the four genes, gene prevalence shown
p_tet_family_heatmap <- ggplot(family_gene_prevalence, aes(x = gene, y = family, fill = prevalence_plot)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = tile_label, colour = prevalence >= 0.40), size = 3.2, fontface = "bold") +
  scale_colour_manual(values = c("FALSE" = "#222222", "TRUE" = "white"), guide = "none") +
  scale_fill_gradient(low = "#D9EAE5", high = "#006D5B", na.value = "grey94", trans = "sqrt",
                       limits = c(0, 0.2), labels = percent_format(accuracy = 1), name = "Gene prevalence among HQ MAGs") +
  scale_y_discrete(drop = FALSE) +
  labs(x = NULL, y = "Bacterial family") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", colour = "black"),
        axis.text.y = element_text(colour = "black"), legend.position = "bottom")

# Companion bar plot: family-level prevalence, faceted by gene, families ordered
# independently within each gene facet
family_gene_bar_data <- family_gene_prevalence %>%
  filter(n_gene_positive_MAGs > 0) %>%
  mutate(
    gene = factor(gene, levels = target_tet_genes),
    family_plot = tidytext::reorder_within(as.character(family), prevalence, gene),
    count_label = paste0(n_gene_positive_MAGs, "/", n_HQ_MAGs)
  )

p_tet_family_barplot <- ggplot(family_gene_bar_data, aes(x = family_plot, y = prevalence)) +
  geom_col(fill = "#4F8A8B", colour = "grey25", linewidth = 0.2, width = 0.7) +
  geom_text(aes(label = count_label), vjust = -0.35, size = 3, colour = "black") +
  facet_grid(cols = vars(gene), scales = "free_x", space = "free_x") +
  tidytext::scale_x_reordered() +
  scale_y_continuous(limits = c(0, 0.1), breaks = c(0, 0.05, 0.1), labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0))) +
  labs(x = "Bacterial family", y = "Prevalence among high-quality MAGs") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey92", colour = "grey30"),
        strip.text = element_text(face = "italic", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"))
```

