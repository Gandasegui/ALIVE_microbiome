# 7. Enterobacteriaceae strain-level (99% ANI) dereplication

This stage dereplicates all high-quality MAGs assigned to family
Enterobacteriaceae at a stricter 99% ANI threshold, to resolve *Escherichia
coli* and related taxa at the strain level, and computes strain-level
recovery prevalence across countries, treatment arms, and timepoints
(Fig. 4A–C, Extended Data Fig. 10 / Supplementary Fig. 3).

## Step 1 — collect Enterobacteriaceae HQ MAGs and run dRep at 99% ANI

`enterobacteriaceae_HQ_MAGs.txt` lists every HQ MAG (not just
representatives) assigned to family Enterobacteriaceae by GTDB-Tk — derived
from the taxonomy table built in the previous step.

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

echo "Strain-level cluster count:"
ls drep_HQ_enterobacteriaceae_099/dereplicated_genomes/*.fasta | wc -l
```

This yields **50 strain-level genomic clusters**, of which **31 belong to
*E. coli***.

## Step 2 — MAG-to-strain mapping and strain labelling

Strain clusters are labelled per species in order of cluster size (e.g.
"*Escherichia coli* s14", matching Fig. 4B's naming convention).

```r
library(tidyverse)

sample_metadata <- read_csv("databases/MAG_to_sample_to_rawfastq.csv", show_col_types = FALSE) # Supplementary Table 13
mag_species <- read_csv("06_MAGS/gtdbtk/gtdbtk_species_assignments.csv", show_col_types = FALSE) # MAG_ID, Species

cdb_99 <- read_csv("06_MAGS/drep/drep_HQ_enterobacteriaceae_099/data_tables/Cdb.csv", show_col_types = FALSE)
wdb_99 <- read_csv("06_MAGS/drep/drep_HQ_enterobacteriaceae_099/data_tables/Wdb.csv", show_col_types = FALSE)

mag_to_strain <- cdb_99 %>%
  left_join(wdb_99 %>% select(cluster = secondary_cluster, representative_genome = genome),
            by = c("secondary_cluster" = "cluster")) %>%
  transmute(
    MAG_ID = str_remove(genome, "\\.fasta$"),
    strain_cluster_rep = str_remove(representative_genome, "\\.fasta$")
  ) %>%
  left_join(mag_species, by = "MAG_ID") %>%
  left_join(sample_metadata %>% select(MAG_ID, Sample_ID), by = "MAG_ID") %>%
  left_join(sample_metadata %>% distinct(Sample_ID, country, treatment_arm, timepoint), by = "Sample_ID")

strain_labels <- mag_to_strain %>%
  distinct(Species, strain_cluster_rep) %>%
  count(Species, strain_cluster_rep, name = "n_dummy") %>%
  group_by(Species) %>%
  mutate(strain_label = paste0(Species, " s", row_number())) %>%
  ungroup() %>%
  select(-n_dummy)

mag_to_strain <- mag_to_strain %>% left_join(strain_labels, by = c("Species", "strain_cluster_rep"))

write_csv(mag_to_strain, "06_MAGS/drep/Enterobacteriaceae_MAG_to_strain_099.csv")
```

## Step 3 — strain-level recovery prevalence (Fig. 4C)

The denominator for prevalence is **all samples** in that
country/arm/timepoint group (from metadata), not only those with a
recovered MAG — so prevalence correctly reflects detection rate across the
whole cohort, not just among samples where any Enterobacteriaceae MAG was
assembled.

```r
sample_denominators <- sample_metadata %>%
  distinct(Sample_ID, country, treatment_arm, timepoint) %>%
  count(country, treatment_arm, timepoint, name = "n_samples_total")

strain_prevalence <- mag_to_strain %>%
  distinct(strain_label, Sample_ID, country, treatment_arm, timepoint) %>%
  count(strain_label, country, treatment_arm, timepoint, name = "n_samples_with_strain") %>%
  left_join(sample_denominators, by = c("country", "treatment_arm", "timepoint")) %>%
  mutate(prevalence = n_samples_with_strain / n_samples_total)

write_csv(strain_prevalence, "06_MAGS/drep/Enterobacteriaceae_strain_prevalence.csv")
```

## Step 4 — number of strain clusters per species (Supplementary Fig. 3A)

```r
strains_per_species <- mag_to_strain %>%
  distinct(Species, strain_label) %>%
  count(Species, name = "n_strain_clusters") %>%
  arrange(desc(n_strain_clusters))
```

## Step 5 — strain-level phylogeny (Supplementary Fig. 3B)

A phylogeny of the 50 representative strain-level genomes is built the
same way as the species-level trees in the previous folder: GTDB-Tk's
`identify`/`align` workflow extracts and aligns conserved marker genes,
and FastTree infers an approximate maximum-likelihood tree under the
Whelan–Goldman model.

```bash
set -euo pipefail
cd 06_MAGS/drep/drep_HQ_enterobacteriaceae_099

module load gtdbtk/2.7.1

gtdbtk identify \
  --genome_dir dereplicated_genomes \
  --out_dir gtdbtk_strain_099 \
  --extension fasta \
  --cpus 16

gtdbtk align \
  --identify_dir gtdbtk_strain_099 \
  --out_dir gtdbtk_strain_099 \
  --cpus 16

module load fasttree/2.1.10

FastTree -wag -gamma \
  < gtdbtk_strain_099/align/gtdbtk.bac120.user_msa.fasta \
  > gtdbtk_strain_099/strain_099.tree
```

The tree is then annotated with each representative's species assignment
and its country-specific MAG recovery prevalence (the outer heatmap ring
in Supplementary Fig. 3B), using the strain-to-sample mapping and
prevalence table built in Steps 2–3 above.

```r
library(ape)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)

tree_strain_099 <- read.tree("06_MAGS/drep/drep_HQ_enterobacteriaceae_099/gtdbtk_strain_099/strain_099.tree")

# One row per representative genome: its species, and its MAG recovery
# prevalence in each country (fraction of all samples in that country from
# which a MAG assigned to this strain cluster was recovered, regardless of
# treatment arm/timepoint — i.e. strain_prevalence from Step 3, collapsed
# to country level)
strain_country_prevalence <- strain_prevalence %>%
  group_by(strain_label, country) %>%
  summarise(
    n_samples_with_strain = sum(n_samples_with_strain),
    n_samples_total = sum(n_samples_total),
    prevalence = n_samples_with_strain / n_samples_total,
    .groups = "drop"
  )

rep_strain_info <- mag_to_strain %>%
  filter(MAG_ID == strain_cluster_rep) %>%
  distinct(MAG_ID, Species, strain_label)

p_strain_tree <- ggtree(tree_strain_099, layout = "circular") %<+% rep_strain_info +
  geom_tippoint(aes(color = Species), size = 2) +
  new_scale_fill() +
  geom_fruit(
    data = strain_country_prevalence,
    geom = geom_tile,
    mapping = aes(y = strain_label, x = country, fill = prevalence),
    width = 0.1,
    offset = 0.05
  ) +
  scale_fill_gradient(name = "MAG recovery\nprevalence", low = "white", high = "darkorange") +
  theme(legend.position = "right")

ggsave("06_MAGS/drep/drep_HQ_enterobacteriaceae_099/strain_tree_prevalence.pdf",
       p_strain_tree, width = 8, height = 8)
```

