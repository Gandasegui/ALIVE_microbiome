# 6. Species-level (95% ANI) MAG dereplication

This stage dereplicates the 22,032 high-quality MAGs into species-level
representative genomes, then characterises their taxonomy, country
distribution, and phylogeny (Fig. 4A, Extended Data Fig. 8–9).

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
done < selected_HQ_mags.txt

find "$PWD/selected_HQ_mags" -type f -name "*.fasta" | sort > selected_HQ_mags_paths.txt
wc -l selected_HQ_mags_paths.txt
```

`genomeInfo_HQ.csv` (dRep's required `--genomeInfo` input: a two-column CSV
of `genome,completeness,contamination`) is generated directly from the
CheckM2 quality report written in step 5:

```r
library(tidyverse)
checkm2 <- read_tsv("../qc_mags/results/checkm2/quality_report.tsv", show_col_types = FALSE)
checkm2 %>%
  transmute(genome = paste0(Name, ".fasta"), completeness = Completeness, contamination = Contamination) %>%
  write_csv("genomeInfo_HQ.csv")
```

```bash
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

dRep writes, among other outputs: `data_tables/Cdb.csv` (genome-to-cluster
assignment), `data_tables/Wdb.csv` (the chosen representative per cluster),
and `dereplicated_genomes/` (the **1,372 representative genome fastas**
reported in the manuscript).

## Step 2 — MAG-to-representative mapping (Supplementary Table 14)

```r
library(tidyverse)

sample_metadata <- read_csv("databases/MAG_to_sample_to_rawfastq.csv", show_col_types = FALSE) # Supplementary Table 13

cdb <- read_csv("06_MAGS/drep/drep_HQ_095/data_tables/Cdb.csv", show_col_types = FALSE)
wdb <- read_csv("06_MAGS/drep/drep_HQ_095/data_tables/Wdb.csv", show_col_types = FALSE)

mag_to_rep <- cdb %>%
  left_join(wdb %>% select(cluster = secondary_cluster, representative_genome = genome),
            by = c("secondary_cluster" = "cluster")) %>%
  transmute(
    MAG_ID = str_remove(genome, "\\.fasta$"),
    Representative_MAG_ID_095 = str_remove(representative_genome, "\\.fasta$")
  ) %>%
  left_join(sample_metadata %>% select(MAG_ID, Sample_ID), by = "MAG_ID")

write_csv(mag_to_rep, "06_MAGS/drep/MAG_to_representative_095.csv") # -> Supplementary Table 14

n_representatives <- mag_to_rep %>% filter(MAG_ID == Representative_MAG_ID_095) %>% nrow()
message("Number of representative (species-level) genomes: ", n_representatives) # expect 1,372
```

## Step 3 — GTDB-Tk taxonomy, country pattern, and novel-species summary (Extended Data Fig. 8–9)

Each representative's "country pattern" (shared / Ethiopia-only /
Mozambique-only) is determined by which country(ies) contributed any MAG
assigned to its cluster.

```r
gtdbtk_bac <- read_tsv("06_MAGS/gtdbtk/gtdbtk.bac120.summary.tsv", show_col_types = FALSE)
gtdbtk_ar  <- read_tsv("06_MAGS/gtdbtk/gtdbtk.ar53.summary.tsv",  show_col_types = FALSE)
gtdbtk_all <- bind_rows(gtdbtk_bac, gtdbtk_ar) %>% rename(MAG_ID = user_genome)

representative_ids <- mag_to_rep %>% filter(MAG_ID == Representative_MAG_ID_095) %>% pull(MAG_ID)

rep_taxonomy <- gtdbtk_all %>%
  filter(MAG_ID %in% representative_ids) %>%
  separate(classification,
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", fill = "right") %>%
  mutate(across(Domain:Species, ~ str_remove(., "^[a-z]__")))

country_pattern <- mag_to_rep %>%
  left_join(sample_metadata %>% select(Sample_ID, country), by = "Sample_ID") %>%
  group_by(Representative_MAG_ID_095) %>%
  summarise(
    countries = paste(sort(unique(country)), collapse = ";"),
    cluster_size = n(),
    .groups = "drop"
  ) %>%
  mutate(country_pattern = case_when(
    str_detect(countries, ";") ~ "Shared",
    countries == "Ethiopia" ~ "Ethiopia-only",
    countries == "Mozambique" ~ "Mozambique-only",
    TRUE ~ NA_character_
  ))

rep_summary <- rep_taxonomy %>%
  left_join(country_pattern, by = c("MAG_ID" = "Representative_MAG_ID_095"))

write_csv(rep_summary, "06_MAGS/drep/representative_genomes_summary.csv")

# Under-classified (putative novel) species/genus — GTDB-Tk leaves the
# Species (or Genus) field blank/placeholder when a representative cannot
# be confidently assigned to a named taxon
novel_species <- rep_summary %>% filter(is.na(Species) | Species == "")
novel_genus   <- rep_summary %>% filter(is.na(Genus) | Genus == "")

message("Representative genomes not assigned to a named species: ", nrow(novel_species))
message("Of which, not assigned even to a named genus: ", nrow(novel_genus))
```

## Step 4 — bacterial and archaeal phylogenies (Extended Data Fig. 8)

Trees are built upstream with GTDB-Tk identify/align + FastTree v2.1.10
(Whelan–Goldman model); this step only loads and annotates them.

```r
library(ape)
library(ggtree)
library(ggtreeExtra)

tree_bac <- read.tree("06_MAGS/gtdbtk/gtdbtk.bac120.classify.tree")
tree_ar  <- read.tree("06_MAGS/gtdbtk/gtdbtk.ar53.classify.tree")

p_bac <- ggtree(tree_bac, layout = "circular") %<+% rep_summary +
  geom_fruit(geom = geom_tile, mapping = aes(fill = Phylum), width = 0.02) +
  new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping = aes(fill = country_pattern), width = 0.02) +
  new_scale_fill() +
  geom_fruit(geom = geom_bar, mapping = aes(y = MAG_ID, x = cluster_size),
             stat = "identity", orientation = "y", width = 0.3)

p_ar <- ggtree(tree_ar) %<+% rep_summary +
  geom_tiplab(aes(color = country_pattern)) +
  geom_tippoint(aes(size = cluster_size))
```
