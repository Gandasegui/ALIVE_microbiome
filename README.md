# ALIVE Microbiome Study — Analysis Code

This repository contains the analysis code supporting:

> **Genome-resolved metagenomics reveals population-dependent gut microbiome responses to albendazole-ivermectin treatment** — Gandasegui, J. et al.

which characterises the gut microbiome, resistome, and functional profile
before and after anthelmintic treatment in the ALIVE phase 2/3 clinical trial
(ClinicalTrials.gov NCT05124691), in participants recruited in Ethiopia and
Mozambique.

Raw sequencing reads are deposited in the European Nucleotide Archive (ENA)
under BioProject **PRJEB79832**. Metagenome-assembled genomes are deposited
under the same BioProject (see Data Availability in the manuscript for the
specific tier/accession details). Sample-to-sequencing-file and
MAG-to-sample mappings are provided in Supplementary Tables 12–14.

## Repository structure

Each numbered folder corresponds to one stage of the analysis pipeline, and
contains a single `.md` walking through that analysis.
All R code chunks depend on the shared functions in
`R/functions/analysis_functions.R`.

| Folder | Stage | Key tools |
|---|---|---|
| `01_fastq_qc/` | Read QC and host-read removal | fastp v0.23.4, BBTools (Clumpify) v39.01, Bowtie2 v2.5.1 |
| `02_taxonomy/` | Taxonomic profiling and analysis | Sylph v0.8.1, sylph-tax v1.2.0, GlobDB r226, GTDB r226 |
| `03_resistome_ARG/` | Resistome (ARG) profiling and analysis | RGI bwt v6.0.8, CARD v4.0.1/v4.0.2 |
| `04_functional_HUMAnN/` | Functional (pathway) profiling and analysis | HUMAnN v3.8, MetaPhlAn v4.0.6 |
| `05_MAG_generation_QC/` | MAG assembly, binning, and quality control | mags_maker v1.1.1, qc_mags v1.0.0, CheckM2, GUNC, GTDB-Tk |
| `06_MAG_dereplication_species_95/` | Species-level (95% ANI) dereplication | dRep v3.5.0 |
| `07_MAG_dereplication_strain_99/` | Enterobacteriaceae strain-level (99% ANI) dereplication | dRep v3.5.0 |
| `08_MAG_characterisation_ARG/` | ARG annotation of MAGs, linked to taxonomy | merged_annotator (CARD/RGI + DeepARG) |

Prevalence/relative-abundance filter thresholds for differential-abundance
testing follow the manuscript Methods (0.5% relative abundance / 5%
prevalence for taxonomy; 5% prevalence for ARGs and pathways) and are set
as the defaults in `R/functions/analysis_functions.R`.

## Assumed inputs

All scripts assume the following are already available (built once, shared
across all analysis stages):

- `sample_list.txt` — one internal sample ID per line
- `databases/metadata_long.csv` — one row per Sample_ID, with country,
  timepoint (pre/post), treatment_arm, Participant_ID, and all baseline
  covariates described in Methods and Table 1
- `02_QC/read_counts.tsv` — per-sample read counts at each QC step (produced
  by `01_fastq_qc/run_qc.sh`)
- `databases/MAG_to_sample_to_rawfastq.csv` — Supplementary Table 13
  (MAG ID ↔ sample ID ↔ raw fastq / ENA accession mapping)

## Software environment

R v4.5.0, with: tidyverse, phyloseq, microbiome, vegan, maaslin3, ANCOMBC,
tableone, ggtree, ggtreeExtra, ape, ComplexUpset, patchwork, ggsci.

See each module's script header for the specific version of each
command-line tool used.

See `R/functions/analysis_functions.R` for the implementation, and the
manuscript Methods section for the full description.
