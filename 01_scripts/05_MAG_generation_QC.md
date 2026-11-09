# 5. Metagenome-assembled genome (MAG) generation and quality control

This stage assembles, bins, refines, taxonomically classifies, and
quality-filters MAGs from every sample's host-depleted reads, producing
the set of 22,032 high-quality MAGs that underlies all genome-resolved
analyses in the manuscript (Fig. 4, Fig. 5D, Extended Data Fig. 8–10).

Two purpose-built Nextflow pipelines are used in sequence:

- **mags_maker** (https://github.com/sanger-pathogens/mags_maker) —
  assembly (MEGAHIT), binning (MetaBAT2, MaxBin2, CONCOCT), and bin
  refinement (metaWRAP bin-refinement module), run independently per
  sample.
- **qc_mags** (https://github.com/sanger-pathogens/qc_mags) — taxonomic
  classification (GTDB-Tk), quality assessment (QUAST, CheckM2, GUNC),
  contaminant-contig removal (MDMCleaner), short-contig removal (seqkit),
  and a final re-assessment of the decontaminated bins.


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

## Step 2 — qc_mags: taxonomy and quality assessment

A manifest is built from mags_maker's best-bin output per sample.

```bash
QC_MANIFEST="${OUT_DIR}/manifest_qc_mags.csv"
echo "ID,mags_dir" > "${QC_MANIFEST}"
find "${OUT_DIR}/mags_maker" -mindepth 2 -maxdepth 2 -type d -name "*_best_bins" | \
  awk -F'/' '{sample=$(NF-1); print sample","$0}' >> "${QC_MANIFEST}"

module load qc_mags/v1.0.0
qc_mags --manifest "${QC_MANIFEST}" --fasta_ext fasta
```

qc_mags assesses bins both before and after running MDMCleaner
(contaminating-contig removal) and seqkit (short-contig removal), and
re-runs QUAST/CheckM2/GUNC on the decontaminated set. All downstream
analyses use the **post-decontamination** bins and their corresponding
quality metrics.

## Step 3 — selecting high-quality (HQ) MAGs

HQ MAGs are defined per MIMAG-consistent criteria (Methods): at least 90%
estimated completeness (CheckM2), no more than 5% estimated contamination
(CheckM2), and a GUNC chimerism pass — assessed on the post-decontamination
bins.

```r
library(tidyverse)

checkm2 <- read_tsv("06_MAGS/qc_mags/results/checkm2/quality_report.tsv", show_col_types = FALSE)
gunc    <- read_tsv("06_MAGS/qc_mags/results/gunc/GUNC.progenomes_2.1.maxCSS_level.tsv", show_col_types = FALSE)

hq_ids <- checkm2 %>%
  filter(Completeness >= 90, Contamination <= 5) %>%
  inner_join(gunc %>% filter(pass.GUNC == TRUE), by = c("Name" = "genome")) %>%
  pull(Name)

writeLines(paste0(hq_ids, ".fasta"), "06_MAGS/drep/selected_HQ_mags.txt")
```

This retains **22,032 high-quality MAGs** from the 60,738 MAGs initially
reconstructed across the 858 stool metagenomes, matching the counts
reported in the manuscript Results. This list (`selected_HQ_mags.txt`)
feeds directly into the species-level (95% ANI) dereplication step in the
next folder.
