# 1. Raw read quality control and host-read removal

This stage takes paired-end fastq.gz files (already retrieved and renamed
to `<sample_id>_1.fastq.gz` / `<sample_id>_2.fastq.gz`; see Supplementary
Table 12 for the sample-to-raw-fastq mapping and ENA accessions) and
produces host-depleted, quality-filtered reads ready for taxonomic,
resistome, functional, and MAG analysis.

Three tools are used in sequence: **fastp** (adapter trimming / quality
filtering), **Clumpify** from BBTools (PCR/optical duplicate removal), and
**Bowtie2** (human host-read removal against GRCh38/hg38). Read counts are
tracked at every step, which is both reported in the manuscript
(Supplementary Fig. 1A) and used to apply the pre-specified exclusion
criterion of a minimum of 20 million host-depleted read pairs per sample.

Software versions: fastp v0.23.4, BBTools v39.01, Bowtie2 v2.5.1.

We deliberately do not include the step that retrieves raw fastq files
from our institution's iRODS storage system here, since that step is
specific to internal Sanger Institute infrastructure and is of no use to
anyone reproducing this pipeline elsewhere. The raw reads are available
directly from ENA under BioProject **PRJEB79832**.

## Step 1 — fastp

Adapter trimming and quality filtering: a Phred quality score of at least
20 is required per base; reads are discarded if more than 30% of their
bases fail this threshold; reads with more than five ambiguous bases or
shorter than 50 bp after trimming are removed; poly-G tails are trimmed;
and overlapping paired reads are used for base correction.

```bash
set -euo pipefail

RAW_DIR="01_RAW/renamed"          # contains <sample_id>_1/2.fastq.gz
QC_DIR="02_QC"
SAMPLE_LIST="sample_list.txt"     # one sample_id per line

mkdir -p "${QC_DIR}/fastp"

module load fastp/0.23.4--hadf994f_2

while read -r sample_id; do
  echo "Fastq QC for ${sample_id} started"

  R1="${RAW_DIR}/${sample_id}_1.fastq.gz"
  R2="${RAW_DIR}/${sample_id}_2.fastq.gz"

  fastp \
    -i "${R1}" -I "${R2}" \
    -o "${QC_DIR}/fastp/${sample_id}_qc_1.fastq.gz" \
    -O "${QC_DIR}/fastp/${sample_id}_qc_2.fastq.gz" \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 30 \
    --n_base_limit 5 \
    --length_required 50 \
    --trim_poly_g \
    --correction \
    --json "${QC_DIR}/fastp/${sample_id}_fastp.json" \
    --html "${QC_DIR}/fastp/${sample_id}_fastp.html"

done < "${SAMPLE_LIST}"
```

## Step 2 — Clumpify (BBTools): PCR/optical duplicate removal

Run with an optical duplicate distance of 12,000, two deduplication passes,
and no substitutions permitted.

```bash
module load bbtools/39.01
mkdir -p "${QC_DIR}/dedup"

while read -r sample_id; do
  clumpify.sh \
    in1="${QC_DIR}/fastp/${sample_id}_qc_1.fastq.gz" \
    in2="${QC_DIR}/fastp/${sample_id}_qc_2.fastq.gz" \
    out1="${QC_DIR}/dedup/${sample_id}_qc_dedup_1.fastq.gz" \
    out2="${QC_DIR}/dedup/${sample_id}_qc_dedup_2.fastq.gz" \
    dedupe=t \
    optical=t \
    dupedist=12000 \
    passes=2 \
    subs=0
done < "${SAMPLE_LIST}"
```

## Step 3 — Bowtie2: human host-read removal

Reads are mapped against the GRCh38/hg38 reference in `--very-sensitive`
mode; unmapped (non-human) read pairs are extracted and retained for all
downstream analysis.

```bash
module load bowtie2/2.5.1--py38he00c5e5_2
module load samtools/1.9--h91753b0_8

HUMAN_REF="human_ref/hg38.fa"     # Bowtie2 index prefix, GRCh38/hg38
mkdir -p "${QC_DIR}/unmapped"

while read -r sample_id; do
  bowtie2 \
    -x "${HUMAN_REF}" \
    -1 "${QC_DIR}/dedup/${sample_id}_qc_dedup_1.fastq.gz" \
    -2 "${QC_DIR}/dedup/${sample_id}_qc_dedup_2.fastq.gz" \
    --very-sensitive \
    --threads 8 \
    -S "${QC_DIR}/${sample_id}_qc_dedup_map.sam"

  # Extract unmapped (non-human) read pairs
  samtools view -bS -f 12 -F 256 "${QC_DIR}/${sample_id}_qc_dedup_map.sam" \
    | samtools sort -n -o "${QC_DIR}/${sample_id}_qc_dedup_unmap.bam" -
  samtools fastq \
    -1 "${QC_DIR}/unmapped/${sample_id}_qc_dedup_unmap_1.fastq.gz" \
    -2 "${QC_DIR}/unmapped/${sample_id}_qc_dedup_unmap_2.fastq.gz" \
    "${QC_DIR}/${sample_id}_qc_dedup_unmap.bam"

  rm -f "${QC_DIR}/${sample_id}_qc_dedup_map.sam" "${QC_DIR}/${sample_id}_qc_dedup_unmap.bam"
  echo "${sample_id} ended"
done < "${SAMPLE_LIST}"
```

## Step 4 — read-count tracking and the 20-million-read-pair exclusion rule

Read counts are recorded at every QC step (input → fastp → Clumpify →
Bowtie2), producing `02_QC/read_counts.tsv`, which both feeds
Supplementary Fig. 1A and is used to apply the manuscript's pre-specified
exclusion criterion: *"All the participants with a read count lower than
20 million pairs of reads were excluded."*

```bash
OUT_COUNTS="${QC_DIR}/read_counts.tsv"
echo -e "sample\tinput_reads\tfastp_qc_reads\tclumpify_dedup_reads\tbowtie_hostrm_reads" > "${OUT_COUNTS}"

while read -r sample_id; do
  input_reads=$(( $(zcat "${RAW_DIR}/${sample_id}_1.fastq.gz" | wc -l) / 4 ))
  fastp_reads=$(( $(zcat "${QC_DIR}/fastp/${sample_id}_qc_1.fastq.gz" | wc -l) / 4 ))
  dedup_reads=$(( $(zcat "${QC_DIR}/dedup/${sample_id}_qc_dedup_1.fastq.gz" | wc -l) / 4 ))
  hostrm_reads=$(( $(zcat "${QC_DIR}/unmapped/${sample_id}_qc_dedup_unmap_1.fastq.gz" | wc -l) / 4 ))
  echo -e "${sample_id}\t${input_reads}\t${fastp_reads}\t${dedup_reads}\t${hostrm_reads}" >> "${OUT_COUNTS}"
done < "${SAMPLE_LIST}"

# Exclude samples with <20,000,000 host-depleted read pairs
awk -F'\t' 'NR>1 && $5 < 20000000 {print $1}' "${OUT_COUNTS}" > "${QC_DIR}/samples_excluded_low_depth.txt"
echo "Samples excluded for <20M host-depleted read pairs:"
cat "${QC_DIR}/samples_excluded_low_depth.txt"
```

`02_QC/read_counts.tsv` is a required input for every downstream analysis
module in this repository (used both for QC reporting and for rescaling
relative-abundance data to pseudo-counts ahead of ANCOM-BC2 testing).

