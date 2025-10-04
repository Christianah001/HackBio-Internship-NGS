# Differential Gene Expression Analysis of Arabidopsis Vasculature Under UV-C Stress

**Author:**  
Funmilayo C. Ligali  
---

## Introduction

Berkowitz et al. (2021) studied tissue-specific responses of *Arabidopsis thaliana* leaves to abiotic stresses. In this project, we focus on the **vasculature tissue** and determine genes responding to UV-C treatment compared to water-treated control samples. RNA-seq data from three biological replicates per condition were analyzed.

UV-C is a high-energy radiation that induces DNA damage and stress response pathways. Understanding the transcriptional response in vasculature helps uncover stress adaptation mechanisms.

---

## Objectives

- Download and process raw RNA-seq data (FASTQ files) for control and UV-C treated samples.  
- Perform quality control, trimming, and alignment of reads to the *Arabidopsis thaliana* genome (TAIR10).  
- Quantify gene expression using **featureCounts**.  
- Perform differential expression analysis between UV-C treated and control samples using **DESeq2**.  
- Identify top 100 differentially expressed genes.  
- Conduct **GO** and **KEGG** enrichment analysis to determine significantly affected biological pathways.

---

## Materials & Data

| Replicate | Control SRR   | UV-C Treated SRR |
|-----------|---------------|----------------|
| 1         | SRR12808527   | SRR12808497    |
| 2         | SRR12808528   | SRR12808498    |
| 3         | SRR12808529   | SRR12808499    |

**Reference genome:** *Arabidopsis thaliana* TAIR10 (FASTA & GFF3 from Ensembl Plants)

## Methods
## RNA-seq Analysis Pipeline

| Step | Tool / Software | Input | Output | Description |
|------|----------------|-------|--------|-------------|
| 1 | Bash | SRR IDs | Raw FASTQ files | Download raw RNA-seq data from ENA |
| 2 | FastQC / MultiQC | Raw FASTQ | QC reports | Assess quality of raw reads |
| 3 | fastp | Raw FASTQ | Trimmed FASTQ, HTML/JSON report | Adapter and quality trimming |
| 4 | FastQC / MultiQC | Trimmed FASTQ | QC reports | Assess quality post-trimming |
| 5 | STAR | Trimmed FASTQ, TAIR10 genome | Sorted BAM files | Align reads to Arabidopsis genome |
| 6 | Samtools | BAM files | Indexed BAM | Index aligned BAM for downstream use |
| 7 | featureCounts | BAM, GFF3 annotation | Counts table | Gene-level read quantification |
| 8 | R / DESeq2 | Counts table, metadata | DESeq2 results | Differential expression analysis |
| 9 | pheatmap (R) | DESeq2 results | Heatmaps | Visualize differentially expressed genes |
| 10 | clusterProfiler | Top DE genes | GO enrichment | Identify significantly enriched pathways |

### Bash Scripts
# RNA-seq Preprocessing and Quantification Pipeline

**Author:** Funmilayo C. Ligali  
**Date:** 2025-10-04  

This pipeline automates **RNA-seq preprocessing** from raw FASTQ downloads to gene-level quantification, specifically for Arabidopsis vasculature UV-C stress study.

---

## Pipeline Overview

This pipeline performs the following steps:

1. Create project directories.  
2. Download raw FASTQ files from ENA using SRR IDs.  
3. Quality control with `FastQC` and `MultiQC`.  
4. Trim adapters and low-quality bases using `fastp`, followed by QC.  
5. Download the reference genome (TAIR10) and annotation (GFF3) and generate STAR genome index.  
6. Align trimmed reads to the genome with `STAR` and index BAM files with `samtools`.  
7. Quantify gene-level counts using `featureCounts`.

---

## Bash Pipeline
# Complete RNA-seq preprocessing and quantification pipeline

0. All project files and outputs were organized into directories
```bash
#!/bin/bash
# rna_seq_pipeline.sh

# ---------------------------
# 0. Set Directories
# ---------------------------
BASE_DIR=$(pwd)
RAW_DIR="$BASE_DIR/raw_fastq"
TRIM_DIR="$BASE_DIR/trimmed_fastq"
QC_DIR="$BASE_DIR/qc_results"
QC_TRIM="$QC_DIR/trimmed"
REF_DIR="$BASE_DIR/reference"
GENOME_DIR="$REF_DIR/genome"
GTF_DIR="$REF_DIR/annotation"
STAR_INDEX_DIR="$REF_DIR/STAR_index"
BAM_DIR="$BASE_DIR/bam"
COUNT_DIR="$BASE_DIR/counts"
LOG_DIR="$BASE_DIR/logs"

mkdir -p "$RAW_DIR" "$TRIM_DIR" "$QC_DIR" "$QC_TRIM" "$GENOME_DIR" \
         "$GTF_DIR" "$STAR_INDEX_DIR" "$BAM_DIR" "$COUNT_DIR" "$LOG_DIR"
```
1. Download Raw RNA-seq FASTQ Files

SRR IDs for control and UV-C treated samples were downloaded from ENA using FTP links.
```bash
# ---------------------------
# 1. Download SRR FASTQs
# ---------------------------
cat > "$BASE_DIR/srrs.txt" <<'EOF'
SRR12808527
SRR12808528
SRR12808529
SRR12808497
SRR12808498
SRR12808499
EOF

echo "=== Downloading FASTQ files ==="
while read SRR; do
    echo "Processing $SRR..."
    curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$SRR&result=read_run&fields=run_accession,fastq_ftp,library_layout" \
    | tail -n +2 | while IFS=$'\t' read -r run_accession fastq_ftp library_layout; do
        for url in $(echo $fastq_ftp | tr ';' ' '); do
            fname=$(basename $url)
            echo "Downloading $fname ($library_layout)"
            curl -L "ftp://$url" -o "$RAW_DIR/$fname"
        done
    done
done < "$BASE_DIR/srrs.txt"
echo "=== FASTQ download completed ==="
```
2. Quality Control (FastQC + MultiQC)

Raw FASTQ files were assessed for quality before trimming.
```bash
# ---------------------------
# 2. QC raw FASTQs
# ---------------------------
echo "=== Running FastQC on raw FASTQs ==="
for fq in "$RAW_DIR"/*.fastq.gz; do
    sample=$(basename "$fq")
    fastqc "$fq" --outdir "$QC_DIR" >"$LOG_DIR/${sample%.fastq.gz}_fastqc.log" 2>&1
done

echo "=== MultiQC summary ==="
multiqc "$QC_DIR" --outdir "$QC_DIR" >"$LOG_DIR/multiqc.log" 2>&1
```
3. Trimming (fastp) and QC

Reads were trimmed to remove adapters and low-quality bases, then re-checked with FastQC.
```bash
# ---------------------------
# 3. Trim with fastp + QC
# ---------------------------
echo "=== Trimming raw FASTQs ==="
for fq in "$RAW_DIR"/*.fastq.gz; do
    base=$(basename "$fq" .fastq.gz)
    fastp -i "$fq" -o "$TRIM_DIR/${base}_trim.fastq.gz" \
          -h "$TRIM_DIR/${base}_fastp.html" \
          -j "$TRIM_DIR/${base}_fastp.json" \
          >>"$LOG_DIR/fastp.log" 2>&1
done

echo "=== FastQC on trimmed FASTQs ==="
for fq in "$TRIM_DIR"/*.fastq.gz; do
    base=$(basename "$fq" .fastq.gz)
    fastqc "$fq" --outdir "$QC_TRIM" >>"$LOG_DIR/fastqc_trimmed.log" 2>&1
done

multiqc "$QC_TRIM" --outdir "$QC_TRIM" >>"$LOG_DIR/multiqc_trimmed.log" 2>&1
echo "=== Trimming + QC completed ==="
```
4. Reference Genome Download and STAR Index

The Arabidopsis TAIR10 genome and GFF3 annotation were downloaded from Ensembl Plants, then indexed for STAR alignment.
```bash
# ---------------------------
# 4. Download reference genome + STAR index
# ---------------------------
echo "=== Downloading Arabidopsis genome & annotation ==="
wget -c https://ftp.ensemblgenomes.org/pub/plants/release-62/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz \
     -O "$GENOME_DIR/TAIR10.fa.gz"
wget -c https://ftp.ensemblgenomes.org/pub/plants/release-62/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.62.gff3.gz \
     -O "$GTF_DIR/TAIR10.gff3.gz"
gunzip -f "$GENOME_DIR/TAIR10.fa.gz"
gunzip -f "$GTF_DIR/TAIR10.gff3.gz"

echo "=== Generating STAR genome index ==="
STAR --runThreadN 12 \
     --runMode genomeGenerate \
     --genomeDir "$STAR_INDEX_DIR" \
     --genomeFastaFiles "$GENOME_DIR/TAIR10.fa" \
     --sjdbGTFfile "$GTF_DIR/TAIR10.gff3" \
     --sjdbOverhang 99
echo "=== STAR genome index ready ==="
```
5. STAR Alignment

Trimmed reads were aligned to TAIR10 with STAR. BAM files were sorted and indexed.
```bash
# ---------------------------
# 5. STAR Alignment
# ---------------------------
echo "=== Aligning trimmed FASTQs with STAR ==="
for fq in "$TRIM_DIR"/*_trim.fastq.gz; do
    base=$(basename "$fq" _trim.fastq.gz)
    bam_file="$BAM_DIR/${base}_Aligned.sortedByCoord.out.bam"
    [ -s "$bam_file" ] && { echo "$base BAM exists, skipping."; continue; }
    
    STAR --runThreadN 12 \
         --genomeDir "$STAR_INDEX_DIR" \
         --readFilesIn "$fq" \
         --readFilesCommand zcat \
         --outFileNamePrefix "$BAM_DIR/${base}_" \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts \
         >"$LOG_DIR/${base}_star.log" 2>&1
         
    if [ -s "$bam_file" ]; then
        samtools index "$bam_file"
        echo "$base alignment + indexing done."
    else
        echo "⚠️ Alignment failed for $base."
    fi
done
echo "=== STAR alignment completed ==="
```
7. Read Counting (featureCounts)

Gene-level quantification was performed using the GFF3 annotation.
```bash
# ---------------------------
# 6. featureCounts
# ---------------------------
echo "=== Counting reads with featureCounts ==="
featureCounts -O -t gene -g ID \
  -a "$GTF_DIR/TAIR10.gff3" \
  -o "$COUNT_DIR/counts.txt" "$BAM_DIR"/*.bam \
  >"$LOG_DIR/featureCounts.log" 2>&1

echo "=== featureCounts done. Counts file at $COUNT_DIR/counts.txt ==="
```
