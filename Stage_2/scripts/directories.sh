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
