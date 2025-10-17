#!/bin/bash

# qc.sh - Run FastQC and MultiQC on raw FASTQ files

RAW_DIR="raw_fastq"
QC_DIR="qc_results"
LOG_DIR="logs"

mkdir -p "$QC_DIR" "$LOG_DIR"

echo "=== Running FastQC on all FASTQs in $RAW_DIR ==="
for fq in "$RAW_DIR"/*.fastq.gz; do
    sample=$(basename "$fq")
    echo ">>> Processing $sample"
    fastqc "$fq" \
        --outdir "$QC_DIR" \
        >"$LOG_DIR/${sample%.fastq.gz}_fastqc.log" 2>&1
done

echo "=== Summarizing with MultiQC ==="
multiqc "$QC_DIR" \
    --outdir "$QC_DIR" \
    >"$LOG_DIR/multiqc.log" 2>&1

echo "=== QC completed! Reports are in $QC_DIR ==="
