#!/bin/bash
# trim_qc.sh - Trim raw FASTQs and run QC
# Usage: bash trim_qc.sh

# Directories
RAW_DIR="raw_fastq"
TRIM_DIR="trimmed_fastq"
QC_TRIM="qc_results/trimmed"
LOG_DIR="logs"

# Create directories if they don't exist
mkdir -p "$TRIM_DIR" "$QC_TRIM" "$LOG_DIR"

echo "=== Trimming with fastp ==="

for fq in "$RAW_DIR"/*.fastq.gz; do
    base=$(basename "$fq" .fastq.gz)
    echo "Trimming $base ..."
    fastp \
        -i "$fq" \
        -o "$TRIM_DIR/${base}_trim.fastq.gz" \
        -h "$TRIM_DIR/${base}_fastp.html" \
        -j "$TRIM_DIR/${base}_fastp.json" \
        >>"$LOG_DIR/fastp.log" 2>&1
done

echo "=== Running FastQC on trimmed FASTQs ==="
fastqc "$TRIM_DIR"/*.fastq.gz \
echo "=== Running FastQC on trimmed FASTQs ==="
for fq in "$TRIM_DIR"/*.fastq.gz; do
    base=$(basename "$fq" .fastq.gz)
    echo "FastQC on $base ..."
    fastqc "$fq" \
        --outdir "$QC_TRIM" \
        >>"$LOG_DIR/fastqc_trimmed.log" 2>&1
done

echo "=== Summarizing with MultiQC (trimmed) ==="
multiqc "$QC_TRIM" \
    --outdir "$QC_TRIM" \
    >>"$LOG_DIR/multiqc_trimmed.log" 2>&1

echo "=== Trimming + QC completed! Reports are in $QC_TRIM ==="
