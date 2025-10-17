#!/bin/bash
# align_star_fixed.sh
# Robust STAR alignment for Arabidopsis trimmed RNA-seq reads

# === Directories ===
TRIM_DIR="trimmed_fastq"
GENOME_DIR="reference/STAR_index"
BAM_DIR="bam"
LOG_DIR="logs/star_logs"

# === Create necessary directories ===
mkdir -p "$BAM_DIR" "$LOG_DIR"

echo "=== STAR Alignment Started ==="

# Loop over all trimmed FASTQ files
for fq in "$TRIM_DIR"/*_trim.fastq.gz; do
    [ -e "$fq" ] || continue  # skip if no files found
    base=$(basename "$fq" _trim.fastq.gz)
    echo ">>> Aligning $base ..."

    # Check if BAM already exists
    bam_file="$BAM_DIR/${base}_Aligned.sortedByCoord.out.bam"
    if [ -s "$bam_file" ]; then
        echo "✅ $base BAM already exists, skipping alignment."
        continue
    fi

    # Run STAR
    STAR --runThreadN 12 \
         --genomeDir "$GENOME_DIR" \
         --readFilesIn "$fq" \
         --readFilesCommand zcat \
         --outFileNamePrefix "$BAM_DIR/${base}_" \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts \
         >"$LOG_DIR/${base}_star.log" 2>&1

    # Check if BAM is generated and non-empty
    if [ -s "$bam_file" ]; then
        echo "✅ $base alignment completed!"
        # Index BAM for downstream use
        samtools index "$bam_file"
        echo "🔹 $base BAM indexed."
    else
        echo "⚠️ Alignment failed for $base. Check logs: $LOG_DIR/${base}_star.log"
    fi
done

echo "=== STAR Alignment Finished ==="
