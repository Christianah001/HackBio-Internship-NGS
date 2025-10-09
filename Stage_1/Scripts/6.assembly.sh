#!/bin/env bash
#6. assembly.sh
# Description: Assemble trimmed paired-end reads into genomes using SPAdes v3.13.1

# Set directories
TRIMMED_DATA_DIR="./selected/trimmed"
ASSEMBLY_DIR="./selected/assembly"

# Create output directory
echo "Creating assembly output directory: $ASSEMBLY_DIR"
mkdir -p "$ASSEMBLY_DIR"

# Check if trimmed data exists
if [ -z "$(ls -A $TRIMMED_DATA_DIR/*_1.trim.fastq.gz 2>/dev/null)" ]; then
    echo "Error: No trimmed FASTQ files found in $TRIMMED_DATA_DIR!"
    exit 1
fi

echo "=== GENOME ASSEMBLY WITH SPADES ==="
echo "Input: $TRIMMED_DATA_DIR"
echo "Output: $ASSEMBLY_DIR"

# Initialize counters
sample_count=0
success_count=0

# Process each paired-end sample
for r1 in "$TRIMMED_DATA_DIR"/*_1.trim.fastq.gz; do
    base_name=$(basename "$r1" _1.trim.fastq.gz)
    r2="${TRIMMED_DATA_DIR}/${base_name}_2.trim.fastq.gz"

    sample_count=$((sample_count + 1))
    echo ""
    echo "Processing sample $sample_count: $base_name"

    # Create output directory for this sample
    sample_outdir="${ASSEMBLY_DIR}/${base_name}"
    mkdir -p "$sample_outdir"

    # Run SPAdes assembly
    echo "Running SPAdes assembly for $base_name..."
    spades.py \
        -1 "$r1" \
        -2 "$r2" \
        -o "$sample_outdir" \
        --careful \
        --cov-cutoff auto \
        -t 4 \
        -m 32 \
        --phred-offset 33

    # Check if assembly was successful
    if [ -f "${sample_outdir}/contigs.fasta" ] && [ -s "${sample_outdir}/contigs.fasta" ]; then
        echo "✓ Assembly successful: ${sample_outdir}/contigs.fasta"
        success_count=$((success_count + 1))
        # Optional: basic statistics
        echo "Assembly stats for $base_name:"
        echo "Number of contigs: $(grep -c '^>' ${sample_outdir}/contigs.fasta)"
        echo "Total length: $(awk '/^>/ {if (seqlen) print seqlen; seqlen=0; next} {seqlen += length($0)} END {print seqlen}' ${sample_outdir}/contigs.fasta) bp"
    else
        echo "✗ Assembly failed for $base_name"
    fi
done

echo ""
echo "=== ASSEMBLY SUMMARY ==="
echo "Total samples processed: $sample_count"
echo "Successful assemblies: $success_count"
echo "Assemblies saved in: $ASSEMBLY_DIR"
echo ""
echo "Next step: Run your downstream QC or annotation scripts."
