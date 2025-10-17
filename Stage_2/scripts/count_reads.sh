#!/usr/bin/bash
# count_reads.sh
# Description: Count reads from STAR-aligned BAM files using featureCounts

# === Directories ===
BAM_DIR="bam"
COUNT_DIR="counts"
ANNOTATION="reference/annotation/TAIR10.gff3"

mkdir -p "$COUNT_DIR"

echo "=== Starting featureCounts ==="
echo "Using annotation file: $ANNOTATION"
echo "BAM files from: $BAM_DIR"

# List BAM files being counted
echo "BAM files to be counted:"
ls -lh "$BAM_DIR"/*.bam
echo "----------------------------"

# === Run featureCounts ===
     echo "📊 Running featureCounts with GFF3 annotation..."
  featureCounts -O -t gene -g ID \
  -a reference/annotation/TAIR10.gff3 \
  -o counts/counts.txt bam/*.bam
  
    echo "✅ featureCounts completed!"
    echo "Counts saved in: $COUNT_DIR/counts.txt"
  # Quick summary of counts
    if [ -f "$COUNT_DIR/counts.txt" ]; then
        echo "Summary of counts per sample:"
        cut -f1,7- "$COUNT_DIR/counts.txt" | head -n 10
        echo "----------------------------"
    
else
    echo "❌ featureCounts failed. Check BAM file permissions and annotation path."
fi

echo "All done!"
