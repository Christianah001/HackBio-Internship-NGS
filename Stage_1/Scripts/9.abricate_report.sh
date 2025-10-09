#!/bin/bash
# Script: abricate_master_pipeline.sh
# Description: Run Abricate (AMR + toxin genes), summarize results, and generate prevalence tables + report.
# Author: Funmilayo Ligali

# -----------------------------
# Setup
# -----------------------------
ASSEMBLY_DIR="./selected/assembly"
ABRICATE_DIR="./results/abricate"
REPORT_DIR="./results/report"

mkdir -p "$ABRICATE_DIR/amr" "$ABRICATE_DIR/toxin" "$ABRICATE_DIR/summary" "$REPORT_DIR"

echo "=== ABRicate Master Pipeline ==="
echo "Input assemblies: $ASSEMBLY_DIR"
echo "Output: $ABRICATE_DIR"

success_count=0
total_count=0

# -----------------------------
# Run Abricate for each sample
# -----------------------------
for assembly_dir in "$ASSEMBLY_DIR"/*; do
    sample_name=$(basename "$assembly_dir")
    contigs_file="$assembly_dir/contigs.fasta"

    total_count=$((total_count + 1))

    if [[ -f "$contigs_file" && -s "$contigs_file" ]]; then
        echo "Processing sample: $sample_name"

        # AMR (CARD)
        abricate --db card --quiet "$contigs_file" > "$ABRICATE_DIR/amr/${sample_name}_amr.tsv"

        # Toxin genes (VFDB)
        abricate --db vfdb --quiet "$contigs_file" > "$ABRICATE_DIR/toxin/${sample_name}_toxin.tsv"

        success_count=$((success_count + 1))
        echo "✓ ABRicate completed for $sample_name"
    else
        echo "✗ Skipping $sample_name (no contigs found)"
    fi
done

# -----------------------------
# Summarize results
# -----------------------------
echo ""
echo "Generating summary tables..."

# Summaries
abricate --summary "$ABRICATE_DIR/amr"/*.tsv > "$ABRICATE_DIR/summary/amr_summary.tsv"
abricate --summary "$ABRICATE_DIR/toxin"/*.tsv > "$ABRICATE_DIR/summary/toxin_summary.tsv"

# Combine all results
cat "$ABRICATE_DIR/amr"/*.tsv > "$ABRICATE_DIR/summary/all_amr_results.tsv"
cat "$ABRICATE_DIR/toxin"/*.tsv > "$ABRICATE_DIR/summary/all_toxin_results.tsv"

# -----------------------------
# Generate prevalence table (CSV)
# -----------------------------
echo "Sample,Gene" > "$ABRICATE_DIR/summary/amr_gene_prevalence.csv"
for f in "$ABRICATE_DIR/amr"/*.tsv; do
    sample=$(basename "$f" | sed 's/_amr.tsv//')
    awk 'NR>1 {print "'"$sample"'," $2}' "$f" >> "$ABRICATE_DIR/summary/amr_gene_prevalence.csv"
done

echo "Sample,Gene" > "$ABRICATE_DIR/summary/toxin_gene_prevalence.csv"
for f in "$ABRICATE_DIR/toxin"/*.tsv; do
    sample=$(basename "$f" | sed 's/_toxin.tsv//')
    awk 'NR>1 {print "'"$sample"'," $2}' "$f" >> "$ABRICATE_DIR/summary/toxin_gene_prevalence.csv"
done

# -----------------------------
# Human-readable report
# -----------------------------
REPORT_FILE="$REPORT_DIR/AMR_Toxin_Report.txt"
echo "AMR and Toxin Gene Report" > "$REPORT_FILE"
echo "=========================" >> "$REPORT_FILE"
echo "" >> "$REPORT_FILE"

echo "Total assemblies checked: $total_count" >> "$REPORT_FILE"
echo "Successful ABRicate analyses: $success_count" >> "$REPORT_FILE"
echo "" >> "$REPORT_FILE"

echo "== AMR Genes (CARD Summary) ==" >> "$REPORT_FILE"
awk 'NR>1 {print $1, $2, $3, $5}' "$ABRICATE_DIR/summary/amr_summary.tsv" >> "$REPORT_FILE"
echo "" >> "$REPORT_FILE"

echo "== Suggested Antibiotic Notes ==" >> "$REPORT_FILE"
echo "- bla genes: avoid beta-lactams (ampicillin, penicillin)" >> "$REPORT_FILE"
echo "- aac/aph genes: avoid aminoglycosides (gentamicin, kanamycin)" >> "$REPORT_FILE"
echo "- erm genes: avoid macrolides (erythromycin, azithromycin)" >> "$REPORT_FILE"
echo "- van genes: avoid vancomycin" >> "$REPORT_FILE"
echo "" >> "$REPORT_FILE"

echo "== Toxin/Virulence Genes (VFDB Summary) ==" >> "$REPORT_FILE"
awk 'NR>1 {print $1, $2, $3, $5}' "$ABRICATE_DIR/summary/toxin_summary.tsv" >> "$REPORT_FILE"

echo ""
echo "=== Pipeline Completed ==="
echo "Summary files created in: $ABRICATE_DIR/summary"
echo "Human-readable report: $REPORT_FILE"
echo "Prevalence tables: "
echo "  - $ABRICATE_DIR/summary/amr_gene_prevalence.csv"
echo "  - $ABRICATE_DIR/summary/toxin_gene_prevalence.csv"
