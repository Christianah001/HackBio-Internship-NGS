#!/bin/bash
# Script: run_quast.sh
# Purpose: Run QUAST on all assemblies

# ==============================
# Input and output directories
# ==============================
ASSEMBLY_DIR="./selected/assembly"
QUAST_DIR="./selected/quast_results"

# Create output directory if it doesn’t exist
mkdir -p $QUAST_DIR

echo "=== STARTING QUAST ASSEMBLY QUALITY CHECK ==="
echo "Input assemblies: $ASSEMBLY_DIR"
echo "Output directory: $QUAST_DIR"
echo ""

# ==============================
# Run QUAST
# ==============================
# Find all contigs.fasta from assemblies and run quast
quast.py -o $QUAST_DIR -t 12 $ASSEMBLY_DIR/*/contigs.fasta

echo ""
echo "=== QUAST COMPLETED SUCCESSFULLY ==="
echo "Results saved in: $QUAST_DIR"
echo "Key file: $QUAST_DIR/report.html (open in browser)"
