## Comparative Genomic Analysis of Listeria monocytogenes Isolates from the South African Polony Outbreak (2017–2018)
# 1. Introduction and Objectives
This report details the bioinformatics analysis of 50 Listeria monocytogenes isolates, a subset of the strains responsible for the devastating 2017–2018 South African listeriosis outbreak. The project utilizes Whole-Genome Sequencing (WGS) data to genetically characterize the pathogen, which is crucial for informing public health response and clinical guidance.

## The objectives are:
a. Confirm the Organism Identity

b. Determine the Antimicrobial Resistance (AMR) Profile

c. Detect Virulence (Toxin) Genes

d. Suggest Evidence-Based Treatment Options

# 2. Methodology: WGS Pipeline Overview
The analysis follows a comprehensive microbial WGS pipeline, structured into four main phases using a suite of bioinformatics tools and custom Bash/Python scripts.

| Phase | Steps | Tools Used | Purpose |
| --- | --- | --- | --- |
| I. Sample Preparation | Selection, Organization | `shuf`, `ln -s` (Bash) | Randomly selects the 50 isolates and uses symbolic links to organize the raw sequencing data (FASTQ files) for processing. |
| II. QC & Trimming | Quality Assessment, Filtering | FastQC, MultiQC, fastp | Assesses read quality, removes sequencing adaptors, and trims low-quality bases to ensure clean input for assembly.|
| III. Assembly & QC | Genome Construction, Assessment | SPAdes, QUAST | Assembles cleaned reads into draft contiguous genome sequences (contigs.fasta). QUAST evaluates assembly quality metrics (N50, contig count). |
| IV. Gene Analysis | Organism ID, Resistance/Virulence Screening | BLASTn, ABRicate, Pandas | Confirms species identity. ABRicate screens assemblies against the CARD (AMR) and VFDB (Toxin) databases. Python calculates prevalence and generates visualizations. |

# 3. Functional Scripts
The following scripts automate the WGS pipeline used for this analysis.
## 3.1. Phase 1: Sample Preparation Scripts
### Identify sample name pattern and list sample prefixes
Bash Script 1: `link_samples.sh` 
```bash
#!/bin/bash
#1. show example files (first 20)
ls -1 | head -n 20
# If files are paired and named SAMPLE_1.fastq.gz and SAMPLE_2.fastq.gz:
ls *_1*.fastq.gz | head
# Create list of prefixes (remove the trailing _1.fastq.gz)
for f in *_1*.fastq.gz; do echo "${f%%_1*}"; done | sort -u > all_samples.txt
# Count how many unique prefixes found
wc -l all_samples.txt
```
### Randomly choose 50 samples
Bash Script 2: `Selected.samples.sh` (Random Selection)

```bash
#!/bin/bash
# 2. Selected.samples.sh
# Set how many samples you want
N=50
# Randomly pick N sample prefixes
shuf -n $N all_samples.txt > selected_samples.txt
# Inspect selected
cat selected_samples.txt | sed -n '1,20p'
```

### Post-Trimming QC
Bash Script 5: `trim_qc.sh`

```bash
#!/usr/bin/bash
#4. trim_qc.sh
# Description: Quality control on trimmed data

TRIMMED_DATA_DIR="/home/maa/Funmilayo/WGS_microbes/selected/trimmed"
QC_DIR="$TRIMMED_DATA_DIR/trimmed_qc"

# Create QC output directory
mkdir -p "$QC_DIR"
echo "Creating QC directory: $QC_DIR"

# Check if there are trimmed FASTQ files
if [ -z "$(ls -A $TRIMMED_DATA_DIR/*_*.trim.fastq.gz 2>/dev/null)" ]; then
    echo "Error: No trimmed FASTQ files found in $TRIMMED_DATA_DIR"
    exit 1
fi

echo "=== Running FastQC on trimmed data ==="
fastqc -o "$QC_DIR" "$TRIMMED_DATA_DIR"/*_*.trim.fastq.gz

echo "=== Running MultiQC on trimmed QC reports ==="
multiqc "$QC_DIR" -o "$QC_DIR"

echo "QC reports saved to: $QC_DIR"
```

## 3.3. Phase 3: Assembly & QC Scripts
Bash Script 6: `assembly.sh` (SPAdes Assembly)

```bash
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
        echo "Total length: $(awk '/^>/ {if (seqlen) print seqlen; seqlen=0; next} {seqlen += length($0)} END {print seqlen}' ${sample_outdir}/contigs.fasta) b>
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
```

