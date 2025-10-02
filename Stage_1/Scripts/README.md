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
```
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
Randomly choose 50 samples
Bash Script 2: `Selected.samples.sh` (Random Selection)

```
#!/bin/bash
# 2. Selected.samples.sh
# Set how many samples you want
N=50
# Randomly pick N sample prefixes
shuf -n $N all_samples.txt > selected_samples.txt
# Inspect selected
cat selected_samples.txt | sed -n '1,20p'
```



