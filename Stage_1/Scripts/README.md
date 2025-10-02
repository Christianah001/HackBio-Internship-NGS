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
| I. Sample Preparation | Selection, Organization | shuf, ln -s (Bash) | Randomly selects the 50 isolates and uses symbolic links to organize the raw sequencing data (FASTQ files) for processing. |
| II. QC & Trimming | Quality Assessment, Filtering | FastQC, MultiQC, fastp | Assesses read quality, removes sequencing adaptors, and trims low-quality bases to ensure clean input for assembly.|
| III. Assembly & QC | Genome Construction, Assessment | SPAdes, QUAST | Assembles cleaned reads into draft contiguous genome sequences (contigs.fasta). QUAST evaluates assembly quality metrics (N50, contig count). |
| IV. Gene Analysis | Organism ID, Resistance/Virulence Screening | BLASTn, ABRicate, Pandas | Confirms species identity. ABRicate screens assemblies against the CARD (AMR) and VFDB (Toxin) databases. Python calculates prevalence and generates visualizations. |

# 3. Functional Scripts
The following scripts automate the WGS pipeline used for this analysis.

## 3.1. Phase 1: Sample Preparation Scripts
Bash Script 1: '''Selected.samples.sh (Random Selection)
