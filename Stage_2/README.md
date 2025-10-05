# Differential Gene Expression Analysis of Arabidopsis Vasculature Under UV-C Stress

**Author:**  
Funmilayo C. Ligali  
---

## Introduction

Berkowitz et al. (2021) studied tissue-specific responses of *Arabidopsis thaliana* leaves to abiotic stresses. In this project, we focus on the **vasculature tissue** and determine genes responding to UV-C treatment compared to water-treated control samples. RNA-seq data from three biological replicates per condition were analyzed.

UV-C is a high-energy radiation that induces DNA damage and stress response pathways. Understanding the transcriptional response in vasculature helps uncover stress adaptation mechanisms.

---

## Objectives

- Download and process raw RNA-seq data (FASTQ files) for control and UV-C treated samples.  
- Perform quality control, trimming, and alignment of reads to the *Arabidopsis thaliana* genome (TAIR10).  
- Quantify gene expression using **featureCounts**.  
- Perform differential expression analysis between UV-C treated and control samples using **DESeq2**.  
- Identify top 100 differentially expressed genes.  
- Conduct **GO** and **KEGG** enrichment analysis to determine significantly affected biological pathways.

---

## Materials & Data

| Replicate | Control SRR   | UV-C Treated SRR |
|-----------|---------------|----------------|
| 1         | SRR12808527   | SRR12808497    |
| 2         | SRR12808528   | SRR12808498    |
| 3         | SRR12808529   | SRR12808499    |

**Reference genome:** *Arabidopsis thaliana* TAIR10 (FASTA & GFF3 from Ensembl Plants)

## Methods
## RNA-seq Analysis Pipeline

| Step | Tool / Software | Input | Output | Description |
|------|----------------|-------|--------|-------------|
| 1 | Bash | SRR IDs | Raw FASTQ files | Download raw RNA-seq data from ENA |
| 2 | FastQC / MultiQC | Raw FASTQ | QC reports | Assess quality of raw reads |
| 3 | fastp | Raw FASTQ | Trimmed FASTQ, HTML/JSON report | Adapter and quality trimming |
| 4 | FastQC / MultiQC | Trimmed FASTQ | QC reports | Assess quality post-trimming |
| 5 | STAR | Trimmed FASTQ, TAIR10 genome | Sorted BAM files | Align reads to Arabidopsis genome |
| 6 | Samtools | BAM files | Indexed BAM | Index aligned BAM for downstream use |
| 7 | featureCounts | BAM, GFF3 annotation | Counts table | Gene-level read quantification |
| 8 | R / DESeq2 | Counts table, metadata | DESeq2 results | Differential expression analysis |
| 9 | pheatmap (R) | DESeq2 results | Heatmaps | Visualize differentially expressed genes |
| 10 | clusterProfiler | Top DE genes | GO enrichment | Identify significantly enriched pathways |

### Bash Scripts
# RNA-seq Preprocessing and Quantification Pipeline

**Author:** Funmilayo C. Ligali  
**Date:** 2025-10-04  

This pipeline automates **RNA-seq preprocessing** from raw FASTQ downloads to gene-level quantification, specifically for Arabidopsis vasculature UV-C stress study.

---

## Pipeline Overview

This pipeline performs the following steps:

1. Create project directories.  
2. Download raw FASTQ files from ENA using SRR IDs.  
3. Quality control with `FastQC` and `MultiQC`.  
4. Trim adapters and low-quality bases using `fastp`, followed by QC.  
5. Download the reference genome (TAIR10) and annotation (GFF3) and generate STAR genome index.  
6. Align trimmed reads to the genome with `STAR` and index BAM files with `samtools`.  
7. Quantify gene-level counts using `featureCounts`.

---

## Bash Pipeline
# Complete RNA-seq preprocessing and quantification pipeline

0. All project files and outputs were organized into directories
```bash
#!/bin/bash
# rna_seq_pipeline.sh

# ---------------------------
# 0. Set Directories
# ---------------------------
BASE_DIR=$(pwd)
RAW_DIR="$BASE_DIR/raw_fastq"
TRIM_DIR="$BASE_DIR/trimmed_fastq"
QC_DIR="$BASE_DIR/qc_results"
QC_TRIM="$QC_DIR/trimmed"
REF_DIR="$BASE_DIR/reference"
GENOME_DIR="$REF_DIR/genome"
GTF_DIR="$REF_DIR/annotation"
STAR_INDEX_DIR="$REF_DIR/STAR_index"
BAM_DIR="$BASE_DIR/bam"
COUNT_DIR="$BASE_DIR/counts"
LOG_DIR="$BASE_DIR/logs"

mkdir -p "$RAW_DIR" "$TRIM_DIR" "$QC_DIR" "$QC_TRIM" "$GENOME_DIR" \
         "$GTF_DIR" "$STAR_INDEX_DIR" "$BAM_DIR" "$COUNT_DIR" "$LOG_DIR"
```
1. Download Raw RNA-seq FASTQ Files

SRR IDs for control and UV-C treated samples were downloaded from ENA using FTP links.
```bash
# ---------------------------
# 1. Download SRR FASTQs
# ---------------------------
cat > "$BASE_DIR/srrs.txt" <<'EOF'
SRR12808527
SRR12808528
SRR12808529
SRR12808497
SRR12808498
SRR12808499
EOF

echo "=== Downloading FASTQ files ==="
while read SRR; do
    echo "Processing $SRR..."
    curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$SRR&result=read_run&fields=run_accession,fastq_ftp,library_layout" \
    | tail -n +2 | while IFS=$'\t' read -r run_accession fastq_ftp library_layout; do
        for url in $(echo $fastq_ftp | tr ';' ' '); do
            fname=$(basename $url)
            echo "Downloading $fname ($library_layout)"
            curl -L "ftp://$url" -o "$RAW_DIR/$fname"
        done
    done
done < "$BASE_DIR/srrs.txt"
echo "=== FASTQ download completed ==="
```
2. Quality Control (FastQC + MultiQC)

Raw FASTQ files were assessed for quality before trimming.
```bash
# ---------------------------
# 2. QC raw FASTQs
# ---------------------------
echo "=== Running FastQC on raw FASTQs ==="
for fq in "$RAW_DIR"/*.fastq.gz; do
    sample=$(basename "$fq")
    fastqc "$fq" --outdir "$QC_DIR" >"$LOG_DIR/${sample%.fastq.gz}_fastqc.log" 2>&1
done

echo "=== MultiQC summary ==="
multiqc "$QC_DIR" --outdir "$QC_DIR" >"$LOG_DIR/multiqc.log" 2>&1
```
3. Trimming (fastp) and QC

Reads were trimmed to remove adapters and low-quality bases, then re-checked with FastQC.
```bash
# ---------------------------
# 3. Trim with fastp + QC
# ---------------------------
echo "=== Trimming raw FASTQs ==="
for fq in "$RAW_DIR"/*.fastq.gz; do
    base=$(basename "$fq" .fastq.gz)
    fastp -i "$fq" -o "$TRIM_DIR/${base}_trim.fastq.gz" \
          -h "$TRIM_DIR/${base}_fastp.html" \
          -j "$TRIM_DIR/${base}_fastp.json" \
          >>"$LOG_DIR/fastp.log" 2>&1
done

echo "=== FastQC on trimmed FASTQs ==="
for fq in "$TRIM_DIR"/*.fastq.gz; do
    base=$(basename "$fq" .fastq.gz)
    fastqc "$fq" --outdir "$QC_TRIM" >>"$LOG_DIR/fastqc_trimmed.log" 2>&1
done

multiqc "$QC_TRIM" --outdir "$QC_TRIM" >>"$LOG_DIR/multiqc_trimmed.log" 2>&1
echo "=== Trimming + QC completed ==="
```
4. Reference Genome Download and STAR Index

The Arabidopsis TAIR10 genome and GFF3 annotation were downloaded from Ensembl Plants, then indexed for STAR alignment.
```bash
# ---------------------------
# 4. Download reference genome + STAR index
# ---------------------------
echo "=== Downloading Arabidopsis genome & annotation ==="
wget -c https://ftp.ensemblgenomes.org/pub/plants/release-62/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz \
     -O "$GENOME_DIR/TAIR10.fa.gz"
wget -c https://ftp.ensemblgenomes.org/pub/plants/release-62/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.62.gff3.gz \
     -O "$GTF_DIR/TAIR10.gff3.gz"
gunzip -f "$GENOME_DIR/TAIR10.fa.gz"
gunzip -f "$GTF_DIR/TAIR10.gff3.gz"

echo "=== Generating STAR genome index ==="
STAR --runThreadN 12 \
     --runMode genomeGenerate \
     --genomeDir "$STAR_INDEX_DIR" \
     --genomeFastaFiles "$GENOME_DIR/TAIR10.fa" \
     --sjdbGTFfile "$GTF_DIR/TAIR10.gff3" \
     --sjdbOverhang 99
echo "=== STAR genome index ready ==="
```
5. STAR Alignment

Trimmed reads were aligned to TAIR10 with STAR. BAM files were sorted and indexed.
```bash
# ---------------------------
# 5. STAR Alignment
# ---------------------------
echo "=== Aligning trimmed FASTQs with STAR ==="
for fq in "$TRIM_DIR"/*_trim.fastq.gz; do
    base=$(basename "$fq" _trim.fastq.gz)
    bam_file="$BAM_DIR/${base}_Aligned.sortedByCoord.out.bam"
    [ -s "$bam_file" ] && { echo "$base BAM exists, skipping."; continue; }
    
    STAR --runThreadN 12 \
         --genomeDir "$STAR_INDEX_DIR" \
         --readFilesIn "$fq" \
         --readFilesCommand zcat \
         --outFileNamePrefix "$BAM_DIR/${base}_" \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts \
         >"$LOG_DIR/${base}_star.log" 2>&1
         
    if [ -s "$bam_file" ]; then
        samtools index "$bam_file"
        echo "$base alignment + indexing done."
    else
        echo "⚠️ Alignment failed for $base."
    fi
done
echo "=== STAR alignment completed ==="
```
7. Read Counting (featureCounts)

Gene-level quantification was performed using the GFF3 annotation.
```bash
# ---------------------------
# 6. featureCounts
# ---------------------------
echo "=== Counting reads with featureCounts ==="
featureCounts -O -t gene -g ID \
  -a "$GTF_DIR/TAIR10.gff3" \
  -o "$COUNT_DIR/counts.txt" "$BAM_DIR"/*.bam \
  >"$LOG_DIR/featureCounts.log" 2>&1

echo "=== featureCounts done. Counts file at $COUNT_DIR/counts.txt ==="
```

# 8. RNA-seq Differential Expression Analysis

```r
# ---------------------------------------------------------------
# 0. Set Working Directory
# ---------------------------------------------------------------
# Set the directory where the project files (counts.txt, Metadata.csv) are stored
setwd("C:/Users/ELITEBOOK 1040 G3/Desktop/HackBio NGS/Stage 2")

# ---------------------------------------------------------------
# 1. Install and Load Required Packages
# ---------------------------------------------------------------
# Bioconductor manager for package installation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# List of required packages for DE analysis and visualization
packages <- c("DESeq2", "pheatmap", "dplyr", "ggplot2",
              "clusterProfiler", "org.At.tair.db")

# Loop through packages: install if missing, then load
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
  library(pkg, character.only = TRUE)
}

# Install additional package for log fold change shrinkage
install.packages("ashr")

# ---------------------------------------------------------------
# 2. Import Data
# ---------------------------------------------------------------
# Load featureCounts output as count matrix
counts <- read.delim("counts.txt")

# Load metadata describing sample conditions
metadata <- read.csv("Metadata.csv")

# Quick view of the first rows
head(counts)
head(metadata)

# ---------------------------------------------------------------
# 3. Clean and Align Metadata
# ---------------------------------------------------------------
# Assign descriptive sample names
metadata$sample <- c("Control_1", "Control_2", "Control_3",
                     "Treated_1", "Treated_2", "Treated_3")

# Set condition as a factor with control as reference
metadata$condition <- factor(metadata$condition, levels = c("control", "treated"))

# Preview metadata
head(metadata)

# ---------------------------------------------------------------
# 4. Prepare Count Matrix
# ---------------------------------------------------------------
# Align column names in count matrix with metadata sample names
colnames(counts)[7:12] <- metadata$sample

# Extract count data for analysis
raw_counts <- counts[, c("Geneid", metadata$sample)]

# Set rownames as Gene IDs
rownames(raw_counts) <- raw_counts$Geneid

# Remove the Geneid column from the data frame
raw_counts <- raw_counts[, -1]

# Verify column names match metadata
all(colnames(raw_counts) == metadata$sample)

# Preview the prepared count matrix
head(raw_counts)

# ---------------------------------------------------------------
# 5. Create DESeq2 Dataset
# ---------------------------------------------------------------
# Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData   = metadata,
                              design    = ~ condition)

# Relevel factor to set control as reference for DE analysis
dds$condition <- relevel(dds$condition, ref = "control")

# ---------------------------------------------------------------
# 6. Run DESeq2 Analysis
# ---------------------------------------------------------------
# Run differential expression analysis
dds <- DESeq(dds)

# Extract DE results (log2 fold change and p-values)
res <- results(dds)

# Shrink log2 fold changes for more accurate estimates
res <- lfcShrink(dds, coef = "condition_treated_vs_control", type = "ashr")

# Preview DE results
head(res)
summary(res)

# ---------------------------------------------------------------
# 7. Volcano Plot (interactive)
# ---------------------------------------------------------------
# Convert results to data frame
res_df <- as.data.frame(res) %>% na.omit()
res_df$GeneID <- rownames(res_df)

# Define significance thresholds
padj_cutoff <- 0.05
lfc_cutoff  <- 1

# Volcano Plot
res_df <- res_df %>%
  mutate(
    DE = case_when(
      log2FoldChange > lfc_cutoff & padj < padj_cutoff ~ "Up Regulated",
      log2FoldChange < -lfc_cutoff & padj < padj_cutoff ~ "Down Regulated",
      TRUE ~ "Not Significant"
    )
  )

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = DE)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("Up Regulated" = "salmon", "Down Regulated" = "lightblue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Volcano Plot: Vasculature UV-C vs Control",
       x = expression(log[2]~"Fold Change"),
       y = expression(-log[10]~"Adjusted P-value"))

# Heatmap for Top 50 DE Genes
ordered_cols <- c("Control_1","Control_2","Control_3", "Treated_1","Treated_2","Treated_3")
sample_group <- data.frame(Treatment = factor(c("Control","Control","Control","UV-C","UV-C","UV-C"), levels=c("Control", "UV-C")))
rownames(sample_group) <- ordered_cols

top_50 <- res_df %>% arrange(padj) %>% head(50)
top_50_genes <- top_50$GeneID

# Extract Normalized/Regularized Log Transformed counts for visualization
vsd <- vst(dds, blind=FALSE)
vsd_mat <- assay(vsd)[top_50_genes, ordered_cols]

# Plot heatmap
pheatmap(vsd_mat,
         scale = "row",
         cluster_cols = FALSE,
         show_rownames = FALSE, 
         annotation_col = sample_group,
         main = "Heatmap of Top 50 DE Genes (VSD Normalized)")
# ---------------------------------------------------------------
# 8. List Top 100 Differentially Expressed Genes
# ---------------------------------------------------------------
# Filter for significant DE genes (padj < 0.05 and |log2FoldChange| > 1)
sig_res_df <- res_df %>%
  filter(padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff)

# Order significant genes by adjusted p-value (smallest first)
ordered_sig_res_df <- sig_res_df[order(sig_res_df$padj), ]

# List the top 100 differentially expressed genes
top_100_de_genes <- head(ordered_sig_res_df, 100)

print("Top 10 Differentially Expressed Genes:")
print(head(top_100_de_genes))

write.csv(top_100_de_genes, "top_100_DE_genes_vasculature_UV.csv", row.names = FALSE)

# ---------------------------------------------------------------
# 9 Functional Enrichment Analysis (GO and KEGG)
significant_genes <- top_100_de_genes$GeneID

# Print a few cleaned IDs to verify the format
print(head(cleaned_significant_genes))

# Get all valid key types for the annotation package
keytypes(org.At.tair.db)

# Get a list of actual TAIR Gene IDs present in the database
head(keys(org.At.tair.db, keytype="TAIR"))

# REMOVE VERSION NUMBER
# If the GeneIDs contain version numbers (e.g., AT1G01010.1), remove them
cleaned_significant_genes_v1 <- sub("\\..*$", "", significant_genes)

# REMOVE "gene:" PREFIX (THIS IS THE CRUCIAL STEP)
# Inspection shows IDs are "gene:AT...", which is not a valid TAIR key.
final_cleaned_genes <- sub("^gene:", "", cleaned_significant_genes_v1)

# Map the final cleaned IDs to ENTREZID
gene_list_mapped <- bitr(final_cleaned_genes,
                         fromType = "TAIR",
                         toType = "ENTREZID",
                         OrgDb = org.At.tair.db,
                         drop = TRUE)

# ---------------------------------------------------------------
# 10. Functional Enrichment Analysis (KEGG and GO)
# ---------------------------------------------------------------

# The list of final cleaned TAIR Locus IDs:
final_cleaned_genes <- sub("^gene:", "", sub("\\..*$", "", significant_genes))

# --- A. KEGG Pathway Enrichment (Use TAIR IDs Directly for 'ath') ---
print("--- Running KEGG Enrichment (Using TAIR Locus IDs) ---")
kegg_results <- enrichKEGG(gene         = final_cleaned_genes, # USE TAIR LOCUS IDs HERE
                           organism     = 'ath', 
                           pvalueCutoff = 0.05)

# --- B. GO Term Enrichment (Use TAIR IDs directly) ---
print("--- Running GO Enrichment (Biological Process) ---")
go_results <- enrichGO(gene          = final_cleaned_genes, # Use the same cleaned TAIR IDs
                       OrgDb         = org.At.tair.db,
                       keyType       = 'TAIR',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.1)

# Simplify GO results
go_results_simplified <- simplify(go_results)


# --- C. Extract Top 5 Enriched Pathways ---
safe_to_df <- function(res, source) {
  if (is.null(res) || nrow(res@result) == 0) {
    return(data.frame(ID=character(), Description=character(), p.adjust=numeric(), Pathway_Source=character(), stringsAsFactors=FALSE))
  }
  as.data.frame(res@result) %>%
    mutate(Pathway_Source = source) %>%
    dplyr::select(ID, Description, p.adjust, Pathway_Source)
}


if (!is.null(go_results_simplified) && nrow(go_results_simplified@result) > 0) {
  print("Generating GO Dot Plot...")
  dotplot(go_results_simplified, 
          showCategory=5, 
          title="GO Biological Process Enrichment (Top 5)")
} else {
  print("Skipping GO Dot Plot: No significant GO terms found.")
}


# Extract and combine results
go_df <- safe_to_df(go_results_simplified, "GO_BP")
kegg_df <- safe_to_df(kegg_results, "KEGG")

all_enriched_pathways <- bind_rows(go_df, kegg_df) %>%
  arrange(p.adjust)


# List of top 5 enriched pathways
top_5_pathways <- head(all_enriched_pathways, 5)

print("Top 5 Enriched Pathways:")
print(top_5_pathways)

if (!is.null(kegg_results) && nrow(kegg_results@result) > 0) {
  print("Generating KEGG Dot Plot...")
  dotplot(kegg_results, 
          showCategory=5, 
          title="KEGG Pathway Enrichment (Top 5)")
} else {
  print("Skipping KEGG Dot Plot: No significant KEGG terms found.")
}

# ---------------------------------------------------------------
# 11. Optional Visualization
# ---------------------------------------------------------------

# Heatmap (using VSD normalized counts)
top_50_genes <- top_100_de_genes$GeneID[1:50]
ordered_cols <- rownames(metadata)
vsd <- vst(dds, blind=FALSE)
vsd_mat <- assay(vsd)[top_50_genes, ordered_cols]
sample_group <- data.frame(Treatment = metadata$condition)
rownames(sample_group) <- ordered_cols

pheatmap(vsd_mat, scale = "row", cluster_cols = FALSE, show_rownames = FALSE, 
         annotation_col = sample_group, main = "Heatmap of Top 50 DE Genes")

```
### Result
```
```
Volcano Plot
![Volcano Plot](https://github.com/Christianah001/HackBio-Internship-NGS/blob/main/Stage_2/Results/Volcano%20res_df.png)
```
```
The plot shows a clear and significant shift in gene expression. The vast majority of genes significantly affected (those scattered above the red dashed line, which represents your padj<0.05 threshold) show either strong up-regulation (red dots to the right of the log 2 FC>1 threshold) or strong down-regulation (blue dots to the left of the log 
2 FC<−1 threshold).

The UV-C treatment caused a large-scale, highly significant transcriptional reprogramming in the vascular tissue. The most significant changes are clustered well away from the center, confirming a genuine biological effect.
```
```
Heatmap of Top 100 Differentially Expressed Genes
![Heatmap of Differentially expressed genes](https://github.com/Christianah001/HackBio-Internship-NGS/blob/main/Stage_2/Results/DE%20Genes.png)
```
```
The Heatmap provides irrefutable visual validation of the experiment's success.
The clustering tree at the top of the plot clearly separates the samples into two distinct branches: Control (cyan bar) and UV-C Treated (pink/red bar). This shows the biological variation between the two conditions is much greater than the variation within the three biological replicates.

Gene Expression Patterns:

The left block (Control samples) shows high expression (red/orange) for one group of genes and low expression (blue) for another.
The right block (Treated samples) shows the inverse pattern: the genes highly expressed in control are now low, and the genes low in control are now highly induced (red/orange).

Key Takeaway: The clear separation validates the entire experiment. The UV-C treatment effectively and consistently switched the transcriptional state of the vascular cells.
It displays the expression values of the most significant genes, where red signifies high expression and blue signifies low expression. The dendrogram (tree) structure at the top shows a perfect, clear segregation of all samples into two major clusters: Control and Treated. This confirms that the UV-C exposure caused a massive and consistent transcriptional shift. Furthermore, the genes cluster into two primary groups: a large set of genes that are highly downregulated (switched off) in the UV-C samples (likely growth/metabolism genes), and a smaller, highly upregulated set that represents the core defense and repair machinery activated by the stress.
```
```
Functional Interpretation (GO and KEGG Enrichment)
![GO](https://github.com/Christianah001/HackBio-Internship-NGS/blob/main/Stage_2/Results/GO%20top5.png)
```
```
Interpretation:
The GO plot reveals the urgent biological actions the Arabidopsis vasculature immediately initiates to cope with the stress:

Damage Control: The most significant terms, "response to oxidative stress" and "response to hypoxia," are directly related to damage control. UV-C light generates harmful molecules (ROS), and the plant is urgently activating defenses to neutralize this toxicity and mitigate the resulting internal cellular distress.

Physical Repair: The term "response to wounding" is also highly enriched. While UV-C is not a physical wound, the cellular damage it causes is so severe that the plant triggers generic damage repair mechanisms to maintain tissue integrity.

In short, the GO analysis confirms the vasculature is in an immediate crisis management state, focusing all available resources on detoxification and cellular repair.
```
```
![KEGG](https://github.com/Christianah001/HackBio-Internship-NGS/blob/main/Stage_2/Results/kegg%20top%205.png)
```
```
Interpretation:

Bubble Size (Count): Indicates how many DE genes belong to that specific pathway. Larger bubbles mean more biological coverage.

Bubble Color (p.adjust): Indicates statistical significance. Bubbles that are red (closer to 0.01) are more highly significant than those that are blue (closer to 0.03).

X-axis (GeneRatio): The proportion of DE genes found in the pathway relative to the total number of genes in that pathway.

Key Biological Findings:
he KEGG plot shows the plant's immediate plan for survival, which is dominated by three crucial actions:

Systemic Signaling (MAPK Pathway): This pathway is the most significant. It acts as the plant's central communication hub, immediately sensing the UV threat and sending "emergency alert" signals across the leaf to coordinate the entire defense response.

Detoxification (Glutathione Metabolism): Upregulation of this pathway is a direct fight against the damage. Glutathione is a potent antioxidant that is rapidly deployed to neutralize the toxic Reactive Oxygen Species (ROS) generated by the UV-C light, protecting the cells from being poisoned.

Defense Repurposing (Pathogen Interaction): The plant quickly turns on generic stress defenses, reusing pathways it normally uses to fight bacteria or fungi to combat the UV-C physical stress, maximizing its ability to survive the attack.

In short, the vasculature acts as the plant's command center, initiating signaling, detoxification, and a strong general defense.

