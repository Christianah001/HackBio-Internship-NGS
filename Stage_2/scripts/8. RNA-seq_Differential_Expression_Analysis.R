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
