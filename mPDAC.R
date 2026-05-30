# Install BiocManager first
install.packages("BiocManager")

# Core analysis packages
BiocManager::install(c(
  "DESeq2",
  "edgeR",
  "limma",
  "TCGAbiolinks",
  "GEOquery",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot",
  "Seurat"
))

# Visualisation & utilities
install.packages(c(
  "ggplot2",
  "pheatmap",
  "ggrepel",
  "dplyr",
  "tidyr",
  "readr",
  "rmarkdown",
  "devtools"
))

# Survival analysis
install.packages(c(
  "survival",
  "survminer"
))

# ============================================
# Script 01 - Download TCGA-PAAD Dataset
# Project: mPDAC Transcriptomic Analysis
# ============================================

# Load the library
library(TCGAbiolinks)

# Query TCGA for pancreatic adenocarcinoma RNA-seq data
query <- GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Download the data (this will take a few minutes)
GDCdownload(query, directory = "data/raw/TCGA-PAAD")

# Prepare the data into a usable format
tcga_data <- GDCprepare(query, directory = "data/raw/TCGA-PAAD")

# Save it so you don't have to re-download
saveRDS(tcga_data, file = "data/processed/tcga_paad_raw.rds")

# Quick check - how many samples and genes do we have?
dim(tcga_data)

dir.exists("data/raw")
library(TCGAbiolinks)
BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
query <- GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query, directory = "data/raw/TCGA-PAAD")
tcga_data <- GDCprepare(query, directory = "data/raw/TCGA-PAAD")
saveRDS(tcga_data, file = "data/processed/tcga_paad_raw.rds")
# Save your session
save.session <- getwd()
cat("Working directory:", save.session)
setwd("C:/Sanjay/mPDAC Project")
tcga_data <- readRDS("data/processed/tcga_paad_raw.rds")
setwd("C:/Sanjay/mPDAC Project")
tcga_data <- readRDS("data/processed/tcga_paad_raw.rds")
# Extract clinical information
clinical_data <- as.data.frame(colData(tcga_data))

# How many patients and variables do we have?
dim(clinical_data)
library(TCGAbiolinks)
clinical_data <- as.data.frame(colData(tcga_data))
library(SummarizedExperiment)
clinical_data <- as.data.frame(colData(tcga_data))
dim(clinical_data)
clinical_data <- as.data.frame(colData(tcga_data))
dim(clinical_data)
colnames(clinical_data)
# Look at vital status and tumour stage
table(clinical_data$vital_status)
table(clinical_data$ajcc_pathologic_stage)
# Extract the unstranded raw counts
count_matrix <- assay(tcga_data, "unstranded")

# Check dimensions
dim(count_matrix)
count_matrix[1:5, 1:3]
# Save count matrix
saveRDS(count_matrix, file = "data/processed/tcga_paad_counts.rds")

# Save clinical data
saveRDS(clinical_data, file = "data/processed/tcga_paad_clinical.rds")

setwd("C:/Sanjay/mPDAC Project")
BiocManager::install("DESeq2")
library(DESeq2)

results_df <- read.csv("results/tables/deseq2_tumor_vs_normal.csv", 
                       row.names = 1)

# Filter significant genes only
sig_genes <- results_df[!is.na(results_df$padj) & 
                          results_df$padj < 0.05, ]

# Number of significant genes
nrow(sig_genes)

write.csv(sig_genes, file = "results/tables/significant_DEGs.csv")

# Filter by both significance AND fold change
sig_genes_strict <- sig_genes[abs(sig_genes$log2FoldChange) > 1, ]

# Split upregulated from downregulated
upregulated <- sig_genes_strict[sig_genes_strict$log2FoldChange > 1, ]
downregulated <- sig_genes_strict[sig_genes_strict$log2FoldChange < -1, ]

# How many in each direction
nrow(upregulated)
nrow(downregulated)

#Save File for UR and DR
write.csv(upregulated, file = "results/tables/upregulated_DEGs.csv")
write.csv(downregulated, file = "results/tables/downregulated_DEGs.csv")


library(ggplot2)
library(ggrepel)

# Preparing data for volcano plot
volcano_data <- results_df
volcano_data$gene <- rownames(results_df)
volcano_data$significance <- "Not Significant"
volcano_data$significance[volcano_data$padj < 0.05 & 
                            volcano_data$log2FoldChange > 1] <- "Upregulated"
volcano_data$significance[volcano_data$padj < 0.05 & 
                            volcano_data$log2FoldChange < -1] <- "Downregulated"

# Top 10 genes to label
top_genes <- c(rownames(upregulated)[1:5], 
               rownames(downregulated)[1:5])

# Create volcano plot
ggplot(volcano_data, aes(x = log2FoldChange, 
                         y = -log10(padj), 
                         color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "#E64B35",
                                "Downregulated" = "#4DBBD5",
                                "Not Significant" = "grey60")) +
  geom_text_repel(data = subset(volcano_data, gene %in% top_genes),
                  aes(label = gene), size = 3, max.overlaps = 20) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title = "Differential Gene Expression in PDAC",
       subtitle = "Tumour vs Normal — TCGA-PAAD",
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-value",
       color = "Expression") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
ggsave("results/figures/volcano_plot_PDAC.png", 
       width = 10, height = 8, dpi = 300)

library(pheatmap)

# Top 25 upregulated and downregulated genes
top25_up <- rownames(upregulated)[1:25]
top25_down <- rownames(downregulated)[1:25]
top50_genes <- c(top25_up, top25_down)

# Extracting counts from the annotated matrix
count_matrix_annotated <- readRDS("data/processed/tcga_paad_counts_annotated.rds")
heatmap_matrix <- count_matrix_annotated[top50_genes, ]

# Check dimensions
dim(heatmap_matrix)

# Log transform counts (adding 1 to avoid log of zero)
heatmap_matrix_log <- log2(heatmap_matrix + 1)

# Scale by row (gene) so each gene's expression is relative to its own mean
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix_log)))

clinical_filtered <- readRDS("data/processed/tcga_paad_clinical_filtered.rds")

# Creating annotation for samples - tumor vs normal
sample_annotation <- data.frame(
  Condition = clinical_filtered$sample_type,
  row.names = colnames(heatmap_matrix)
)

# Plot heatmap
pheatmap(heatmap_matrix_scaled,
         annotation_col = sample_annotation,
         show_colnames = FALSE,
         fontsize_row = 7,
         color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
         main = "Top 50 DEGs — PDAC Tumour vs Normal",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         border_color = NA)
png("results/figures/heatmap_top50_DEGs.png", 
    width = 10, height = 12, units = "in", res = 300)

pheatmap(heatmap_matrix_scaled,
         annotation_col = sample_annotation,
         show_colnames = FALSE,
         fontsize_row = 7,
         color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
         main = "Top 50 DEGs — PDAC Tumour vs Normal",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         border_color = NA)

dev.off()