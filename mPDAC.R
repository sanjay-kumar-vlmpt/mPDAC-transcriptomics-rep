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

#Installing Clusterprofiler package
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
#Loading Library
library(clusterProfiler)
library(org.Hs.eg.db)

setwd("C:/Sanjay/mPDAC Project")

#Uploading Tables of UR and DR DEGs
upregulated <- read.csv("results/tables/upregulated_DEGs.csv", row.names = 1)
downregulated <- read.csv("results/tables/downregulated_DEGs.csv", row.names = 1)

# Combining upregulated and downregulated gene symbols
all_sig_genes <- c(rownames(upregulated), rownames(downregulated))

# Converting gene symbols to Entrez IDs for clusterProfiler
gene_entrez <- bitr(all_sig_genes, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)

head(gene_entrez)
nrow(gene_entrez)

#Identifying Biological Processes using Gene Ontology 
go_bp <- enrichGO(gene = gene_entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

head(go_bp@result, 10)

library(enrichplot)

# Dot plot of top 15 enriched terms
dotplot(go_bp, showCategory = 15) + 
  ggtitle("GO Biological Process Enrichment - PDAC DEGs")


library(ggplot2)

# Fewer categories + larger plot + text wrapping
dotplot(go_bp, showCategory = 10) + 
  ggtitle("GO Biological Process Enrichment - PDAC DEGs") +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40))

ggsave("results/figures/GO_BP_dotplot.png", 
       width = 11, height = 9, dpi = 300)

# Convert upregulated genes
up_entrez <- bitr(rownames(upregulated), 
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = org.Hs.eg.db)

# Convert downregulated genes
down_entrez <- bitr(rownames(downregulated), 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)

# Check counts
nrow(up_entrez)
nrow(down_entrez)

# GO enrichment for upregulated genes
go_bp_up <- enrichGO(gene = up_entrez$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     readable = TRUE)

# GO enrichment for downregulated genes
go_bp_down <- enrichGO(gene = down_entrez$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)

#Dot plot for Upregulated genes
dotplot(go_bp_up, showCategory = 10) + 
  ggtitle("GO Biological Process - Upregulated Genes (Tumour vs Normal)") +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40))

ggsave("results/figures/GO_BP_upregulated.png", 
       width = 11, height = 9, dpi = 300)

#Dot plot for Downregulated genes
dotplot(go_bp_down, showCategory = 10) + 
  ggtitle("GO Biological Process - Downregulated Genes (Tumour vs Normal)") +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40))

ggsave("results/figures/GO_BP_downregulated.png", 
       width = 11, height = 9, dpi = 300)

# KEGG enrichment for upregulated genes
kegg_up <- enrichKEGG(gene = up_entrez$ENTREZID,
                      organism = "hsa",
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH")

head(kegg_up@result[, c("Description", "GeneRatio", "p.adjust", "Count")], 10)

# KEGG enrichment for downregulated genes
kegg_down <- enrichKEGG(gene = down_entrez$ENTREZID,
                        organism = "hsa",
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH")

head(kegg_down@result[, c("Description", "GeneRatio", "p.adjust", "Count")], 10)

#Dot plot for KEGG_DOWN
dotplot(kegg_down, showCategory = 10) + 
  ggtitle("KEGG Pathways - Downregulated Genes (Tumour vs Normal)") +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40))

ggsave("results/figures/KEGG_downregulated.png", 
       width = 11, height = 9, dpi = 300)

#Dot plot for KEGG_UP
dotplot(kegg_up, showCategory = 10) + 
  ggtitle("KEGG Pathways - Upregulated Genes (Tumour vs Normal)") +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40))

ggsave("results/figures/KEGG_upregulated.png", 
       width = 11, height = 9, dpi = 300)

library(survival)
library(survminer)
setwd("C:/Sanjay/mPDAC Project")
clinical_filtered <- readRDS("data/processed/tcga_paad_clinical_filtered.rds")
count_matrix_annotated <- readRDS("data/processed/tcga_paad_counts_annotated.rds")

# Survival time and event status
clinical_filtered$time <- ifelse(
  clinical_filtered$vital_status == "Dead",
  clinical_filtered$days_to_death,
  clinical_filtered$days_to_last_follow_up
)

clinical_filtered$event <- ifelse(
  clinical_filtered$vital_status == "Dead", 1, 0
)

# Checking for missing values
sum(is.na(clinical_filtered$time))
table(clinical_filtered$event)

# Filtering tumor samples only for survival analysis
tumor_samples <- clinical_filtered$sample_type == "Primary Tumor"
clinical_tumor <- clinical_filtered[tumor_samples, ]
count_matrix_tumor <- count_matrix_annotated[, tumor_samples]

dim(clinical_tumor)
dim(count_matrix_tumor)
table(clinical_tumor$event)

# Extraction of PSCA expression across tumor samples
# Since PSCA is well established in literature as a marker for PDAC
psca_expr <- count_matrix_tumor["PSCA", ]

# Median splitting for High and Low
psca_group <- ifelse(psca_expr > median(psca_expr), "High", "Low")

# Including to clinical data
clinical_tumor$psca_group <- psca_group

clinical_tumor$psca_group <- psca_group

# Checking the median split
table(clinical_tumor$psca_group)

# Creating survival object
surv_object <- Surv(time = clinical_tumor$time, event = clinical_tumor$event)

# Fitting of Kaplan-Meier curves by PSCA group
fit_psca <- survfit(surv_object ~ psca_group, data = clinical_tumor)

print(fit_psca)

# Log-rank test
survdiff(surv_object ~ psca_group, data = clinical_tumor)

ggsurvplot(fit_psca, 
           data = clinical_tumor,
           pval = TRUE,
           risk.table = TRUE,
           conf.int = TRUE,
           xlab = "Time (days)",
           ylab = "Survival Probability",
           title = "Kaplan-Meier Survival Curve - PSCA Expression in PDAC",
           legend.title = "PSCA Expression",
           legend.labs = c("High", "Low"),
           palette = c("#E64B35", "#4DBBD5"))

png("results/figures/KM_PSCA_survival.png", 
    width = 9, height = 8, units = "in", res = 300)

print(ggsurvplot(fit_psca, 
                 data = clinical_tumor,
                 pval = TRUE,
                 risk.table = TRUE,
                 conf.int = TRUE,
                 xlab = "Time (days)",
                 ylab = "Survival Probability",
                 title = "Kaplan-Meier Survival Curve - PSCA Expression in PDAC",
                 legend.title = "PSCA Expression",
                 legend.labs = c("High", "Low"),
                 palette = c("#E64B35", "#4DBBD5")))

dev.off()