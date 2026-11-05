# ============================================================
# 02 - Sample filtering, low-count gene filtering, gene annotation
# Input : data/processed/tcga_paad_counts.rds, tcga_paad_clinical.rds
# Output: data/processed/tcga_paad_counts_annotated.rds
#         data/processed/tcga_paad_clinical_filtered.rds
# ============================================================

setwd("C:/Sanjay/mPDAC Project")

library(org.Hs.eg.db)
library(AnnotationDbi)

count_matrix  <- readRDS("data/processed/tcga_paad_counts.rds")
clinical_data <- readRDS("data/processed/tcga_paad_clinical.rds")

# 1. Keep primary tumours and solid tissue normals (drops the 1 metastatic sample)
keep_samples <- clinical_data$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")
count_matrix_filtered <- count_matrix[, keep_samples]
clinical_filtered     <- clinical_data[keep_samples, ]

# 2. Keep genes with at least 10 counts in at least 10 samples
keep_genes <- rowSums(count_matrix_filtered >= 10) >= 10
count_matrix_filtered <- count_matrix_filtered[keep_genes, ]

# 3. Ensembl IDs -> HGNC symbols (strip version suffix first; drop unmapped genes)
gene_ids_clean <- gsub("\\..*", "", rownames(count_matrix_filtered))
gene_symbols <- mapIds(org.Hs.eg.db, keys = gene_ids_clean, column = "SYMBOL",
                       keytype = "ENSEMBL", multiVals = "first")
valid_genes <- !is.na(gene_symbols)
count_matrix_annotated <- count_matrix_filtered[valid_genes, ]
rownames(count_matrix_annotated) <- gene_symbols[valid_genes]
# Note: a few symbols map from >1 Ensembl ID; DESeq2 later renames these duplicates.

dim(count_matrix_annotated)   # 20,829 genes x 182 samples

saveRDS(count_matrix_annotated, "data/processed/tcga_paad_counts_annotated.rds")
saveRDS(clinical_filtered,      "data/processed/tcga_paad_clinical_filtered.rds")
