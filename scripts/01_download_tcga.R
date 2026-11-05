# ============================================================
# 01 - Download TCGA-PAAD RNA-seq counts (STAR - Counts, GRCh38)
# Output: data/processed/tcga_paad_raw.rds
#         data/processed/tcga_paad_counts.rds
#         data/processed/tcga_paad_clinical.rds
# ============================================================

setwd("C:/Sanjay/mPDAC Project")   # change to your own project folder
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

library(TCGAbiolinks)
library(SummarizedExperiment)

query <- GDCquery(
  project       = "TCGA-PAAD",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query, directory = "data/raw/TCGA-PAAD")
tcga_data <- GDCprepare(query, directory = "data/raw/TCGA-PAAD")
saveRDS(tcga_data, "data/processed/tcga_paad_raw.rds")

# Raw (unstranded) counts and sample-level clinical data
count_matrix  <- assay(tcga_data, "unstranded")
clinical_data <- as.data.frame(colData(tcga_data))

dim(count_matrix)                 # genes x samples
table(clinical_data$sample_type)  # 178 primary tumour, 4 normal, 1 metastatic

saveRDS(count_matrix,  "data/processed/tcga_paad_counts.rds")
saveRDS(clinical_data, "data/processed/tcga_paad_clinical.rds")
