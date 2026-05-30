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

# Download data 
GDCdownload(query, directory = "data/raw/TCGA-PAAD")

# Prepare data into a usable format
tcga_data <- GDCprepare(query, directory = "data/raw/TCGA-PAAD")

# Save 
saveRDS(tcga_data, file = "data/processed/tcga_paad_raw.rds")

# Number of samples and genes
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
# Save 
save.session <- getwd()
cat("Working directory:", save.session)
setwd("C:/Sanjay/mPDAC Project")
tcga_data <- readRDS("data/processed/tcga_paad_raw.rds")
setwd("C:/Sanjay/mPDAC Project")
tcga_data <- readRDS("data/processed/tcga_paad_raw.rds")
# Extract clinical information
clinical_data <- as.data.frame(colData(tcga_data))

# Number of Patients
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
# Extract unstranded raw counts
count_matrix <- assay(tcga_data, "unstranded")

# Check dimensions
dim(count_matrix)
count_matrix[1:5, 1:3]
# Save count matrix
saveRDS(count_matrix, file = "data/processed/tcga_paad_counts.rds")

# Save clinical data
saveRDS(clinical_data, file = "data/processed/tcga_paad_clinical.rds")
