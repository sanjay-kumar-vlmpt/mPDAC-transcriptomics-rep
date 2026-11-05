# ============================================================
# 00 - Install packages (run ONCE per computer, not every session)
# Project: Transcriptomic deconstruction of mPDAC
# ============================================================

install.packages("BiocManager")

BiocManager::install(c(
  "TCGAbiolinks", "SummarizedExperiment", "DESeq2",
  "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "enrichplot"
))

install.packages(c(
  "ggplot2", "ggrepel", "pheatmap", "stringr",
  "survival", "survminer"
))
