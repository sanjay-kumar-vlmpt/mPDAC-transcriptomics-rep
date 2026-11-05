# Transcriptomic Deconstruction of Metastatic Pancreatic Ductal Adenocarcinoma (mPDAC)

An ongoing computational project using public bulk and single-cell RNA-seq data to identify candidate biomarkers for metastatic PDAC. Full methods, results and interpretation are in the [project wiki](../../wiki).

## Current status

**Phase 1 (TCGA-PAAD, bulk RNA-seq): primary PDAC vs normal pancreas.** This stage is *not* metastasis-specific. TCGA-PAAD has only one metastatic sample. It gives a baseline of tumour-associated genes for the later metastasis analyses.

- Cohort: **150 confirmed PDAC tumours** (the TCGA 2017 paper's curated set) **vs 4 normal tissues**. An earlier version used all 178 tumours labelled "Primary Tumor", which included non-PDAC tumours (see the wiki's correction log).
- **1,090 high-confidence DEGs** (padj < 0.05, |log2FC| > 1): 664 up, 426 down.
- Upregulated: gut-like (digestion, gastrointestinal epithelium) and epidermal programs, plus cell division. Downregulated: mainly immune-associated genes, which in bulk data cannot be separated from lower immune-cell content.
- SI and CT83 have extreme fold changes because they are undetected in the 4 normals (on/off, not effect size).
- Survival (150 tumours, 82 deaths): PSCA shows a consistent but modest trend (median split p = 0.06; age-adjusted Cox HR 1.10 per log2 unit, p = 0.024). SI (p = 0.8) and KLRF1 (p = 0.4) show no association.

**Limitations:** only 4 normal samples; bulk tissue (purity and composition confound immune genes); no independent validation yet.

## Next steps

1. Validate in independent bulk cohorts (GSE71729, GSE21501)
2. Primary vs metastatic comparison at single-cell level (GSE154778, Seurat)
3. Pan-cancer and normal-tissue specificity (TCGA Pan-Cancer, GTEx)

## Repository structure

```
scripts/
  00_install_packages.R      one-time package installation
  01_download_tcga.R         TCGA-PAAD download (TCGAbiolinks, STAR - Counts)
  02_preprocess.R            sample + low-count filtering, Ensembl -> gene symbol
  03_deseq2_confirmed150.R   cohort definition, DESeq2, volcano, heatmap
  04_enrichment.R            GO BP and KEGG over-representation
  05_subtype_purity.R        Moffitt subtype and tumour-purity checks
  06_survival.R              Kaplan-Meier, log-rank and Cox models
results/
  tables/                    DESeq2 results, DEG lists, enrichment tables
  figures/                   volcano plot, heatmap
legacy/
  mPDAC_original_178tumours.R   first-pass analysis, superseded (kept for transparency)
```

Raw and processed data (`data/`) are not tracked. Run scripts 01 and 02 to regenerate them. Set the `setwd()` path at the top of each script to your own project folder.

## Environment

R 4.5 on Windows, RStudio. Key packages: TCGAbiolinks, DESeq2, clusterProfiler, org.Hs.eg.db, enrichplot, survival, survminer, ggplot2, pheatmap.

## Data source

TCGA-PAAD via the NCI Genomic Data Commons. Curated sample annotations from The Cancer Genome Atlas Research Network, *Cancer Cell* 2017; 32(2):185–203.
