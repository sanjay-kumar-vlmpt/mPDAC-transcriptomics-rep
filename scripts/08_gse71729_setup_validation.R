# ============================================================
# 08 - GSE71729 (Moffitt 2015, Agilent two-colour microarray):
#      download, sample table, and validation of the TCGA tumour-vs-normal genes
#
# Data notes (from the authors' processing description):
#   - values are log2 Cy5 (sample channel only; the reference RNA is not used)
#   - negative signals were set to 0, so 0 means "below detection"
#   - probes already collapsed to gene symbols; samples quantile normalised
#   - gene symbols are ~2013 names, so they are updated before comparing with TCGA
#
# Output: data/processed/gse71729_eset.rds, gse71729_samples.rds
#         results/tables/GSE71729_limma_tumor_vs_normal.csv
# ============================================================

setwd("C:/Sanjay/mPDAC Project")
dir.create("data/raw/GEO", recursive = TRUE, showWarnings = FALSE)

library(GEOquery)
library(limma)
library(org.Hs.eg.db)
options(timeout = 600)

# 1. Download the series matrix
gse  <- getGEO("GSE71729", GSEMatrix = TRUE, destdir = "data/raw/GEO")
eset <- gse[[1]]
dim(exprs(eset))      # 19,749 genes x 357 samples
annotation(eset)      # GPL20769
saveRDS(eset, "data/processed/gse71729_eset.rds")

# 2. Clean sample table from the parsed ":ch2" columns (channel 2 = the real sample)
pd <- pData(eset)
pd$data_processing[1]
samples <- data.frame(
  id      = pd$geo_accession,
  source  = pd$source_name_ch2,
  tissue  = pd$`tissue type:ch2`,                          # NA = cell lines
  subtype = pd$`tumor_subtype_0na_1classical_2basal:ch2`,  # 0 = NA, 1 = classical, 2 = basal
  surv_mo = as.numeric(pd$`survival_months:ch2`),
  death   = as.numeric(pd$`death_event_1death_0censor:ch2`),
  row.names = pd$geo_accession
)
table(samples$tissue, useNA = "ifany")                   # 61 metastasis, 134 normal, 145 primary, 17 NA
table(samples$tissue, samples$subtype, useNA = "ifany")
saveRDS(samples, "data/processed/gse71729_samples.rds")

# 3. limma: 145 primary tumours vs 46 normal pancreas
expr <- exprs(eset)
stopifnot(all(colnames(expr) == rownames(samples)))
keep_s <- samples$source %in% c("Pancreas_Primary", "Pancreas_Normal")
x   <- expr[, keep_s]
grp <- factor(ifelse(samples$source[keep_s] == "Pancreas_Primary", "Tumor", "Normal"),
              levels = c("Normal", "Tumor"))
keep_g <- rowSums(x > 1) >= 23                    # above 1 in at least half the normal group
fit <- eBayes(lmFit(x[keep_g, ], model.matrix(~ grp)))
tt  <- topTable(fit, coef = "grpTumor", number = Inf)
summary(decideTests(fit, p.value = 0.05, lfc = 1))   # 285 up / 175 down

# 4. Update old gene symbols so they match TCGA
tt$symbol <- alias2SymbolTable(rownames(tt), species = "Hs")
tt <- tt[!is.na(tt$symbol) & !duplicated(tt$symbol), ]
rownames(tt) <- tt$symbol
write.csv(tt, "results/tables/GSE71729_limma_tumor_vs_normal.csv")

# 5. Do the TCGA high-confidence genes replicate?
up_tcga   <- rownames(read.csv("results/tables/upregulated_DEGs_confirmed150.csv",   row.names = 1))
down_tcga <- rownames(read.csv("results/tables/downregulated_DEGs_confirmed150.csv", row.names = 1))
up_in   <- intersect(up_tcga,   rownames(tt))
down_in <- intersect(down_tcga, rownames(tt))
length(up_in); length(down_in)                                               # 531 / 324 measured
sum(tt[up_in,   "adj.P.Val"] < 0.05 & tt[up_in,   "logFC"] > 0)              # 343 replicate (65%)
sum(tt[down_in, "adj.P.Val"] < 0.05 & tt[down_in, "logFC"] < 0)              # 71 replicate (22%)

tcga <- read.csv("results/tables/deseq2_tumor_vs_normal_confirmed150.csv", row.names = 1)
both <- intersect(c(up_in, down_in), rownames(tcga))
cor.test(tcga[both, "log2FoldChange"], tt[both, "logFC"], method = "spearman")   # rho 0.40

tt[intersect(c("PSCA", "SI", "CT83", "KLRF1", "AGR2", "TFF2", "KRT5"), rownames(tt)),
   c("logFC", "adj.P.Val")]

# 6. Is basal more common among metastases than primaries? (labels)
tab <- table(samples$tissue, samples$subtype)[c("Primary", "Metastasis"), c("1", "2")]
tab
fisher.test(tab)   # OR 2.7 (1.4-5.3), p = 0.002
