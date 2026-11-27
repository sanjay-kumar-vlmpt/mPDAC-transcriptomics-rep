# ============================================================
# 10 - Portable basal-classical score and metastasis vs primary (GSE71729)
#
# Score = mean z-score of basal signature genes - mean z-score of classical genes
# (Moffitt et al. 2015 signature; each gene's direction confirmed in TCGA and GSE71729).
# Works in datasets without subtype labels. + = basal-like, - = classical-like.
#
# Caveats: sample titles do not identify patients (metastases may share patients);
# metastases likely represent later-stage disease than resected primaries; bulk tissue.
#
# Output: data/processed/subtype_score_genes.rds
#         results/tables/GSE71729_subtype_scores.csv
#         results/figures/GSE71729_score_by_organ.png, GSE71729_score_primary_vs_met.png
# ============================================================

setwd("C:/Sanjay/mPDAC Project")
while (!is.null(dev.list())) dev.off()

library(Biobase)

eset    <- readRDS("data/processed/gse71729_eset.rds")
samples <- readRDS("data/processed/gse71729_samples.rds")
expr    <- exprs(eset)
tcga_bc <- read.csv("results/tables/deseq2_basal_vs_classical.csv",         row.names = 1)
geo_bc  <- read.csv("results/tables/GSE71729_limma_basal_vs_classical.csv", row.names = 1)
stopifnot(all(colnames(expr) == rownames(samples)))

# 1. Moffitt signature; direction of each gene decided by the data
moffitt_genes <- c("MYO1A","TNS4","CYP3A7","SPRR3","TFF2","SCEL","KRT17","CST6","FGFBP1",
  "VGLL1","LYZ","AGR2","CEACAM6","PRR15L","LY6D","KRT6C","SERPINB4","ANXA8L1","UCA1",
  "PLA2G10","SPRR1B","ST6GALNAC1","S100A2","KRT15","CTSV","KRT6A","VSIG2","TSPAN8",
  "LGALS4","FAM83A","REG4","LEMD1","AGR3","TFF3","SLC2A1","TFF1","KRT7","CDH17","DHRS9",
  "KRT20","SPINK4","CTSE","ANXA10","AREG","CLRN3","BTNL8","SERPINB3","FAM3D","GPR87")
sig <- Reduce(intersect, list(moffitt_genes, rownames(tcga_bc), rownames(geo_bc)))
dir_tcga <- sign(tcga_bc[sig, "log2FoldChange"])     # +1 = higher in basal
dir_geo  <- sign(geo_bc[sig, "logFC"])
table(TCGA = dir_tcga, GEO = dir_geo)                # 25 basal / 24 classical, full agreement

agree <- dir_tcga == dir_geo
basal_genes     <- sig[agree & dir_tcga > 0]
classical_genes <- sig[agree & dir_tcga < 0]

# GEO matrix uses ~2013 names; genes renamed since are dropped (ANXA8L1, CTSV and one classical gene)
basal_use <- basal_genes[basal_genes %in% rownames(expr)]
class_use <- classical_genes[classical_genes %in% rownames(expr)]
setdiff(c(basal_genes, classical_genes), c(basal_use, class_use))
length(basal_use); length(class_use)                 # 23 / 23
saveRDS(list(basal = basal_use, classical = class_use), "data/processed/subtype_score_genes.rds")

# 2. Score every non-cell-line sample
keep <- !is.na(samples$tissue)
z <- t(scale(t(expr[c(basal_use, class_use), keep])))
s <- samples[keep, ]
s$score <- colMeans(z[basal_use, , drop = FALSE]) - colMeans(z[class_use, , drop = FALSE])
s$title <- pData(eset)[rownames(s), "title"]
write.csv(s[, c("title", "source", "tissue", "subtype", "score")],
          "results/tables/GSE71729_subtype_scores.csv")

# 3. Sanity check: score separates labelled subtypes in primaries
p   <- s$tissue == "Primary"
lab <- ifelse(s$subtype[p] == "2", "Basal", "Classical")
tapply(s$score[p], lab, median)                      # 0.81 vs -0.42
wilcox.test(s$score[p] ~ lab)                        # p < 2.2e-16 (expected: labels from this dataset)

# 4. Organ background: score in normal tissues
nrm <- s$tissue == "Normal"
round(sort(tapply(s$score[nrm], s$source[nrm], median)), 2)

png("results/figures/GSE71729_score_by_organ.png", width = 10, height = 6, units = "in", res = 300)
par(mar = c(10, 4, 3, 1))
boxplot(score ~ source, data = s[s$tissue %in% c("Primary", "Normal"), ],
        las = 2, cex.axis = 0.8, main = "Basal-classical score: primaries and normal organs",
        ylab = "Score (+ = basal-like)", xlab = "")
abline(h = 0, lty = 2)
dev.off()

# 5. Metastases vs primaries
pm <- s[s$tissue %in% c("Primary", "Metastasis"), ]
tapply(pm$score, pm$tissue, median)                  # met 0.10 vs primary -0.19
wilcox.test(score ~ tissue, data = pm)               # p = 0.015

met <- s[s$tissue == "Metastasis", ]; met$organ <- sub("_Metastasis", "", met$source)
nrm_df <- s[nrm, ];                   nrm_df$organ <- sub("_Normal", "", nrm_df$source)
data.frame(
  n_met        = as.vector(table(met$organ)),
  met_median   = round(tapply(met$score, met$organ, median), 2),
  organ_normal = round(tapply(nrm_df$score, nrm_df$organ, median)[names(table(met$organ))], 2)
)

liver <- s[s$source %in% c("Pancreas_Primary", "Liver_Metastasis"), ]
tapply(liver$score, liver$source, median)            # 0.28 vs -0.19 (liver background -0.12)
wilcox.test(score ~ source, data = liver)            # p = 0.010

grp <- ifelse(s$tissue == "Primary", "Primary",
       ifelse(s$source == "Liver_Metastasis", "Liver met",
       ifelse(s$tissue == "Metastasis", "Other met",
       ifelse(s$source == "Liver_Normal", "Liver normal",
       ifelse(s$source == "Pancreas_Normal", "Pancreas normal", NA)))))
grp <- factor(grp, levels = c("Pancreas normal", "Primary", "Liver normal", "Liver met", "Other met"))
png("results/figures/GSE71729_score_primary_vs_met.png", width = 9, height = 6, units = "in", res = 300)
boxplot(s$score ~ grp, col = c("grey80", "#4DBBD5", "grey80", "#E64B35", "#F39B7F"),
        ylab = "Basal-classical score (+ = basal-like)", xlab = "",
        main = "GSE71729: primaries vs metastases")
abline(h = 0, lty = 2)
dev.off()
