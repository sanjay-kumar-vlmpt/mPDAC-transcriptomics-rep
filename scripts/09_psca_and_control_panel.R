# ============================================================
# 09 - PSCA follow-up and the positive-control gene panel
#   (a) PSCA survival in GSE71729 (independent cohort)
#   (b) Is PSCA subtype-linked? (TCGA and GSE71729)
#   (c) GSE71729 basal vs classical (limma, primaries)
#   (d) Positive-control panel: up in tumours in both datasets AND
#       subtype-neutral in both datasets AND well expressed
#
# Input : outputs of scripts 03, 07 and 08
# Output: results/figures/KM_PSCA_GSE71729.png
#         results/tables/GSE71729_limma_basal_vs_classical.csv
#         results/tables/positive_control_panel.csv
# ============================================================

setwd("C:/Sanjay/mPDAC Project")
while (!is.null(dev.list())) dev.off()

library(Biobase)
library(limma)
library(org.Hs.eg.db)
library(DESeq2)
library(survival)
library(survminer)

eset    <- readRDS("data/processed/gse71729_eset.rds")
samples <- readRDS("data/processed/gse71729_samples.rds")
expr    <- exprs(eset)
stopifnot(all(colnames(expr) == rownames(samples)))

# ---- (a) PSCA survival, GSE71729 primaries ----
# Verdict rules set beforehand: replicates = HR > 1 with p < 0.05;
# inconclusive = same direction, not significant; fails = no effect / opposite.
use <- !is.na(samples$surv_mo) & !is.na(samples$death) & samples$tissue %in% "Primary"
d <- samples[use, ]
d$psca <- as.numeric(expr["PSCA", use])
nrow(d); sum(d$death)                                          # 125 patients, 84 deaths

d$psca_group <- ifelse(d$psca > median(d$psca), "High", "Low")
survdiff(Surv(surv_mo, death) ~ psca_group, data = d)          # p = 0.3

fit_geo <- survfit(Surv(surv_mo, death) ~ psca_group, data = d)
png("results/figures/KM_PSCA_GSE71729.png", width = 9, height = 8, units = "in", res = 300)
print(ggsurvplot(fit_geo, data = d, pval = TRUE, risk.table = TRUE, conf.int = TRUE,
                 xlab = "Time (months)", ylab = "Survival Probability",
                 title = "PSCA Expression and Survival - GSE71729 Primary Tumours",
                 legend.title = "PSCA Expression", legend.labs = c("High", "Low"),
                 palette = c("#E64B35", "#4DBBD5")))
dev.off()

summary(coxph(Surv(surv_mo, death) ~ psca, data = d))                   # HR 1.08, p = 0.29
d$subtype_name <- factor(ifelse(d$subtype == "2", "Basal", "Classical"),
                         levels = c("Classical", "Basal"))
summary(coxph(Surv(surv_mo, death) ~ psca + subtype_name, data = d))    # PSCA HR 1.03; basal HR 1.86
survdiff(Surv(surv_mo, death) ~ subtype_name, data = d)                 # basal worse, p = 0.005

# ---- (b) Is PSCA subtype-linked? ----
dds2 <- readRDS("data/processed/dds2_confirmed150.rds")
clin <- readRDS("data/processed/tcga_paad_clinical_filtered.rds")[colnames(dds2), ]
tum  <- dds2$condition == "Tumor"
psca_tcga <- log2(counts(dds2, normalized = TRUE)["PSCA", tum] + 1)
mof       <- clin$`paper_mRNA.Moffitt.clusters..All.150.Samples..1basal..2classical`[tum]
sub_tcga  <- factor(ifelse(mof == 1, "Basal", "Classical"), levels = c("Classical", "Basal"))
tapply(psca_tcga, sub_tcga, median)
wilcox.test(psca_tcga ~ sub_tcga)                              # higher in basal, p = 0.005

prim     <- samples$tissue %in% "Primary"
psca_geo <- expr["PSCA", prim]
sub_geo  <- factor(ifelse(samples$subtype[prim] == "2", "Basal", "Classical"),
                   levels = c("Classical", "Basal"))
tapply(psca_geo, sub_geo, median)
wilcox.test(psca_geo ~ sub_geo)                                # higher in basal, p = 0.07

# ---- (c) GSE71729 basal vs classical (primaries) ----
xb <- expr[, prim]
keep_g <- rowSums(xb > 1) >= 21
fit_bc <- eBayes(lmFit(xb[keep_g, ], model.matrix(~ sub_geo)))
tt_bc  <- topTable(fit_bc, coef = "sub_geoBasal", number = Inf)
summary(decideTests(fit_bc, p.value = 0.05, lfc = 1))          # 40 up / 98 down in basal
tt_bc$symbol <- alias2SymbolTable(rownames(tt_bc), species = "Hs")
tt_bc <- tt_bc[!is.na(tt_bc$symbol) & !duplicated(tt_bc$symbol), ]
rownames(tt_bc) <- tt_bc$symbol
write.csv(tt_bc, "results/tables/GSE71729_limma_basal_vs_classical.csv")

# ---- (d) Positive-control panel ----
tcga_tn <- read.csv("results/tables/deseq2_tumor_vs_normal_confirmed150.csv", row.names = 1)
geo_tn  <- read.csv("results/tables/GSE71729_limma_tumor_vs_normal.csv",     row.names = 1)
tcga_bc <- read.csv("results/tables/deseq2_basal_vs_classical.csv",          row.names = 1)
geo_bc  <- tt_bc

g <- Reduce(intersect, list(rownames(tcga_tn), rownames(geo_tn), rownames(tcga_bc), rownames(geo_bc)))
f1 <- !is.na(tcga_tn[g, "padj"]) & tcga_tn[g, "padj"] < 0.05 & tcga_tn[g, "log2FoldChange"] > 2
f2 <- geo_tn[g, "adj.P.Val"] < 0.05 & geo_tn[g, "logFC"] > 1
f3 <- abs(tcga_bc[g, "log2FoldChange"]) < 0.5
f4 <- tcga_tn[g, "baseMean"] >= 100
f5 <- abs(geo_bc[g, "logFC"]) < 0.5
c(up_both = sum(f1 & f2), plus_TCGA_neutral = sum(f1 & f2 & f3),
  plus_expressed = sum(f1 & f2 & f3 & f4), all_five = sum(f1 & f2 & f3 & f4 & f5))   # 78 / 24 / 23 / 20

panel <- g[f1 & f2 & f3 & f4 & f5]
panel_df <- data.frame(
  TCGA_log2FC      = round(tcga_tn[panel, "log2FoldChange"], 2),
  GEO_logFC        = round(geo_tn[panel, "logFC"], 2),
  TCGA_basal_class = round(tcga_bc[panel, "log2FoldChange"], 2),
  GEO_basal_class  = round(geo_bc[panel, "logFC"], 2),
  baseMean         = round(tcga_tn[panel, "baseMean"]),
  row.names        = panel
)
panel_df <- panel_df[order(-panel_df$GEO_logFC), ]
panel_df
write.csv(panel_df, "results/tables/positive_control_panel.csv")

# Filter check: PSCA should fail neutrality
c(TCGA_bc = tcga_bc["PSCA", "log2FoldChange"], GEO_bc = geo_bc["PSCA", "logFC"])   # 0.68 / 0.51
