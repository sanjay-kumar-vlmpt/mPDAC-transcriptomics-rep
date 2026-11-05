# ============================================================
# 05 - Do the upregulated programs track Moffitt subtype or tumour purity?
# Input : data/processed/dds2_confirmed150.rds, go_up2_confirmed150.rds,
#         tcga_paad_clinical_filtered.rds
# ============================================================

setwd("C:/Sanjay/mPDAC Project")

library(DESeq2)

dds2              <- readRDS("data/processed/dds2_confirmed150.rds")
go_up2            <- readRDS("data/processed/go_up2_confirmed150.rds")
clinical_filtered <- readRDS("data/processed/tcga_paad_clinical_filtered.rds")

# Line clinical rows up with the samples in dds2 (150 tumours + 4 normals)
clin <- clinical_filtered[colnames(dds2), ]
stopifnot(all(rownames(clin) == colnames(dds2)))
tum  <- dds2$condition == "Tumor"
vsd_mat <- assay(vst(dds2, blind = TRUE))

moffitt <- clin$`paper_mRNA.Moffitt.clusters..All.150.Samples..1basal..2classical`
mof <- factor(ifelse(moffitt == 1, "Basal", "Classical"))[tum]
keep_m <- !is.na(mof)

# 1. SI expressed (>= 10 raw counts) vs Moffitt subtype
si_raw   <- counts(dds2)["SI", tum]
si_group <- ifelse(si_raw >= 10, "Expressed", "Not expressed")
table(si_group, mof)
fisher.test(table(si_group, mof))

# 2. Gut-type program score (mean VST of 8 upregulated gut genes)
gut <- c("SI", "AGR2", "ERN2", "TFF2", "MUC13", "MUC5B", "CAPN8", "AQP5")
score <- colMeans(vsd_mat[gut, tum])
tapply(score[keep_m], mof[keep_m], median)
wilcox.test(score[keep_m] ~ mof[keep_m])

# 3. Epidermal program score (genes from GO:0008544 epidermis development)
epi <- intersect(strsplit(go_up2@result["GO:0008544", "geneID"], "/")[[1]], rownames(vsd_mat))
epi_score <- colMeans(vsd_mat[epi, tum])
tapply(epi_score[keep_m], mof[keep_m], median)
wilcox.test(epi_score[keep_m] ~ mof[keep_m])

# 4. Immune genes vs tumour purity (bulk composition confound)
purity <- clin$paper_ABSOLUTE.Purity[tum]
cor.test(counts(dds2)["KLRF1", tum], purity, method = "spearman")
cor.test(counts(dds2)["PSCA",  tum], purity, method = "spearman")
