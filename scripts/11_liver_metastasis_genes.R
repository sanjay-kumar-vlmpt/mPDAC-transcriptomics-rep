# ============================================================
# 11 - Liver metastasis genes beyond subtype, tumour content and tissue background
#      (GSE71729: 25 liver metastases, 145 primaries, 27 normal liver, 46 normal pancreas)
#
# Logic: a bulk sample is a mixture (cancer cells + surrounding tissue), and a mixture
# cannot exceed both of its ingredients. A gene is kept only if it is
#   (1) different in liver mets vs primaries after adjusting for the basal-classical
#       score and a tumour-content score (mean z of the positive-control panel), AND
#   (2) different in the same direction in liver mets vs normal liver, AND
#   (3) (down genes) not a normal-pancreas tissue gene, since primaries contain pancreas.
# Islet genes (PPY, CPE) survive filter 3 because islets persist inside tumours,
# so they are removed by hand.
#
# Caveats: metastases likely from rapid autopsy (post-mortem ischaemia can induce
# hypoxia genes); patient structure unknown; bulk tissue cannot show which cells change.
#
# Input : outputs of scripts 08-10
# Output: results/tables/GSE71729_liver_met_genes.csv
#         results/tables/GSE71729_liver_met_stage_counts.csv
#         results/tables/GO_BP_liver_met_up.csv, GO_BP_liver_met_down.csv
# ============================================================

setwd("C:/Sanjay/mPDAC Project")

library(Biobase)
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)

eset    <- readRDS("data/processed/gse71729_eset.rds")
samples <- readRDS("data/processed/gse71729_samples.rds")
expr    <- exprs(eset)
sc      <- read.csv("results/tables/GSE71729_subtype_scores.csv",        row.names = 1)
panel   <- rownames(read.csv("results/tables/positive_control_panel.csv", row.names = 1))
geo_tn  <- read.csv("results/tables/GSE71729_limma_tumor_vs_normal.csv",  row.names = 1)
stopifnot(all(colnames(expr) == rownames(samples)))

# 1. Tumour-content score (17 of 20 panel genes carry current names in GEO)
panel_use <- intersect(panel, rownames(expr)); length(panel_use)          # 17
tum <- samples$source %in% c("Pancreas_Primary", "Liver_Metastasis")
content <- colMeans(t(scale(t(expr[panel_use, tum]))))
tapply(content, samples$source[tum], median)                             # 0.44 met vs 0.19 primary
wilcox.test(content ~ samples$source[tum])                               # p = 0.025

# 2. Liver mets vs primaries: naive and adjusted
x1 <- expr[, tum]
g1 <- factor(ifelse(samples$source[tum] == "Liver_Metastasis", "LiverMet", "Primary"),
             levels = c("Primary", "LiverMet"))
keep1    <- rowSums(x1 > 1) >= 13
subscore <- sc[colnames(x1), "score"]
t0 <- topTable(eBayes(lmFit(x1[keep1, ], model.matrix(~ g1))),
               coef = "g1LiverMet", number = Inf)
t1 <- topTable(eBayes(lmFit(x1[keep1, ], model.matrix(~ g1 + subscore + content))),
               coef = "g1LiverMet", number = Inf)

# 3. Liver mets vs normal liver
s2 <- samples$source %in% c("Liver_Metastasis", "Liver_Normal")
x2 <- expr[, s2]
g2 <- factor(ifelse(samples$source[s2] == "Liver_Metastasis", "LiverMet", "LiverNormal"),
             levels = c("LiverNormal", "LiverMet"))
keep2 <- rowSums(x2 > 1) >= 13
t2 <- topTable(eBayes(lmFit(x2[keep2, ], model.matrix(~ g2))),
               coef = "g2LiverMet", number = Inf)

# 4. Combine (mixture rule)
sig_up   <- function(t) rownames(t)[t$adj.P.Val < 0.05 & t$logFC >  1]
sig_down <- function(t) rownames(t)[t$adj.P.Val < 0.05 & t$logFC < -1]
up_final   <- intersect(sig_up(t1),   rownames(t2)[t2$adj.P.Val < 0.05 & t2$logFC > 0])
down_liver <- intersect(sig_down(t1), rownames(t2)[t2$adj.P.Val < 0.05 & t2$logFC < 0])

# 5. Remove normal-pancreas tissue genes from the down list
res <- data.frame(
  logFC_vs_primary = round(t1[c(up_final, down_liver), "logFC"], 2),
  padj_vs_primary  = signif(t1[c(up_final, down_liver), "adj.P.Val"], 2),
  logFC_vs_liver   = round(t2[c(up_final, down_liver), "logFC"], 2),
  padj_vs_liver    = signif(t2[c(up_final, down_liver), "adj.P.Val"], 2),
  row.names = c(up_final, down_liver)
)
res$symbol <- alias2SymbolTable(rownames(res), species = "Hs")
pancreas_genes <- rownames(geo_tn)[geo_tn$adj.P.Val < 0.05 & geo_tn$logFC < -1]   # 171
islet_genes    <- c("PPY", "PPY2P", "CPE")                                         # persist in tumours
res$pancreas_tissue <- res$logFC_vs_primary < 0 &
  ((!is.na(res$symbol) & res$symbol %in% c(pancreas_genes, islet_genes)) |
   rownames(res) %in% c("PPY", "PPY2", "CPE"))
res[res$pancreas_tissue, c("logFC_vs_primary", "logFC_vs_liver", "symbol")]
down_final <- rownames(res)[res$logFC_vs_primary < 0 & !res$pancreas_tissue]

stage_counts <- data.frame(
  stage = c("Naive (no adjustment)", "Adjusted for subtype + tumour content",
            "Beyond liver background", "Also excluding pancreas tissue genes"),
  up    = c(length(sig_up(t0)),   length(sig_up(t1)),   length(up_final),   length(up_final)),
  down  = c(length(sig_down(t0)), length(sig_down(t1)), length(down_liver), length(down_final))
)
stage_counts                                       # 445/608, 465/534, 32/87, 32/~79
write.csv(stage_counts, "results/tables/GSE71729_liver_met_stage_counts.csv", row.names = FALSE)

res <- res[order(-abs(res$logFC_vs_primary)), ]
head(res[res$logFC_vs_primary > 0, ], 20)
head(res[res$logFC_vs_primary < 0 & !res$pancreas_tissue, ], 20)
write.csv(res, "results/tables/GSE71729_liver_met_genes.csv")

# 6. GO enrichment of the final lists
go_list <- function(genes) {
  sym <- na.omit(alias2SymbolTable(genes, species = "Hs"))
  ent <- bitr(sym, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  enrichGO(ent$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
           pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
}
go_up   <- go_list(up_final)       # cut-offs relaxed to 1 so the full table is saved;
go_down <- go_list(down_final)     # judge significance by p.adjust < 0.05
head(go_up@result[,   c("Description", "GeneRatio", "p.adjust")], 10)   # hypoxia, padj ~0.11 (not significant)
head(go_down@result[, c("Description", "GeneRatio", "p.adjust")], 10)   # MHC class II, padj ~1e-15
write.csv(go_up@result,   "results/tables/GO_BP_liver_met_up.csv",   row.names = FALSE)
write.csv(go_down@result, "results/tables/GO_BP_liver_met_down.csv", row.names = FALSE)

# 7. Hypoxia-response genes (interpret with the autopsy caveat)
res[intersect(c("HK2", "VEGFA", "NDRG1", "SLC5A3", "SLC45A4"), rownames(res)),
    c("logFC_vs_primary", "padj_vs_primary", "logFC_vs_liver", "padj_vs_liver")]
