# ============================================================
# 07 - Basal vs classical (Moffitt) within the 150 confirmed tumours
#
# Tumour-only comparison: avoids the 4-normal limitation and asks what
# separates the two PDAC cell states. Classical = reference, so a positive
# log2FC means higher in basal tumours.
#
# Caveat: Moffitt labels were derived from expression of these samples, so
# subtype-defining genes are expected to rank highly (circularity).
#
# Input : data/processed/dds2_confirmed150.rds, tcga_paad_clinical_filtered.rds
# Output: results/tables/deseq2_basal_vs_classical.csv
#         results/tables/up_in_basal.csv            (280 genes)
#         results/tables/up_in_classical.csv        (535 genes)
#         results/tables/GO_BP_up_in_basal.csv, GO_BP_up_in_classical.csv
#         results/figures/volcano_basal_vs_classical.png
#         data/processed/dds_basal_vs_classical.rds
# ============================================================

setwd("C:/Sanjay/mPDAC Project")

library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)

dds2              <- readRDS("data/processed/dds2_confirmed150.rds")
clinical_filtered <- readRDS("data/processed/tcga_paad_clinical_filtered.rds")

# 1. Tumours only: raw counts + matching clinical rows
tum <- dds2$condition == "Tumor"
counts_t <- counts(dds2)[, tum]
clin_t <- clinical_filtered[colnames(counts_t), ]
stopifnot(all(rownames(clin_t) == colnames(counts_t)))

# 2. Subtype labels (1 = basal, 2 = classical)
mof <- clin_t$`paper_mRNA.Moffitt.clusters..All.150.Samples..1basal..2classical`
table(mof, useNA = "ifany")   # 65 basal, 85 classical
info_t <- data.frame(
  subtype = factor(ifelse(mof == 1, "Basal", "Classical"), levels = c("Classical", "Basal")),
  row.names = colnames(counts_t)
)
keep_s <- !is.na(info_t$subtype)

# 3. DESeq2
dds_sub <- DESeqDataSetFromMatrix(countData = counts_t[, keep_s],
                                  colData   = info_t[keep_s, , drop = FALSE],
                                  design    = ~ subtype)
dds_sub <- DESeq(dds_sub)
saveRDS(dds_sub, "data/processed/dds_basal_vs_classical.rds")

res_sub <- results(dds_sub, contrast = c("subtype", "Basal", "Classical"), alpha = 0.05)
summary(res_sub)   # padj < 0.05: 2,130 higher in basal, 2,529 higher in classical

# 4. High-confidence genes (padj < 0.05, |log2FC| > 1)
res_sub_df <- as.data.frame(res_sub)
sig_sub <- res_sub_df[!is.na(res_sub_df$padj) & res_sub_df$padj < 0.05 &
                      abs(res_sub_df$log2FoldChange) > 1, ]
up_basal     <- sig_sub[sig_sub$log2FoldChange >  1, ]
up_classical <- sig_sub[sig_sub$log2FoldChange < -1, ]
nrow(up_basal); nrow(up_classical)   # 280 / 535

head(up_basal[order(up_basal$padj), c("log2FoldChange", "padj")], 15)
head(up_classical[order(up_classical$padj), c("log2FoldChange", "padj")], 15)

write.csv(res_sub_df,   "results/tables/deseq2_basal_vs_classical.csv")
write.csv(up_basal,     "results/tables/up_in_basal.csv")
write.csv(up_classical, "results/tables/up_in_classical.csv")

# 5. Purity check: does one subtype simply contain more tumour cells?
purity <- clin_t$paper_ABSOLUTE.Purity[keep_s]
tapply(purity, info_t$subtype[keep_s], median, na.rm = TRUE)   # 0.315 basal vs 0.340 classical
wilcox.test(purity ~ info_t$subtype[keep_s])                    # p = 0.77, no difference

# 6. Overlap with the Moffitt et al. 2015 subtype signature (49 genes; WikiPathways WP5390)
moffitt_genes <- c("MYO1A","TNS4","CYP3A7","SPRR3","TFF2","SCEL","KRT17","CST6","FGFBP1",
  "VGLL1","LYZ","AGR2","CEACAM6","PRR15L","LY6D","KRT6C","SERPINB4","ANXA8L1","UCA1",
  "PLA2G10","SPRR1B","ST6GALNAC1","S100A2","KRT15","CTSV","KRT6A","VSIG2","TSPAN8",
  "LGALS4","FAM83A","REG4","LEMD1","AGR3","TFF3","SLC2A1","TFF1","KRT7","CDH17","DHRS9",
  "KRT20","SPINK4","CTSE","ANXA10","AREG","CLRN3","BTNL8","SERPINB3","FAM3D","GPR87")
sum(rownames(up_basal)     %in% moffitt_genes)   # 17
sum(rownames(up_classical) %in% moffitt_genes)   # 19

# 7. GO BP enrichment for each list
b_entrez <- bitr(rownames(up_basal),     fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
c_entrez <- bitr(rownames(up_classical), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_basal <- enrichGO(b_entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                     pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
go_class <- enrichGO(c_entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                     pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
head(go_basal@result[, c("Description", "GeneRatio", "p.adjust")], 10)   # epidermis development, keratinization
head(go_class@result[, c("Description", "GeneRatio", "p.adjust")], 10)   # digestion, xenobiotic metabolism
write.csv(go_basal@result, "results/tables/GO_BP_up_in_basal.csv",     row.names = FALSE)
write.csv(go_class@result, "results/tables/GO_BP_up_in_classical.csv", row.names = FALSE)

# 8. Volcano plot (Moffitt signature genes labelled in black, others by colour)
vd <- res_sub_df[!is.na(res_sub_df$padj), ]
vd$gene <- rownames(vd)
vd$group <- "Not significant"
vd$group[vd$padj < 0.05 & vd$log2FoldChange >  1] <- "Higher in basal"
vd$group[vd$padj < 0.05 & vd$log2FoldChange < -1] <- "Higher in classical"
lab <- c(rownames(up_basal)[order(up_basal$padj)][1:8],
         rownames(up_classical)[order(up_classical$padj)][1:8])

p <- ggplot(vd, aes(log2FoldChange, -log10(padj), color = group)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("Higher in basal" = "#E64B35",
                                "Higher in classical" = "#3C5488",
                                "Not significant" = "grey70")) +
  geom_text_repel(data = subset(vd, gene %in% lab),
                  aes(label = gene, fontface = ifelse(gene %in% moffitt_genes, "bold", "plain")),
                  size = 3.3, max.overlaps = 30, box.padding = 0.5,
                  color = "black", show.legend = FALSE) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Basal vs classical PDAC (Moffitt subtypes)",
       subtitle = "65 basal vs 85 classical tumours, TCGA-PAAD; bold labels = Moffitt signature genes",
       x = "log2 Fold Change (basal / classical)", y = "-log10 Adjusted P-value", color = NULL) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("results/figures/volcano_basal_vs_classical.png", plot = p,
       width = 10, height = 8, dpi = 300)
