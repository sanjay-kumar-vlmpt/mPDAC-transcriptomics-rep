# ============================================================
# 03 - Differential expression: 150 confirmed PDAC tumours vs 4 normals
#
# Why "confirmed": sample_type does not separate PDAC from other pancreatic
# tumours. Tumours are kept only if they are in the TCGA 2017 paper's curated
# set (non-missing paper_Histological.type.by.DCC) -> 150 tumours.
#
# Output: results/tables/deseq2_tumor_vs_normal_confirmed150.csv
#         results/tables/upregulated_DEGs_confirmed150.csv   (664 genes)
#         results/tables/downregulated_DEGs_confirmed150.csv (426 genes)
#         results/figures/volcano_confirmed150.png
#         results/figures/heatmap_top50_confirmed150.png
#         data/processed/dds2_confirmed150.rds
# ============================================================

setwd("C:/Sanjay/mPDAC Project")
dir.create("results/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)

count_matrix_annotated <- readRDS("data/processed/tcga_paad_counts_annotated.rds")
clinical_filtered      <- readRDS("data/processed/tcga_paad_clinical_filtered.rds")

# 1. Build the confirmed cohort
is_tumor  <- clinical_filtered$sample_type == "Primary Tumor"
is_normal <- clinical_filtered$sample_type == "Solid Tissue Normal"
confirmed <- is_tumor & !is.na(clinical_filtered$paper_Histological.type.by.DCC)
keep <- confirmed | is_normal
table(clinical_filtered$sample_type[keep])   # 150 tumour + 4 normal

counts_conf <- count_matrix_annotated[, keep]
info_conf <- data.frame(
  condition = factor(ifelse(is_normal[keep], "Normal", "Tumor"),
                     levels = c("Normal", "Tumor")),   # Normal = reference
  row.names = colnames(counts_conf)
)

# 2. DESeq2
dds2 <- DESeqDataSetFromMatrix(countData = counts_conf, colData = info_conf,
                               design = ~ condition)
dds2 <- DESeq(dds2)
saveRDS(dds2, "data/processed/dds2_confirmed150.rds")

res2 <- results(dds2, contrast = c("condition", "Tumor", "Normal"), alpha = 0.05)
summary(res2)
res2_df <- as.data.frame(res2)
write.csv(res2_df, "results/tables/deseq2_tumor_vs_normal_confirmed150.csv")

# 3. High-confidence DEGs: padj < 0.05 and |log2FC| > 1
sig_df2 <- res2_df[!is.na(res2_df$padj) & res2_df$padj < 0.05 &
                   abs(res2_df$log2FoldChange) > 1, ]
up2   <- sig_df2[sig_df2$log2FoldChange >  1, ]
down2 <- sig_df2[sig_df2$log2FoldChange < -1, ]
nrow(up2); nrow(down2)   # 664 / 426
write.csv(up2,   "results/tables/upregulated_DEGs_confirmed150.csv")
write.csv(down2, "results/tables/downregulated_DEGs_confirmed150.csv")

# 4. Zero-count check: genes with |log2FC| > 10 are usually ~0 in the 4 normals
norm_counts <- counts(dds2, normalized = TRUE)
big <- rownames(sig_df2)[abs(sig_df2$log2FoldChange) > 10]
data.frame(gene          = big,
           log2FC        = round(sig_df2[big, "log2FoldChange"], 1),
           normal_mean   = round(rowMeans(norm_counts[big, is_normal[keep]]), 1),
           tumour_median = round(apply(norm_counts[big, !is_normal[keep]], 1, median), 1))
# -> SI and CT83: 0 in normals; fold change reflects on/off, not magnitude

# 5. Volcano plot
vd <- res2_df[!is.na(res2_df$padj), ]
vd$gene <- rownames(vd)
vd$significance <- "Not Significant"
vd$significance[vd$padj < 0.05 & vd$log2FoldChange >  1] <- "Upregulated"
vd$significance[vd$padj < 0.05 & vd$log2FoldChange < -1] <- "Downregulated"
top_genes <- c(rownames(up2)[order(up2$padj)][1:5],
               rownames(down2)[order(down2$padj)][1:5])

p_volcano <- ggplot(vd, aes(log2FoldChange, -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "#E64B35",
                                "Downregulated" = "#4DBBD5",
                                "Not Significant" = "grey60")) +
  geom_text_repel(data = subset(vd, gene %in% top_genes),
                  aes(label = gene), size = 3.5, max.overlaps = 30,
                  box.padding = 0.6, point.padding = 0.4,
                  show.legend = FALSE) +   # stops the stray "a" in the legend
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Differential Gene Expression in PDAC",
       subtitle = "150 confirmed tumours vs 4 normal - TCGA-PAAD",
       x = "log2 Fold Change", y = "-log10 Adjusted P-value", color = "Expression") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
ggsave("results/figures/volcano_confirmed150.png", plot = p_volcano,
       width = 10, height = 8, dpi = 300)

# 6. Heatmap: top 25 up + top 25 down (by padj), variance-stabilised, row z-scores
vsd_mat <- assay(vst(dds2, blind = TRUE))
top50 <- c(rownames(up2)[order(up2$padj)][1:25],
           rownames(down2)[order(down2$padj)][1:25])
hm <- t(scale(t(vsd_mat[top50, ])))

moffitt_all <- clinical_filtered$`paper_mRNA.Moffitt.clusters..All.150.Samples..1basal..2classical`[keep]
ann <- data.frame(
  Condition = info_conf$condition,
  Moffitt   = factor(ifelse(is.na(moffitt_all), "Not assigned",
                            ifelse(moffitt_all == 1, "Basal", "Classical"))),
  row.names = colnames(vsd_mat)
)

while (!is.null(dev.list())) dev.off()   # close any leaked graphics device first
png("results/figures/heatmap_top50_confirmed150.png",
    width = 10, height = 12, units = "in", res = 300)
pheatmap(hm, annotation_col = ann, show_colnames = FALSE, fontsize_row = 8,
         color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
         breaks = seq(-3, 3, length.out = 101),
         main = "Top 50 DEGs - 150 confirmed PDAC tumours vs 4 normal",
         border_color = NA)
dev.off()
