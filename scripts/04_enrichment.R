# ============================================================
# 04 - GO Biological Process and KEGG over-representation
#      (up- and downregulated gene sets analysed separately)
# Input : results/tables/up/downregulated_DEGs_confirmed150.csv
# Output: results/tables/GO_BP_*_confirmed150.csv, KEGG_*_confirmed150.csv
#         results/figures/GO_BP_*_confirmed150.png, KEGG_*_confirmed150.png
# Caveat: uses clusterProfiler's default background (all annotated genes),
#         not the ~20,800 genes actually expressed in this dataset.
# ============================================================

setwd("C:/Sanjay/mPDAC Project")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

up2   <- read.csv("results/tables/upregulated_DEGs_confirmed150.csv",   row.names = 1)
down2 <- read.csv("results/tables/downregulated_DEGs_confirmed150.csv", row.names = 1)

up2_entrez   <- bitr(rownames(up2),   fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
down2_entrez <- bitr(rownames(down2), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

go_up2   <- enrichGO(up2_entrez$ENTREZID,   OrgDb = org.Hs.eg.db, ont = "BP",
                     pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
go_down2 <- enrichGO(down2_entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                     pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
kegg_up2   <- enrichKEGG(up2_entrez$ENTREZID,   organism = "hsa", pvalueCutoff = 0.05)
kegg_down2 <- enrichKEGG(down2_entrez$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)

write.csv(go_up2@result,     "results/tables/GO_BP_up_confirmed150.csv",   row.names = FALSE)
write.csv(go_down2@result,   "results/tables/GO_BP_down_confirmed150.csv", row.names = FALSE)
write.csv(kegg_up2@result,   "results/tables/KEGG_up_confirmed150.csv",    row.names = FALSE)
write.csv(kegg_down2@result, "results/tables/KEGG_down_confirmed150.csv",  row.names = FALSE)

plot_enrich <- function(obj, title, file) {
  p <- dotplot(obj, showCategory = 10) + ggtitle(title) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40))
  ggsave(file, plot = p, width = 11, height = 9, dpi = 300)
}
plot_enrich(go_up2,     "GO BP - Upregulated (150 confirmed PDAC vs normal)",   "results/figures/GO_BP_up_confirmed150.png")
plot_enrich(go_down2,   "GO BP - Downregulated (150 confirmed PDAC vs normal)", "results/figures/GO_BP_down_confirmed150.png")
plot_enrich(kegg_up2,   "KEGG - Upregulated (150 confirmed PDAC vs normal)",    "results/figures/KEGG_up_confirmed150.png")
plot_enrich(kegg_down2, "KEGG - Downregulated (150 confirmed PDAC vs normal)",  "results/figures/KEGG_down_confirmed150.png")

saveRDS(go_up2, "data/processed/go_up2_confirmed150.rds")   # used by script 05
