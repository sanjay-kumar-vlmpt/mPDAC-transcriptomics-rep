# ============================================================
# 12 - GSE21501 (Stratford 2010): processing and overlap check with GSE71729
#
# GSE21501: Agilent 4x44K two-colour arrays (GPL4133), UNC. Values are Lowess-normalised
# log2 ratios (sample Cy5 / Stratagene reference Cy3). All 132 samples are primary tumours;
# survival is available for 102 resected patients (Johns Hopkins, NorthShore, Northwestern).
#
# Why the overlap check: same lab and platform as GSE71729 (Moffitt 2015), which pooled
# samples from several sources. Without shared patient IDs, overlap is detected by expression:
# each gene is centred within its own dataset (this removes the ratio vs intensity difference,
# because the common reference is a constant per gene), samples are correlated across
# datasets, and a pair counts as the same patient if r > 0.6 AND survival months are identical.
#
# Result: 89 of 102 survival patients are confirmed duplicates of GSE71729 primaries, so
# GSE21501 is NOT an independent survival cohort; the ~10 remaining patients are too few.
#
# Output: data/processed/gse21501_eset.rds, gse21501_samples.rds, gse21501_genes.rds
#         results/tables/GSE21501_overlap_with_GSE71729.csv
#         results/figures/GSE21501_overlap_histogram.png
# ============================================================

setwd("C:/Sanjay/mPDAC Project")
while (!is.null(dev.list())) dev.off()

library(GEOquery)
library(Biobase)
options(timeout = 600)

# 1. Download and inspect
gse   <- getGEO("GSE21501", GSEMatrix = TRUE, destdir = "data/raw/GEO")
eset2 <- gse[[1]]
dim(exprs(eset2))                    # 45,220 probes x 132 samples
annotation(eset2)                    # GPL4133
pd2 <- pData(eset2); fd2 <- fData(eset2); ex2 <- exprs(eset2)
pd2$data_processing[1]               # log(2) ratio
summary(as.vector(ex2))              # centred on 0, -11 to +11.7
saveRDS(eset2, "data/processed/gse21501_eset.rds")

# 2. Sample table (survival time is in months: range 1-59, median 14.5)
s2 <- data.frame(
  id       = pd2$geo_accession,
  title    = pd2$title,
  provider = pd2$biomaterial_provider_ch2,
  os_time  = as.numeric(pd2$`os time:ch2`),
  os_event = as.numeric(pd2$`os event:ch2`),
  t_stage  = pd2$`t stage:ch2`,
  n_stage  = pd2$`n stage:ch2`,
  row.names = pd2$geo_accession
)
table(s2$provider, has_survival = !is.na(s2$os_time))   # 102 with survival; autopsy + UNC have none
s2$flag <- ifelse(grepl("N Replicate", s2$title), "possible normal",
           ifelse(grepl("IPMN", s2$title), "IPMN", ""))
table(s2$flag[!is.na(s2$os_time)])                       # P271N, P188T IPMN
saveRDS(s2, "data/processed/gse21501_samples.rds")

# 3. Probes -> genes: drop control spots and unannotated probes, average per gene
table(fd2$CONTROL_TYPE)
keep_p <- fd2$CONTROL_TYPE == "FALSE" & !is.na(fd2$GENE_SYMBOL) & fd2$GENE_SYMBOL != ""
sym    <- fd2$GENE_SYMBOL[keep_p]
g_ex   <- rowsum(ex2[keep_p, ], sym) / as.vector(table(sym)[sort(unique(sym))])
dim(g_ex)                                                # 19,749 genes (same as GSE71729)
saveRDS(g_ex, "data/processed/gse21501_genes.rds")

# 4. Overlap with GSE71729 primaries
eset    <- readRDS("data/processed/gse71729_eset.rds")
samples <- readRDS("data/processed/gse71729_samples.rds")
expr    <- exprs(eset)
prim    <- samples$tissue %in% "Primary"

genes <- intersect(rownames(g_ex), rownames(expr))
a <- g_ex[genes, !is.na(s2$os_time)]                    # 102 GSE21501 samples with survival
b <- expr[genes, prim]                                  # 145 GSE71729 primaries
keep_g <- apply(a, 1, sd) > 0 & apply(b, 1, sd) > 0
a <- a[keep_g, ] - rowMeans(a[keep_g, ])
b <- b[keep_g, ] - rowMeans(b[keep_g, ])
r <- cor(a, b)

summary(as.vector(r))                                   # unrelated tumours: median ~0
best   <- apply(r, 1, max)
best_j <- colnames(r)[apply(r, 1, which.max)]

overlap <- data.frame(
  GSE21501_title  = s2[colnames(a), "title"],
  best_r          = round(best, 2),
  GSE71729_match  = best_j,
  months_GSE21501 = s2[colnames(a), "os_time"],
  months_GSE71729 = samples[best_j, "surv_mo"],
  row.names       = colnames(a)
)
overlap$same_months <- !is.na(overlap$months_GSE71729) &
                       overlap$months_GSE21501 == overlap$months_GSE71729
overlap$confirmed_duplicate <- overlap$best_r > 0.6 & overlap$same_months
table(strong_match = overlap$best_r > 0.6, same_survival = overlap$same_months)
sum(overlap$confirmed_duplicate)                        # 89 of 102
overlap[!overlap$confirmed_duplicate, ]
write.csv(overlap[order(-overlap$best_r), ], "results/tables/GSE21501_overlap_with_GSE71729.csv")

png("results/figures/GSE21501_overlap_histogram.png", width = 8, height = 5, units = "in", res = 300)
hist(best, breaks = 30, col = "grey80",
     main = "Best match of each GSE21501 sample in GSE71729",
     xlab = "Correlation with closest GSE71729 primary (gene-centred)")
abline(v = 0.6, lty = 2)
dev.off()
