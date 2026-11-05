# ============================================================
# 06 - Survival analysis on the 150 confirmed tumours
#      Kaplan-Meier + log-rank (median split on raw counts) and Cox models
# Input : data/processed/tcga_paad_counts_annotated.rds, tcga_paad_clinical_filtered.rds
# Output: results/figures/KM_PSCA_confirmed150.png
# ============================================================

setwd("C:/Sanjay/mPDAC Project")
while (!is.null(dev.list())) dev.off()

library(survival)
library(survminer)

clinical_filtered      <- readRDS("data/processed/tcga_paad_clinical_filtered.rds")
count_matrix_annotated <- readRDS("data/processed/tcga_paad_counts_annotated.rds")

conf_idx <- clinical_filtered$sample_type == "Primary Tumor" &
            !is.na(clinical_filtered$paper_Histological.type.by.DCC)
conf_tumor <- clinical_filtered[conf_idx, ]
cm <- count_matrix_annotated[, conf_idx]
stopifnot(all(colnames(cm) == rownames(conf_tumor)))

# Time = days to death if dead, otherwise days to last follow-up (censored)
conf_tumor$time  <- ifelse(conf_tumor$vital_status == "Dead",
                           conf_tumor$days_to_death, conf_tumor$days_to_last_follow_up)
conf_tumor$event <- ifelse(conf_tumor$vital_status == "Dead", 1, 0)
table(conf_tumor$event)   # 82 deaths

# 1. Median-split log-rank tests for each candidate gene
for (g in c("PSCA", "SI", "KLRF1")) {
  expr <- as.numeric(cm[g, ])
  grp  <- ifelse(expr > median(expr), "High", "Low")
  cat("\n=====", g, "=====\n")
  print(table(grp))
  print(survdiff(Surv(conf_tumor$time, conf_tumor$event) ~ grp))
}
# Results: PSCA p = 0.06, SI p = 0.8, KLRF1 p = 0.4

# 2. PSCA Kaplan-Meier plot
psca <- as.numeric(cm["PSCA", ])
conf_tumor$psca_group <- ifelse(psca > median(psca), "High", "Low")
fit_conf <- survfit(Surv(time, event) ~ psca_group, data = conf_tumor)

png("results/figures/KM_PSCA_confirmed150.png", width = 9, height = 8, units = "in", res = 300)
print(ggsurvplot(fit_conf, data = conf_tumor, pval = TRUE,
                 risk.table = TRUE, conf.int = TRUE,
                 xlab = "Time (days)", ylab = "Survival Probability",
                 title = "PSCA Expression and Survival - 150 Confirmed PDAC Tumours",
                 legend.title = "PSCA Expression", legend.labs = c("High", "Low"),
                 palette = c("#E64B35", "#4DBBD5")))
dev.off()

# 3. Cox models, PSCA as continuous log2(count + 1)
conf_tumor$psca_log <- log2(psca + 1)
summary(coxph(Surv(time, event) ~ psca_log + age_at_index, data = conf_tumor))
# -> HR ~1.10 per log2 unit, p = 0.024

conf_tumor$moffitt_group <- ifelse(
  conf_tumor$`paper_mRNA.Moffitt.clusters..All.150.Samples..1basal..2classical` == 1,
  "Basal", "Classical")
summary(coxph(Surv(time, event) ~ psca_log + moffitt_group +
                paper_ABSOLUTE.Purity + age_at_index, data = conf_tumor))
# -> HR ~1.08, p = 0.065 (complete cases only)

# 4. Survival by Moffitt subtype
survdiff(Surv(time, event) ~ moffitt_group, data = conf_tumor)   # p = 0.08
