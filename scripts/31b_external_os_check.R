# 31b: Supplementary OS validation in GSE72970
Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8")
library(survival)
setwd(".")

df <- read.csv("results/tables/nomogram/gse72970_risk_scores.csv", stringsAsFactors = FALSE)
clin <- read.csv("results/tables/GSE72970_clinical_data.csv", stringsAsFactors = FALSE, row.names = 1)

df$os_event <- ifelse(as.numeric(clin[["os.censored.ch1"]]) == 0, 1, 0)
df$os_time <- as.numeric(clin[["os.ch1"]])
df$regimen <- clin[["regimen.ch1"]]

cat("=== Event Counts ===\n")
cat(sprintf("All (n=%d): PFS events=%d, OS events=%d\n", nrow(df), sum(df$pfs_event), sum(df$os_event)))
folfox <- df[df$regimen == "FOLFOX", ]
cat(sprintf("FOLFOX (n=%d): PFS events=%d, OS events=%d\n", nrow(folfox), sum(folfox$pfs_event), sum(folfox$os_event)))
folfiri <- df[df$regimen == "FOLFIRI", ]
cat(sprintf("FOLFIRI (n=%d): PFS events=%d, OS events=%d\n", nrow(folfiri), sum(folfiri$pfs_event), sum(folfiri$os_event)))

# OS all patients
df$risk_score_z <- scale(df$risk_score)[, 1]
cox_os <- coxph(Surv(os_time, os_event) ~ risk_score_z, data = df)
s_os <- summary(cox_os)
cat("\n=== OS (all patients) ===\n")
cat(sprintf("HR = %.4f (%.4f - %.4f), p = %.6f\n",
    s_os$conf.int[1, 1], s_os$conf.int[1, 3], s_os$conf.int[1, 4],
    s_os$coefficients[1, "Pr(>|z|)"]))
cat(sprintf("C-index = %.4f\n\n", concordance(cox_os)$concordance))

# OS FOLFOX
folfox$risk_score_z <- scale(folfox$risk_score)[, 1]
cox_os_f <- coxph(Surv(os_time, os_event) ~ risk_score_z, data = folfox)
s_os_f <- summary(cox_os_f)
cat("=== OS (FOLFOX only) ===\n")
cat(sprintf("N = %d\n", nrow(folfox)))
cat(sprintf("HR = %.4f (%.4f - %.4f), p = %.4f\n",
    s_os_f$conf.int[1, 1], s_os_f$conf.int[1, 3], s_os_f$conf.int[1, 4],
    s_os_f$coefficients[1, "Pr(>|z|)"]))
cat(sprintf("C-index = %.4f\n\n", concordance(cox_os_f)$concordance))

# PFS FOLFOX events table
cat("=== PFS events table (FOLFOX) ===\n")
print(table(folfox$pfs_event))
cat(sprintf("\nTotal FOLFOX PFS events: %d / %d\n", sum(folfox$pfs_event), nrow(folfox)))
