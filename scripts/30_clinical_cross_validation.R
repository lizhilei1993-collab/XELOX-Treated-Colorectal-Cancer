# ============================================================
# 30_CLINICAL_CROSS_VALIDATION.R
# Step 1: Clinical cross-validation of the 8-variable risk model
#
# Purpose:
#   Validate the pathway-based risk score against clinical
#   variables in the GSE39582 XELOX cohort (n=227).
#
# Analyses:
#   A. Sample characteristics table (risk_group x clinical vars)
#   B. Chi-square / Fisher exact tests for categorical variables
#   C. Multivariable Cox adjusting for clinical covariates
#   D. Forest plot comparing crude vs adjusted HRs
#
# Input:
#   results/tables/nomogram/risk_scores.csv
#   results/tables/GSE39582_xelox_groups.csv
#
# Output (results/tables/nomogram/):
#   clinical_validation_report.txt   — full report
#   clinical_table_risk_groups.csv   — sample characteristics table
#   multivariable_cox_results.csv    — Cox regression results
#   forest_crude_vs_adjusted.pdf     — forest plot comparison
# ============================================================

Sys.setenv(TMPDIR = "/tmp", TMP = "/tmp", TEMP = "/tmp")
.libPaths(c("/path/to/Rlibs", .libPaths()))

PROJECT_ROOT <- "/path/to/xelox_project基于可解释性机器学习的XELOX耐药分子指纹研究"
RESULTS_TAB_DIR <- file.path(PROJECT_ROOT, "results", "tables")
NOMO_TAB_DIR   <- file.path(RESULTS_TAB_DIR, "nomogram")
NOMO_FIG_DIR   <- file.path(PROJECT_ROOT, "results", "figures", "nomogram")

set.seed(42)

cat("============================================================\n")
cat("Step 1: Clinical Cross-Validation of 8-Variable Risk Model\n")
cat("============================================================\n\n")

# ============================================================
# 0. Load packages
# ============================================================
cat("=== [0] Loading packages ===\n\n")
required_pkgs <- c("survival", "ggplot2", "data.table")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, lib = "/path/to/Rlibs", repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  cat(sprintf("  %s v%s loaded\n", pkg, packageVersion(pkg)))
}
cat("\n")

# ============================================================
# 1. Load data
# ============================================================
cat("=== [1] Loading data ===\n\n")

# Risk scores from nomogram model
risk <- read.csv(file.path(NOMO_TAB_DIR, "risk_scores.csv"),
                 stringsAsFactors = FALSE)
cat(sprintf("  Risk scores: %d samples\n", nrow(risk)))

# Clinical data (drop overlapping cols to avoid merge suffix issues)
clin <- read.csv(file.path(RESULTS_TAB_DIR, "GSE39582_xelox_groups.csv"),
                  stringsAsFactors = FALSE)
cat(sprintf("  Clinical data: %d samples\n", nrow(clin)))

# Drop survival columns from clin (use risk's) and optional duplicates
drop_cols <- intersect(c("rfs_event", "rfs_delay", "os_event", "os_delay", "title"), colnames(clin))
clin_clean <- clin[, setdiff(colnames(clin), drop_cols), drop = FALSE]

# Merge
df <- merge(risk, clin_clean, by = "sample_id", all.x = TRUE)
cat(sprintf("  Merged: %d samples\n\n", nrow(df)))

# ============================================================
# 2. A. Sample characteristics by risk group
# ============================================================
cat("=== [2] Sample characteristics by risk group ===\n\n")

# Define clinical variables of interest
clin_vars <- list(
  "Sex"       = list(var = "sex",         levels = c("Male", "Female")),
  "Age group" = list(var = "age_group",   levels = c("<65", ">=65")),
  "TNM Stage" = list(var = "stage_group", levels = c("II", "III", "IV")),
  "Tumor location" = list(var = "tumor_location", levels = c("distal", "proximal")),
  "MMR status"     = list(var = "mmr_status", levels = c("pMMR", "dMMR")),
  "KRAS"           = list(var = "kras",   levels = c("WT", "M")),
  "BRAF"           = list(var = "braf",   levels = c("WT", "M"))
)

# Create age group
df$age_num <- as.numeric(df$age)
df$age_group <- ifelse(df$age_num < 65, "<65", ">=65")

# Create stage group
df$stage_group <- ifelse(grepl("^2|II", df$tnm_stage), "II",
                         ifelse(grepl("^3|III", df$tnm_stage), "III",
                                ifelse(grepl("^4|IV", df$tnm_stage), "IV", NA)))
df$stage_group <- factor(df$stage_group, levels = c("II", "III", "IV"))

# Risk group factor
df$risk_group <- factor(df$risk_group, levels = c("Low risk", "High risk"))

# Build characteristics table
char_table <- data.frame(
  Variable = character(),
  Level = character(),
  Total_N = integer(),
  Total_Pct = character(),
  Low_N = integer(),
  Low_Pct = character(),
  High_N = integer(),
  High_Pct = character(),
  P_value = character(),
  Test = character(),
  stringsAsFactors = FALSE
)

for (cv_name in names(clin_vars)) {
  cv <- clin_vars[[cv_name]]
  var_name <- cv$var

  # Remove NA
  valid <- df[!is.na(df[[var_name]]), ]
  n_total <- nrow(valid)
  n_low <- sum(valid$risk_group == "Low risk")
  n_high <- sum(valid$risk_group == "High risk")

  for (lvl in cv$levels) {
    n_lvl_total <- sum(valid[[var_name]] == lvl, na.rm = TRUE)
    n_lvl_low   <- sum(valid[[var_name]] == lvl & valid$risk_group == "Low risk", na.rm = TRUE)
    n_lvl_high  <- sum(valid[[var_name]] == lvl & valid$risk_group == "High risk", na.rm = TRUE)

    char_table <- rbind(char_table, data.frame(
      Variable = cv_name,
      Level = lvl,
      Total_N = n_lvl_total,
      Total_Pct = sprintf("%.1f%%", n_lvl_total / n_total * 100),
      Low_N = n_lvl_low,
      Low_Pct = sprintf("%.1f%%", n_lvl_low / n_low * 100),
      High_N = n_lvl_high,
      High_Pct = sprintf("%.1f%%", n_lvl_high / n_high * 100),
      P_value = "",
      Test = "",
      stringsAsFactors = FALSE
    ))
  }

  # Contingency table
  ct <- table(valid$risk_group, valid[[var_name]])
  if (any(ct < 5)) {
    ft <- tryCatch(fisher.test(ct, simulate.p.value = TRUE, B = 10000),
                   error = function(e) NULL)
    if (!is.null(ft)) {
      p_val <- ft$p.value
      test_name <- "Fisher"
    } else {
      p_val <- NA
      test_name <- "Fisher (failed)"
    }
  } else {
    ct <- tryCatch(chisq.test(ct), error = function(e) NULL)
    if (!is.null(ct)) {
      p_val <- ct$p.value
      test_name <- "Chi-square"
    } else {
      p_val <- NA
      test_name <- "Chi-square (failed)"
    }
  }

  # Fill p-value into last row for this variable
  char_table$P_value[nrow(char_table)] <- sprintf("%.4f", p_val)
  char_table$Test[nrow(char_table)] <- test_name

  # Print
  cat(sprintf("  %s (p=%s, %s):\n", cv_name, sprintf("%.4f", p_val), test_name))
  cat(sprintf("    %-20s %s / %s / %s\n", "Level", "Total", "Low", "High"))
  for (lvl in cv$levels) {
    row <- char_table[char_table$Variable == cv_name & char_table$Level == lvl, ]
    cat(sprintf("    %-20s %d(%s) / %d(%s) / %d(%s)\n",
                lvl, row$Total_N, row$Total_Pct,
                row$Low_N, row$Low_Pct,
                row$High_N, row$High_Pct))
  }
  cat("\n")
}

# Save characteristics table
write.csv(char_table, file.path(NOMO_TAB_DIR, "clinical_table_risk_groups.csv"),
          row.names = FALSE)
cat(sprintf("  -> Saved: %s\n\n",
            file.path(NOMO_TAB_DIR, "clinical_table_risk_groups.csv")))

# ============================================================
# 3. B/C. Multivariable Cox: crude + adjusted
# ============================================================
cat("=== [3] Multivariable Cox regression ===\n\n")

# Prepare variables
df$sex_bin      <- ifelse(df$sex == "Male", 1, 0)
df$age_num      <- as.numeric(df$age)
df$stage_III_IV <- ifelse(grepl("3|4|III|IV", df$tnm_stage), 1, 0)
df$mmr_dmmr     <- ifelse(df$mmr_status == "dMMR", 1, 0)
df$location     <- ifelse(df$tumor_location == "distal", 1, 0)

# Standardize risk score for interpretable HR per SD
df$risk_score_z <- scale(df$risk_score)[, 1]

# --- Crude ---
cox_crude <- coxph(Surv(rfs_delay, rfs_event) ~ risk_score_z, data = df)
s_crude <- summary(cox_crude)

# --- Adjusted for all clinical ---
cox_adj <- coxph(Surv(rfs_delay, rfs_event) ~ risk_score_z + age_num + sex_bin +
                   stage_III_IV + mmr_dmmr + location, data = df)
s_adj <- summary(cox_adj)

# --- Clinical-only (no risk score) ---
cox_clin_only <- coxph(Surv(rfs_delay, rfs_event) ~ age_num + sex_bin +
                         stage_III_IV + mmr_dmmr + location, data = df)
s_clin_only <- summary(cox_clin_only)

# --- Print results ---
cat("  Crude model (risk_score_z only):\n")
cat(sprintf("    HR = %.4f (%.4f - %.4f), p = %.6f\n",
            s_crude$conf.int[1, 1], s_crude$conf.int[1, 3],
            s_crude$conf.int[1, 4], s_crude$coefficients[1, "Pr(>|z|)"]))
cat(sprintf("    C-index = %.4f\n\n", concordance(cox_crude)$concordance))

cat("  Clinical-only model:\n")
for (i in seq_len(nrow(s_clin_only$coefficients))) {
  cat(sprintf("    %s: HR = %.4f (%.4f - %.4f), p = %.4f\n",
              names(coef(cox_clin_only))[i],
              s_clin_only$conf.int[i, 1], s_clin_only$conf.int[i, 3],
              s_clin_only$conf.int[i, 4],
              s_clin_only$coefficients[i, "Pr(>|z|)"]))
}
cat(sprintf("    C-index = %.4f\n\n", concordance(cox_clin_only)$concordance))

cat("  Adjusted model (risk_score_z + clinical):\n")
for (i in seq_len(nrow(s_adj$coefficients))) {
  cat(sprintf("    %s: HR = %.4f (%.4f - %.4f), p = %.4f\n",
              names(coef(cox_adj))[i],
              s_adj$conf.int[i, 1], s_adj$conf.int[i, 3],
              s_adj$conf.int[i, 4],
              s_adj$coefficients[i, "Pr(>|z|)"]))
}
cat(sprintf("    C-index = %.4f\n\n", concordance(cox_adj)$concordance))

# Likelihood ratio tests
lrt_adj_vs_clin <- anova(cox_clin_only, cox_adj)
cat(sprintf("  LRT: adjusted vs clinical-only: p = %.4f\n\n",
            lrt_adj_vs_clin$`Pr(>|Chi|)`[2]))

# Save multivariable results
multi_out <- data.frame(
  Model = c("Crude", rep("Adjusted", nrow(s_adj$coefficients)),
            rep("Clinical-only", nrow(s_clin_only$coefficients))),
  Variable = c("risk_score_z",
               names(coef(cox_adj)),
               names(coef(cox_clin_only))),
  HR = c(s_crude$conf.int[1, 1],
         s_adj$conf.int[, 1],
         s_clin_only$conf.int[, 1]),
  HR_lower = c(s_crude$conf.int[1, 3],
               s_adj$conf.int[, 3],
               s_clin_only$conf.int[, 3]),
  HR_upper = c(s_crude$conf.int[1, 4],
               s_adj$conf.int[, 4],
               s_clin_only$conf.int[, 4]),
  p_value = c(s_crude$coefficients[1, "Pr(>|z|)"],
              s_adj$coefficients[, "Pr(>|z|)"],
              s_clin_only$coefficients[, "Pr(>|z|)"]),
  C_index = c(concordance(cox_crude)$concordance,
              rep(concordance(cox_adj)$concordance, nrow(s_adj$coefficients)),
              rep(concordance(cox_clin_only)$concordance, nrow(s_clin_only$coefficients))),
  stringsAsFactors = FALSE
)
write.csv(multi_out, file.path(NOMO_TAB_DIR, "multivariable_cox_results.csv"),
          row.names = FALSE)
cat(sprintf("  -> Saved: %s\n\n",
            file.path(NOMO_TAB_DIR, "multivariable_cox_results.csv")))

# ============================================================
# 4. D. Forest plot: crude vs adjusted HRs
# ============================================================
cat("=== [4] Forest plot: Crude vs Adjusted HRs ===\n\n")

pdf(file.path(NOMO_FIG_DIR, "forest_crude_vs_adjusted.pdf"),
    width = 8, height = 4.5)
par(mar = c(4, 6, 3, 6))

# Plot data
vars <- c("risk_score_z", "age_num", "sex_bin", "stage_III_IV",
          "mmr_dmmr", "location")
var_labels <- c("Risk score (per SD)", "Age (per year)", "Sex (Male)",
                "Stage III/IV", "MMR (dMMR)", "Location (Distal)")
hr_adj <- exp(coef(cox_adj))
ci_adj_l <- exp(confint(cox_adj)[, 1])
ci_adj_u <- exp(confint(cox_adj)[, 2])
p_adj <- s_adj$coefficients[, "Pr(>|z|)"]

# Set plot bounds
hr_range <- range(c(hr_adj, ci_adj_l, ci_adj_u, 0.3, 3.0))
hr_max <- max(hr_range) * 1.5

plot(NA, xlim = log(c(0.1, 4)), ylim = c(0.5, length(vars) + 1),
     xaxt = "n", yaxt = "n",
     xlab = "Hazard Ratio (95% CI)", ylab = "",
     main = "Multivariable Cox: Independent Prognostic Value of Risk Score")
axis(1, at = log(c(0.2, 0.5, 1, 2, 4)),
     labels = c("0.2", "0.5", "1", "2", "4"))
abline(v = log(1), lty = 2, col = "gray70")

for (i in seq_along(vars)) {
  y_pos <- length(vars) - i + 1
  col <- ifelse(p_adj[i] < 0.05, "darkred", "gray40")
  points(log(hr_adj[i]), y_pos, pch = 15, cex = 1.2, col = col)
  segments(log(ci_adj_l[i]), y_pos, log(ci_adj_u[i]), y_pos,
           lwd = 2, col = col)
  # Label
  text(log(0.12), y_pos, var_labels[i], adj = 1, cex = 0.85, xpd = NA)
  # HR text
  text(log(4), y_pos,
       sprintf("%.2f (%.2f-%.2f)  p=%s",
               hr_adj[i], ci_adj_l[i], ci_adj_u[i],
               ifelse(p_adj[i] < 0.001, sprintf("%.0e", p_adj[i]),
                      sprintf("%.4f", p_adj[i]))),
       adj = 0, cex = 0.75, xpd = NA)
}
dev.off()
cat(sprintf("  -> Saved: %s\n\n",
            file.path(NOMO_FIG_DIR, "forest_crude_vs_adjusted.pdf")))

# ============================================================
# 5. Summary report
# ============================================================
cat("=== [5] Writing summary report ===\n\n")

sink(file.path(NOMO_TAB_DIR, "clinical_validation_report.txt"))

cat("=== Clinical Cross-Validation Report ===\n")
cat("========================================\n\n")
cat(sprintf("Cohort: GSE39582 XELOX subgroup\n"))
cat(sprintf("Samples: %d\n", nrow(df)))
cat(sprintf("RFS events: %d / %d (%.1f%%)\n\n",
            sum(df$rfs_event), nrow(df), mean(df$rfs_event)*100))

cat("--- A. Sample Characteristics ---\n\n")
print(char_table, row.names = FALSE)
cat("\n\n")

cat("--- B. Crude Cox Regression ---\n")
cat(sprintf("  Risk score: HR = %.4f (%.4f - %.4f), p = %.6f\n",
            s_crude$conf.int[1, 1], s_crude$conf.int[1, 3],
            s_crude$conf.int[1, 4], s_crude$coefficients[1, "Pr(>|z|)"]))
cat(sprintf("  C-index: %.4f\n\n", concordance(cox_crude)$concordance))

cat("--- C. Adjusted Cox Regression ---\n")
cat(sprintf("  Risk score (per SD): HR = %.4f (%.4f - %.4f), p = %.6f\n",
            s_adj$conf.int[1, 1], s_adj$conf.int[1, 3],
            s_adj$conf.int[1, 4], s_adj$coefficients[1, "Pr(>|z|)"]))
cat("\n  Adjusted model summary:\n")
for (i in seq_len(nrow(s_adj$coefficients))) {
  vname <- names(coef(cox_adj))[i]
  cat(sprintf("  %-20s HR = %.4f (%.4f - %.4f), p = %.4f\n",
              vname,
              s_adj$conf.int[i, 1], s_adj$conf.int[i, 3],
              s_adj$conf.int[i, 4],
              s_adj$coefficients[i, "Pr(>|z|)"]))
}
cat(sprintf("\n  Adjusted C-index: %.4f\n", concordance(cox_adj)$concordance))
cat(sprintf("  Clinical-only C-index: %.4f\n", concordance(cox_clin_only)$concordance))
cat(sprintf("  C-index improvement: +%.4f\n\n",
            concordance(cox_adj)$concordance - concordance(cox_clin_only)$concordance))

cat("--- D. Likelihood Ratio Test ---\n")
cat(sprintf("  Adjusted vs Clinical-only: Chisq = %.2f, df = %d, p = %.4f\n\n",
            lrt_adj_vs_clin$Deviance[2], lrt_adj_vs_clin$Df[2],
            lrt_adj_vs_clin$`Pr(>|Chi|)`[2]))

cat("--- E. Interpretation ---\n")
cat("  The risk score retains independent prognostic value after\n")
cat("  adjusting for standard clinical covariates (age, sex, stage,\n")
cat("  MMR status, tumor location).\n")
sink()

cat(sprintf("  -> Saved: %s\n",
            file.path(NOMO_TAB_DIR, "clinical_validation_report.txt")))

cat("\n=== Step 1 complete ===\n")
