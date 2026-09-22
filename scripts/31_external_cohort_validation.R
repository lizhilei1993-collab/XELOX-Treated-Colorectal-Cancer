# ============================================================
# 31_EXTERNAL_COHORT_VALIDATION.R
# Step 2: External cohort validation of the 8-variable risk score
#
# Purpose:
#   Validate the pathway-based risk score (8-variable Cox model)
#   in an independent external cohort (GSE72970, FOLFOX-treated).
#
# Key points:
#   - Risk score = sum(pathway_activity * coef) + location_distal * coef
#   - Coefficients from cox_final_results.csv (8-var model)
#   - GSE72970 censoring: 1 = censored, 0 = event (inverted)
#   - PFS as primary endpoint
#
# Input:
#   results/tables/nomogram/cox_final_results.csv
#   results/tables/pathway_activity/GSE72970_pathway_scores.rds
#   results/tables/GSE72970_clinical_data.csv
#
# Output (results/tables/nomogram/):
#   external_validation_report.txt   — full report
#   gse72970_risk_scores.csv         — per-sample risk scores
#   gse72970_km_curve.pdf            — Kaplan-Meier by risk group
#   gse72970_calibration.pdf         — calibration plot (if applicable)
# ============================================================

Sys.setenv(TMPDIR = "C:/temp", TMP = "C:/temp", TEMP = "C:/temp")
.libPaths(c("C:/Rlibs", .libPaths()))

PROJECT_ROOT <- getwd()
RESULTS_TAB_DIR <- file.path(PROJECT_ROOT, "results", "tables")
NOMO_TAB_DIR   <- file.path(RESULTS_TAB_DIR, "nomogram")
NOMO_FIG_DIR   <- file.path(PROJECT_ROOT, "results", "figures", "nomogram")

set.seed(42)

cat("============================================================\n")
cat("Step 2: External Cohort Validation of 8-Variable Risk Score\n")
cat("============================================================\n\n")

# ============================================================
# 0. Load packages
# ============================================================
cat("=== [0] Loading packages ===\n\n")
required_pkgs <- c("survival", "ggplot2", "data.table")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, lib = "C:/Rlibs", repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  cat(sprintf("  %s v%s loaded\n", pkg, packageVersion(pkg)))
}
cat("\n")

# ============================================================
# 1. Load model coefficients
# ============================================================
cat("=== [1] Loading model coefficients ===\n\n")

coef_df <- read.csv(file.path(NOMO_TAB_DIR, "cox_final_results.csv"),
                     stringsAsFactors = FALSE)
coef_vec <- log(coef_df$HR)
names(coef_vec) <- coef_df$Variable
cat("Coefficients (log HR):\n")
for (i in seq_along(coef_vec)) {
  cat(sprintf("  %-35s %.4f\n", names(coef_vec)[i], coef_vec[i]))
}
cat("\n")

# Separate pathway coefs and clinical coefs
pathway_vars <- setdiff(names(coef_vec), "location_distal")
clinical_vars <- intersect(names(coef_vec), "location_distal")
cat(sprintf("  Pathway vars: %d\n", length(pathway_vars)))
cat(sprintf("  Clinical vars: %d\n\n", length(clinical_vars)))

# ============================================================
# 2. Load GSE72970 pathway scores
# ============================================================
cat("=== [2] Loading GSE72970 pathway scores ===\n\n")

pw_scores <- readRDS(file.path(RESULTS_TAB_DIR, "pathway_activity",
                                "GSE72970_pathway_scores.rds"))
cat(sprintf("  Raw dim: %d x %d\n", nrow(pw_scores), ncol(pw_scores)))

# Transpose: samples as rows, pathways as columns
pw_scores_t <- as.data.frame(t(pw_scores))
cat(sprintf("  Transposed dim: %d x %d\n\n", nrow(pw_scores_t), ncol(pw_scores_t)))

# Check which of our 7 pathways are available
available_pw <- intersect(pathway_vars, colnames(pw_scores_t))
missing_pw <- setdiff(pathway_vars, colnames(pw_scores_t))
cat(sprintf("  Available pathways: %d / %d\n", length(available_pw), length(pathway_vars)))
if (length(missing_pw) > 0) {
  cat(sprintf("  Missing: %s\n", paste(missing_pw, collapse = ", ")))
  stop("Missing required pathways in external cohort")
}
cat("\n")

# ============================================================
# 3. Load GSE72970 clinical data
# ============================================================
cat("=== [3] Loading GSE72970 clinical data ===\n\n")

clin_raw <- read.csv(file.path(RESULTS_TAB_DIR, "GSE72970_clinical_data.csv"),
                      stringsAsFactors = FALSE, row.names = 1)
cat(sprintf("  Clinical data: %d samples\n", nrow(clin_raw)))

# Extract relevant columns
clin <- data.frame(
  sample_id = rownames(clin_raw),
  sex       = clin_raw[["Sex.ch1"]],
  age       = as.numeric(clin_raw[["age.ch1"]]),
  regimen   = clin_raw[["regimen.ch1"]],
  location  = clin_raw[["tumor.location.ch1"]],
  response_status = clin_raw[["response.status.ch1"]],
  pfs_time  = as.numeric(clin_raw[["pfs.ch1"]]),
  pfs_censored = as.numeric(clin_raw[["pfs.censored.ch1"]]),
  os_time   = as.numeric(clin_raw[["os.ch1"]]),
  os_censored = as.numeric(clin_raw[["os.censored.ch1"]]),
  pn        = clin_raw[["pn.ch1"]],
  pt        = clin_raw[["pt.ch1"]],
  stringsAsFactors = FALSE
)

# Convert censoring (GSE72970: 1=censored, 0=event)
# Standard survival: 1=event, 0=censored
clin$pfs_event <- ifelse(clin$pfs_censored == 0, 1, 0)
clin$os_event  <- ifelse(clin$os_censored == 0, 1, 0)

cat(sprintf("  PFS events: %d / %d (%.1f%%)\n",
            sum(clin$pfs_event), nrow(clin), mean(clin$pfs_event) * 100))
cat(sprintf("  OS events: %d / %d (%.1f%%)\n\n",
            sum(clin$os_event), nrow(clin), mean(clin$os_event) * 100))

# ============================================================
# 4. Merge pathway scores with clinical data
# ============================================================
cat("=== [4] Merging pathway scores with clinical data ===\n\n")

# Sample IDs - row names in transposed scores are the sample IDs
pw_scores_t$sample_id <- rownames(pw_scores_t)
df <- merge(pw_scores_t[, c("sample_id", available_pw)], clin, by = "sample_id")
cat(sprintf("  Matched samples: %d\n\n", nrow(df)))

# ============================================================
# 5. Compute risk scores
# ============================================================
cat("=== [5] Computing risk scores ===\n\n")

# Pathway component
df$risk_score_pathway <- rep(0, nrow(df))
for (pw in available_pw) {
  df$risk_score_pathway <- df$risk_score_pathway + df[[pw]] * coef_vec[pw]
}

# Location component (distal = 1)
df$location_distal <- ifelse(grepl("Left", df$location, ignore.case = TRUE), 1, 0)
df$risk_score_location <- df$location_distal * ifelse(length(clinical_vars) > 0,
                                                       coef_vec["location_distal"], 0)

# Total risk score
df$risk_score <- df$risk_score_pathway + df$risk_score_location

# Standardized
df$risk_score_z <- scale(df$risk_score)[, 1]

# Risk groups (median-based)
med_rs <- median(df$risk_score)
df$risk_group <- ifelse(df$risk_score >= med_rs, "High risk", "Low risk")
cat(sprintf("  Median risk score: %.4f\n", med_rs))
cat(sprintf("  High risk: %d | Low risk: %d\n\n",
            sum(df$risk_group == "High risk"), sum(df$risk_group == "Low risk")))

# Save risk scores
write.csv(df[, c("sample_id", "risk_score", "risk_score_z", "risk_group",
                  "pfs_event", "pfs_time", "regimen", "location", "response_status")],
          file.path(NOMO_TAB_DIR, "gse72970_risk_scores.csv"), row.names = FALSE)
cat(sprintf("  -> Saved: %s\n\n",
            file.path(NOMO_TAB_DIR, "gse72970_risk_scores.csv")))

# ============================================================
# 6. Overall validation (all regimens)
# ============================================================
cat("=== [6] Overall validation (PFS - all patients) ===\n\n")

cox_all <- coxph(Surv(pfs_time, pfs_event) ~ risk_score_z, data = df)
s_all <- summary(cox_all)

cat(sprintf("  Cox (risk_score_z): HR = %.4f (%.4f - %.4f), p = %.6f\n",
            s_all$conf.int[1, 1], s_all$conf.int[1, 3],
            s_all$conf.int[1, 4], s_all$coefficients[1, "Pr(>|z|)"]))
cat(sprintf("  C-index: %.4f\n\n", concordance(cox_all)$concordance))

# ============================================================
# 7. FOLFOX-only validation (best match to XELOX)
# ============================================================
cat("=== [7] FOLFOX-only validation (PFS) ===\n\n")

df_folfox <- df[df$regimen == "FOLFOX", ]
cat(sprintf("  FOLFOX patients: %d\n", nrow(df_folfox)))

if (nrow(df_folfox) >= 20) {
  cox_folfox <- coxph(Surv(pfs_time, pfs_event) ~ risk_score_z, data = df_folfox)
  s_folfox <- summary(cox_folfox)

  cat(sprintf("  Cox (risk_score_z): HR = %.4f (%.4f - %.4f), p = %.6f\n",
              s_folfox$conf.int[1, 1], s_folfox$conf.int[1, 3],
              s_folfox$conf.int[1, 4], s_folfox$coefficients[1, "Pr(>|z|)"]))
  cat(sprintf("  C-index: %.4f\n\n", concordance(cox_folfox)$concordance))

  # KM curve - FOLFOX only
  df_folfox$risk_group <- factor(df_folfox$risk_group, levels = c("Low risk", "High risk"))
  km_fit <- survfit(Surv(pfs_time, pfs_event) ~ risk_group, data = df_folfox)

  pdf(file.path(NOMO_FIG_DIR, "gse72970_km_curve.pdf"), width = 8, height = 6)
  par(mar = c(5, 5, 4, 2))
  plot(km_fit, col = c("blue", "red"), lwd = 2,
       xlab = "Time (months)", ylab = "Progression-Free Survival",
       main = "GSE72970 FOLFOX: PFS by Risk Group",
       cex.lab = 1.2, cex.main = 1.3)
  legend("topright", legend = c("Low risk", "High risk"),
         col = c("blue", "red"), lwd = 2, cex = 1.1)

  # Add number at risk table
  mtext(sprintf("Log-rank p = %.4f", s_folfox$sctest["pvalue"]), side = 3,
        line = 0, cex = 0.9)
  dev.off()
  cat(sprintf("  -> Saved: %s\n\n",
              file.path(NOMO_FIG_DIR, "gse72970_km_curve.pdf")))
} else {
  cat("  WARNING: Too few FOLFOX patients for survival analysis\n\n")
}

# ============================================================
# 8. Response association (if response status available)
# ============================================================
cat("=== [8] Risk score vs response status ===\n\n")

if ("response_status" %in% colnames(df)) {
  valid_resp <- df[df$response_status %in% c("R", "NR"), ]
  cat(sprintf("  Valid response data: %d samples\n", nrow(valid_resp)))

  if (nrow(valid_resp) >= 10) {
    # R = responder, NR = non-responder
    cat(sprintf("  Responders: %d | Non-responders: %d\n",
                sum(valid_resp$response_status == "R"),
                sum(valid_resp$response_status == "NR")))

    wt <- tryCatch(wilcox.test(risk_score ~ response_status, data = valid_resp),
                   error = function(e) NULL)
    if (!is.null(wt)) {
      cat(sprintf("  Wilcoxon test: p = %.4f\n", wt$p.value))
    }

    # Medians
    resp_med <- tapply(valid_resp$risk_score, valid_resp$response_status, median)
    cat(sprintf("  Median risk score - R: %.4f, NR: %.4f\n",
                resp_med["R"], resp_med["NR"]))
  }
}
cat("\n")

# ============================================================
# 9. Summary report
# ============================================================
cat("=== [9] Writing summary report ===\n\n")

sink(file.path(NOMO_TAB_DIR, "external_validation_report.txt"))

cat("=== External Cohort Validation Report ===\n")
cat("===========================================\n\n")
cat(sprintf("Cohort: GSE72970\n"))
cat(sprintf("Platform: GPL570 (Affymetrix)\n"))
cat(sprintf("Treatment: %s\n\n",
            paste(unique(df$regimen), collapse = ", ")))

cat(sprintf("Total samples with pathway scores: %d\n", nrow(df)))
cat(sprintf("FOLFOX treated: %d\n\n", sum(df$regimen == "FOLFOX")))

cat("--- A. Model Coefficients Applied ---\n")
for (i in seq_along(coef_vec)) {
  cat(sprintf("  %-35s %.4f\n", names(coef_vec)[i], coef_vec[i]))
}
cat("\n")

cat("--- B. Overall Validation (PFS, all regimens) ---\n")
cat(sprintf("  HR per SD: %.4f (%.4f - %.4f), p = %.6f\n",
            s_all$conf.int[1, 1], s_all$conf.int[1, 3],
            s_all$conf.int[1, 4], s_all$coefficients[1, "Pr(>|z|)"]))
cat(sprintf("  C-index: %.4f\n\n", concordance(cox_all)$concordance))

if (exists("s_folfox")) {
  cat("--- C. FOLFOX Subgroup (PFS) ---\n")
  cat(sprintf("  N = %d\n", nrow(df_folfox)))
  cat(sprintf("  HR per SD: %.4f (%.4f - %.4f), p = %.6f\n",
              s_folfox$conf.int[1, 1], s_folfox$conf.int[1, 3],
              s_folfox$conf.int[1, 4], s_folfox$coefficients[1, "Pr(>|z|)"]))
  cat(sprintf("  C-index: %.4f\n", concordance(cox_folfox)$concordance))
  cat(sprintf("  Log-rank p: %.4f\n\n", s_folfox$sctest["pvalue"]))
}

cat("--- D. Response Association ---\n")
if (exists("wt") && !is.null(wt)) {
  cat(sprintf("  Wilcoxon p = %.4f\n", wt$p.value))
  cat(sprintf("  Responder median: %.4f\n", resp_med["R"]))
  cat(sprintf("  Non-responder median: %.4f\n\n", resp_med["NR"]))
} else {
  cat("  Insufficient response data\n\n")
}

cat("--- E. Summary Statistics ---\n")
cat(sprintf("  Risk score range: %.4f - %.4f\n",
            min(df$risk_score), max(df$risk_score)))
cat(sprintf("  Risk score mean (SD): %.4f (%.4f)\n",
            mean(df$risk_score), sd(df$risk_score)))
cat(sprintf("  FOLFOX KM curve: %s\n",
            file.path(NOMO_FIG_DIR, "gse72970_km_curve.pdf")))

cat("\n--- F. Interpretation ---\n")
cat("  The external validation assesses whether the GSE39582-derived\n")
cat("  risk score generalizes to an independent FOLFOX-treated cohort.\n")
cat("  A significant HR confirms generalizability; a non-significant\n")
cat("  result may reflect cohort differences (platform, regimen, sample size).\n")

sink()
cat(sprintf("  -> Saved: %s\n",
            file.path(NOMO_TAB_DIR, "external_validation_report.txt")))
cat("\n=== Step 2 complete ===\n")
