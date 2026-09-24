# ============================================================
# 32_DCA_AND_TDROC.R (v2)
# Decision Curve Analysis + Time-dependent ROC
#
# Uses GSE39582 XELOX-like subcohort (pre-processed groups)
# Fix: uses GSE39582_xelox_groups.csv instead of raw clinical data
# ============================================================

Sys.setenv(TMPDIR = "/tmp", TMP = "/tmp", TEMP = "/tmp")
.libPaths(c("/home/user/AppData/Local/Temp/R_libs", .libPaths()))

PROJECT_ROOT <- Sys.getenv("XELOX_ROOT",
  "/path/to/xelox_project")

RESULTS_TAB_DIR  <- file.path(PROJECT_ROOT, "results", "tables")
RESULTS_FIG_DIR  <- file.path(PROJECT_ROOT, "results", "figures")
DCA_FIG_DIR      <- file.path(RESULTS_FIG_DIR, "dca")
DCA_TAB_DIR      <- file.path(RESULTS_TAB_DIR, "dca")

suppressWarnings(dir.create(DCA_FIG_DIR, showWarnings = FALSE, recursive = TRUE))
suppressWarnings(dir.create(DCA_TAB_DIR, showWarnings = FALSE, recursive = TRUE))

set.seed(42)

cat("========================================\n")
cat("Script 32: DCA + Time-dependent ROC (v2)\n")
cat("========================================\n\n")

# ============================================================
# 0. Load packages
# ============================================================
cat("Loading packages...\n")
required_pkgs <- c("survival", "rms", "timeROC", "ggplot2", "dplyr", "ggsci")
for (pkg in required_pkgs) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  cat(sprintf("  %s v%s\n", pkg, packageVersion(pkg)))
}
cat("\n")

# ============================================================
# 1. Load data
# ============================================================
cat("=== [1] Loading data ===\n\n")

groups <- read.table(
  file.path(RESULTS_TAB_DIR, "GSE39582_xelox_groups.csv"),
  sep = ",", header = TRUE, stringsAsFactors = FALSE,
  quote = "\"", comment.char = ""
)
cat(sprintf("XELOX groups: %d samples\n", nrow(groups)))

ps <- readRDS(file.path(RESULTS_TAB_DIR, "pathway_activity", "GSE39582_pathway_scores.rds"))
cat(sprintf("Pathway scores: %d pathways x %d samples\n", nrow(ps), ncol(ps)))

# ============================================================
# 2. Build analysis dataset
# ============================================================
cat("\n=== [2] Building analysis dataset ===\n\n")

common_ids <- intersect(groups$sample_id, colnames(ps))
cat(sprintf("Matched XELOX samples with pathway scores: %d\n", length(common_ids)))

surv_df <- data.frame(
  rfs_time  = as.numeric(groups$rfs_delay[match(common_ids, groups$sample_id)]),
  rfs_event = as.numeric(groups$rfs_event[match(common_ids, groups$sample_id)]),
  group     = groups$group[match(common_ids, groups$sample_id)],
  row.names = common_ids
)

pathway_mat <- t(as.matrix(ps[, common_ids]))

nomogram_pathways <- c(
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "KEGG_ECM_RECEPTOR_INTERACTION",
  "KEGG_TGF_BETA_SIGNALING_PATHWAY",
  "KEGG_PATHWAYS_IN_CANCER",
  "HALLMARK_MYC_TARGETS_V2",
  "KEGG_COLORECTAL_CANCER"
)
available_pwys <- intersect(nomogram_pathways, colnames(pathway_mat))

for (pwy in available_pwys) {
  surv_df[[pwy]] <- pathway_mat[common_ids, pwy]
}

surv_df$tumor_location <- ifelse(grepl("distal",
  groups$tumor_location[match(common_ids, groups$sample_id)], ignore.case = TRUE), 1, 0)

surv_df <- surv_df[!is.na(surv_df$rfs_time) & !is.na(surv_df$rfs_event), ]

# ============================================================
# 3. Fit model & compute risk score
# ============================================================
cat(sprintf("Analysis sample: n = %d, events = %d\n",
            nrow(surv_df), sum(surv_df$rfs_event)))

formula_str <- paste("Surv(rfs_time, rfs_event) ~",
                     paste(c(available_pwys, "tumor_location"), collapse = " + "))
cox_fit <- coxph(as.formula(formula_str), data = surv_df, x = TRUE)
surv_df$risk_score <- predict(cox_fit, type = "lp")

C_index <- summary(cox_fit)$concordance[1]
cat(sprintf("C-index: %.4f\n", C_index))

# ============================================================
# 4. Time-dependent ROC (1, 3, 5 year)
# ============================================================
cat("\n=== [4] Time-dependent ROC ===\n\n")

time_points <- c(12, 36, 60)
tdroc_results <- list()
tdroc_auc <- data.frame(
  time_month = time_points,
  AUC = NA_real_,
  se = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(time_points)) {
  t <- time_points[i]

  roc_obj <- tryCatch({
    timeROC(
      T = surv_df$rfs_time,
      delta = surv_df$rfs_event,
      marker = surv_df$risk_score,
      cause = 1,
      weighting = "marginal",
      times = t,
      ROC = TRUE,
      iid = TRUE
    )
  }, error = function(e) {
    cat(sprintf("  timeROC failed at t=%d: %s\n", t, e$message))
    return(NULL)
  })

  if (!is.null(roc_obj)) {
    tdroc_results[[as.character(t)]] <- roc_obj
    auc_val <- roc_obj$AUC[2]
    tdroc_auc$AUC[i] <- auc_val
    cat(sprintf("  t=%d months: AUC = %.4f\n", t, auc_val))
  }
}

write.csv(tdroc_auc, file.path(DCA_TAB_DIR, "tdroc_auc.csv"), row.names = FALSE)

# Plot time-dependent ROC curves
pdf(file.path(DCA_FIG_DIR, "tdroc_curve.pdf"), width = 7, height = 6)

plot(NULL, xlim = c(1, 0), ylim = c(0, 1),
     xlab = "1 - Specificity", ylab = "Sensitivity",
     main = "Time-Dependent ROC Curves\n(Pathway Nomogram, GSE39582 XELOX)")
abline(0, 1, lty = 2, col = "gray70")

colors <- c("#E41A1C", "#377EB8", "#4DAF4A")
for (i in seq_along(tdroc_results)) {
  roc_obj <- tdroc_results[[i]]
  t <- time_points[i]
  lines(1 - roc_obj$SP[, 2], roc_obj$Se[, 2],
        col = colors[i], lwd = 2.5)
}

legend("bottomright",
       legend = sprintf("%d-Month (AUC=%.3f)", time_points, tdroc_auc$AUC),
       col = colors[seq_along(tdroc_results)],
       lwd = 2.5, bty = "n")
dev.off()
cat("Time-dependent ROC plot saved.\n")

# ============================================================
# 5. Decision Curve Analysis (DCA)
# ============================================================
cat("\n=== [5] Decision Curve Analysis ===\n\n")

thresholds <- seq(0.01, 0.50, by = 0.01)
dca_results <- data.frame(
  threshold = thresholds,
  nb_treat_all = NA_real_,
  nb_treat_none = 0,
  nb_model = NA_real_,
  stringsAsFactors = FALSE
)

# Predict 1-year survival probability
cox_surv <- survfit(cox_fit, newdata = surv_df)
# Get 12-month survival for each patient
surv_12m <- summary(cox_surv, times = 12)
# Match survival probabilities to patients
# Use linear interpolation for safety
if (length(surv_12m$surv) == nrow(surv_df)) {
  surv_12m_prob <- surv_12m$surv
} else {
  # Fall back to predictSurvProb
  surv_12m_prob <- as.numeric(summary(survfit(cox_fit, newdata = surv_df), times = 12)$surv)
}

for (i in seq_along(thresholds)) {
  pt <- thresholds[i]

  # High-risk = predicted 12-month survival < (1 - pt)
  high_risk <- surv_12m_prob < (1 - pt)

  # Actual event within 12 months
  actual_event <- surv_df$rfs_event == 1 & surv_df$rfs_time <= 12

  # Net benefit
  TP <- sum(high_risk & actual_event, na.rm = TRUE)
  FP <- sum(high_risk & !actual_event, na.rm = TRUE)
  n  <- nrow(surv_df)

  nb_model <- (TP / n) - (FP / n) * (pt / (1 - pt))

  # Treat all
  nb_all <- mean(actual_event, na.rm = TRUE) -
            (1 - mean(actual_event, na.rm = TRUE)) * (pt / (1 - pt))

  dca_results$nb_model[i] <- nb_model
  dca_results$nb_treat_all[i] <- nb_all
}

# Write net benefit table
write.csv(dca_results, file.path(DCA_TAB_DIR, "dca_net_benefit.csv"), row.names = FALSE)

# Plot DCA
pdf(file.path(DCA_FIG_DIR, "dca_curve.pdf"), width = 8, height = 6)

plot(dca_results$threshold, dca_results$nb_treat_all,
     type = "l", lty = 2, col = "gray50", lwd = 2,
     xlab = "Threshold Probability",
     ylab = "Net Benefit",
     main = "Decision Curve Analysis\n(Pathway Nomogram, GSE39582 XELOX)",
     ylim = range(c(dca_results$nb_model, dca_results$nb_treat_all, 0), na.rm = TRUE),
     xlim = c(0, 0.5))

lines(dca_results$threshold, dca_results$nb_model,
      col = "#E41A1C", lwd = 3)

lines(dca_results$threshold, rep(0, nrow(dca_results)),
      lty = 2, col = "gray30", lwd = 1)

legend("topright",
       legend = c("Pathway Nomogram", "Treat All", "Treat None"),
       col = c("#E41A1C", "gray50", "gray30"),
       lty = c(1, 2, 2),
       lwd = c(3, 2, 1), bty = "n")
dev.off()
cat("DCA plot saved.\n")

# Summary
cat("\n========================================\n")
cat("Script 32: DCA + TDROC COMPLETE\n")
cat("========================================\n")
cat("\nKey Results:\n")
cat(sprintf("  C-index: %.4f\n", C_index))
cat(sprintf("  TDROC AUC (12M): %.3f\n", tdroc_auc$AUC[1]))
cat(sprintf("  TDROC AUC (36M): %.3f\n", tdroc_auc$AUC[2]))
cat(sprintf("  TDROC AUC (60M): %.3f\n", tdroc_auc$AUC[3]))
cat(sprintf("\nOutput: %s\n", DCA_FIG_DIR))
cat(sprintf("Output: %s\n", DCA_TAB_DIR))
