# ============================================================
# 31_CALIBRATION_ANALYSIS.R (v2)
# Calibration Curve + Calibration Statistics + Brier Score
#
# Uses GSE39582 XELOX-like subcohort (pre-processed groups)
# Fix: uses GSE39582_xelox_groups.csv instead of raw clinical data
# ============================================================

Sys.setenv(TMPDIR = "C:/temp", TMP = "C:/temp", TEMP = "C:/temp")
.libPaths(c("/home/user/AppData/Local/Temp/R_libs", .libPaths()))

PROJECT_ROOT <- Sys.getenv("XELOX_ROOT",
  "/path/to/xelox_project")

RESULTS_TAB_DIR  <- file.path(PROJECT_ROOT, "results", "tables")
RESULTS_FIG_DIR  <- file.path(PROJECT_ROOT, "results", "figures")
CALIB_FIG_DIR    <- file.path(RESULTS_FIG_DIR, "calibration")
CALIB_TAB_DIR    <- file.path(RESULTS_TAB_DIR, "calibration")

suppressWarnings(dir.create(CALIB_FIG_DIR, showWarnings = FALSE, recursive = TRUE))
suppressWarnings(dir.create(CALIB_TAB_DIR, showWarnings = FALSE, recursive = TRUE))

set.seed(42)

cat("========================================\n")
cat("Script 31: Calibration Analysis (v2)\n")
cat("========================================\n\n")

# ============================================================
# 0. Load packages
# ============================================================
cat("Loading packages...\n")
required_pkgs <- c("survival", "rms", "pec", "ggplot2", "riskRegression")
for (pkg in required_pkgs) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  cat(sprintf("  %s v%s\n", pkg, packageVersion(pkg)))
}
cat("\n")

# ============================================================
# 1. Load data — use pre-processed XELOX groups and pathway scores
# ============================================================
cat("=== [1] Loading data ===\n\n")

# Read XELOX groups (using read.table with explicit quote handling)
groups <- read.table(
  file.path(RESULTS_TAB_DIR, "GSE39582_xelox_groups.csv"),
  sep = ",", header = TRUE, stringsAsFactors = FALSE,
  quote = "\"", comment.char = ""
)
cat(sprintf("XELOX groups: %d samples\n", nrow(groups)))
cat(sprintf("  Sensitive: %d | Resistant: %d | Intermediate: %d\n",
            sum(groups$group == "sensitive"),
            sum(groups$group == "resistant"),
            sum(groups$group == "intermediate")))

# Pathway scores (44 pathways x 585 samples)
ps <- readRDS(file.path(RESULTS_TAB_DIR, "pathway_activity", "GSE39582_pathway_scores.rds"))
cat(sprintf("Pathway scores: %d pathways x %d samples\n", nrow(ps), ncol(ps)))

# ============================================================
# 2. Build analysis dataset
# ============================================================
cat("\n=== [2] Building analysis dataset ===\n\n")

# Match samples
common_ids <- intersect(groups$sample_id, colnames(ps))
cat(sprintf("Matched XELOX samples with pathway scores: %d\n", length(common_ids)))

# Create survival data frame
surv_df <- data.frame(
  rfs_time  = as.numeric(groups$rfs_delay[match(common_ids, groups$sample_id)]),
  rfs_event = as.numeric(groups$rfs_event[match(common_ids, groups$sample_id)]),
  group     = groups$group[match(common_ids, groups$sample_id)],
  row.names = common_ids
)

# Add pathway scores (transpose: samples as rows, pathways as cols)
pathway_mat <- t(as.matrix(ps[, common_ids]))

# Define the 7 nomogram pathways
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
cat(sprintf("Available nomogram pathways: %d / %d\n",
            length(available_pwys), length(nomogram_pathways)))

for (pwy in available_pwys) {
  surv_df[[pwy]] <- pathway_mat[common_ids, pwy]
}

# Add tumor location (from groups)
surv_df$tumor_location <- ifelse(grepl("distal",
  groups$tumor_location[match(common_ids, groups$sample_id)], ignore.case = TRUE), 1, 0)

# Remove samples with missing RFS
surv_df <- surv_df[!is.na(surv_df$rfs_time) & !is.na(surv_df$rfs_event), ]
cat(sprintf("Final analysis sample: n = %d, events = %d\n",
            nrow(surv_df), sum(surv_df$rfs_event)))

# ============================================================
# 3. Fit Cox model
# ============================================================
cat("\n=== [3] Fitting Cox model ===\n\n")

formula_str <- paste("Surv(rfs_time, rfs_event) ~",
                     paste(c(available_pwys, "tumor_location"), collapse = " + "))
cox_fit <- coxph(as.formula(formula_str), data = surv_df, x = TRUE, y = TRUE)
print(summary(cox_fit))

# ============================================================
# 4. Calibration curve (rms-based)
# ============================================================
cat("\n=== [4] Calibration curve ===\n\n")

dd <- datadist(surv_df)
options(datadist = "dd")

cph_formula <- as.formula(paste("Surv(rfs_time, rfs_event) ~",
                                 paste(c(available_pwys, "tumor_location"), collapse = " + ")))
cph_fit <- cph(cph_formula, data = surv_df, x = TRUE, y = TRUE, surv = TRUE)

# Time points
time_points <- c(12, 36, 60)

# Calibration plot
pdf(file.path(CALIB_FIG_DIR, "calibration_curve.pdf"), width = 10, height = 4)
par(mfrow = c(1, 3))
for (t in time_points) {
  cal <- calibrate(cph_fit, cmethod = "KM", method = "boot",
                   u = t, B = 200, m = min(50, round(nrow(surv_df) / 3)))
  plot(cal, xlab = sprintf("Predicted %d-Month Survival", t),
       ylab = sprintf("Observed %d-Month Survival", t),
       main = sprintf("Calibration at %d Months", t),
       subtitles = FALSE)
  abline(0, 1, col = "gray", lty = 2)
}
dev.off()
cat("Calibration curve saved.\n")

# ============================================================
# 5. Calibration statistics (using rms::val.surv)
# ============================================================
cat("\n=== [5] Calibration statistics ===\n\n")

cal_stats <- data.frame(
  time_point = time_points,
  calibration_slope = NA_real_,
  calibration_intercept = NA_real_,
  brier_score = NA_real_,
  stringsAsFactors = FALSE
)

# Use rms's validate.cph to get calibration slope and Dxy
# Then compute Brier score via pec
set.seed(42)
val <- validate(cph_fit, B = 100)
cat("Model validation summary (B=100):\n")
print(val)

# Extract calibration slope from validation
cal_stats$calibration_slope[1] <- val["Slope", "index.corrected"]
cat(sprintf("  Calibration slope (optimism-corrected): %.4f\n", cal_stats$calibration_slope[1]))

# Also compute Dxy (discrimination)
cat(sprintf("  Dxy (training): %.4f\n", val["Dxy", "training"]))
cat(sprintf("  Dxy (optimism-corrected): %.4f\n", val["Dxy", "index.corrected"]))

# Write calibration stats (partial)
write.csv(cal_stats, file.path(CALIB_TAB_DIR, "calibration_stats.csv"), row.names = FALSE)
cat("Partial calibration statistics saved.\n")

# ============================================================
# 6. Time-dependent Brier Score (via pec)
# ============================================================
cat("\n=== [6] Brier Score ===\n\n")

tryCatch({
  # Simpler approach: use pec with coxph and Kaplan-Meier reference
  brier_scores <- pec(
    object = list("Pathway Nomogram" = cox_fit),
    formula = Surv(rfs_time, rfs_event) ~ 1,
    data = surv_df,
    times = time_points,
    start = 1,
    maxtime = max(time_points),
    exact = FALSE,
    cens.model = "marginal",
    verbose = FALSE
  )

  cat("pec object structure:\n")
  print(names(brier_scores))

  # Extract Brier scores - check structure
  if (!is.null(brier_scores$AppErr)) {
    cat("AppErr names:", names(brier_scores$AppErr), "\n")
    cat("Brier for nomogram:", brier_scores$AppErr[["Pathway Nomogram"]], "\n")
    cal_stats$brier_score <- as.numeric(brier_scores$AppErr[["Pathway Nomogram"]])

    # Reference (null model)
    if (!is.null(brier_scores$AppErr$Reference)) {
      cal_stats$brier_score_null <- as.numeric(brier_scores$AppErr$Reference)
      cat("Brier for Reference:", brier_scores$AppErr$Reference, "\n")
    }
  }

  # Update stats file
  write.csv(cal_stats, file.path(CALIB_TAB_DIR, "calibration_stats.csv"), row.names = FALSE)

  cat("\nBrier scores:\n")
  for (i in seq_along(time_points)) {
    null_str <- ifelse(is.null(cal_stats$brier_score_null), "N/A",
                       sprintf("%.4f", cal_stats$brier_score_null[i]))
    cat(sprintf("  t=%dM: Brier=%.4f (null=%s)\n",
                time_points[i], cal_stats$brier_score[i], null_str))
  }
}, error = function(e) {
  cat("Brier score computation failed:", e$message, "\n")
  cat("Saving partial results.\n")
  write.csv(cal_stats, file.path(CALIB_TAB_DIR, "calibration_stats.csv"), row.names = FALSE)
})

cat("\n========================================\n")
cat("Script 31: Calibration Analysis COMPLETE\n")
cat("========================================\n")
cat(sprintf("Output: %s\n", CALIB_FIG_DIR))
cat(sprintf("Output: %s\n", CALIB_TAB_DIR))
