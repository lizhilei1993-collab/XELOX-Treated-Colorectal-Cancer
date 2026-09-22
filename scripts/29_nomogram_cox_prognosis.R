# ============================================================
# 29_NOMOGRAM_COX_PROGNOSIS.R
# Phase IV: Pathway-based Cox Prognostic Model & Nomogram
#
# Purpose:
#   Build a survival-based risk stratification model for
#   XELOX-treated mCRC patients using GSE39582 RFS data.
#   Instead of binary resistance prediction (which failed),
#   this uses time-to-event analysis to identify pathways
#   associated with relapse-free survival.
#
# Key differences from previous PRS approach:
#   - Endpoint: RFS (time-to-event) vs binary Resistant/Sensitive
#   - Feature selection: LASSO-Cox vs LASSO-Logistic
#   - Target: Risk stratification (prognostic) vs drug response (predictive)
#   - Samples: ALL 227 with valid RFS (intermediate included)
#
# Input:
#   results/tables/GSE39582_xelox_groups.csv     — clinical + survival
#   results/tables/pathway_activity/GSE39582_pathway_scores.rds  — 44 pathways
#
# Output (results/tables/nomogram/ + results/figures/nomogram/):
#   - Univariate Cox screening results
#   - Multivariable Cox (pathway + clinical) results
#   - Clinical-only vs pathway-model comparison
#   - Nomogram (12/36/60-month RFS)
#   - Calibration plots
#   - KM curves (high/low risk)
#   - Time-dependent ROC
#   - Bootstrap internal validation
#   - Analysis report
#
# Note: LASSO-Cox selected 0 pathways (signal too weak for L1 penalty).
# Switched to: univariate p<0.05 filter + backward AIC selection.
# This is standard for moderate-N clinical prediction modeling.
# ============================================================

Sys.setenv(TMPDIR = "C:/temp", TMP = "C:/temp", TEMP = "C:/temp")
.libPaths(c("C:/Rlibs", .libPaths()))

PROJECT_ROOT <- "/path/to/xelox_project"

RESULTS_TAB_DIR <- file.path(PROJECT_ROOT, "results", "tables")
RESULTS_FIG_DIR <- file.path(PROJECT_ROOT, "results", "figures")
PA_DIR         <- file.path(RESULTS_TAB_DIR, "pathway_activity")
NOMO_TAB_DIR   <- file.path(RESULTS_TAB_DIR, "nomogram")
NOMO_FIG_DIR   <- file.path(RESULTS_FIG_DIR, "nomogram")

suppressWarnings(dir.create(NOMO_TAB_DIR, showWarnings = FALSE, recursive = TRUE))
suppressWarnings(dir.create(NOMO_FIG_DIR, showWarnings = FALSE, recursive = TRUE))

set.seed(42)

cat("============================================================\n")
cat("Phase IV: Pathway-based Cox Prognostic Model & Nomogram\n")
cat("Cohort: GSE39582 XELOX subgroup (RFS endpoint)\n")
cat("============================================================\n\n")

# ============================================================
# 0. Load packages
# ============================================================
cat("=== [0] Loading packages ===\n\n")
required_pkgs <- c("glmnet", "survival", "survminer", "rms",
                    "ggplot2", "data.table", "timeROC")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, lib = "C:/Rlibs", repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  cat(sprintf("  %s v%s loaded\n", pkg, packageVersion(pkg)))
}
cat("\n")

# ============================================================
# 1. Load clinical + pathway data
# ============================================================
cat("=== [1] Loading data ===\n\n")

# --- 1a. Clinical ---
clin <- read.csv(file.path(RESULTS_TAB_DIR, "GSE39582_xelox_groups.csv"),
                  stringsAsFactors = FALSE)
cat(sprintf("  Clinical data: %d samples loaded\n", nrow(clin)))

# --- 1b. Pathway scores ---
ps <- readRDS(file.path(PA_DIR, "GSE39582_pathway_scores.rds"))
cat(sprintf("  Pathway scores: %d pathways x %d samples\n",
            nrow(ps), ncol(ps)))

# --- 1c. Merge ---
# Pathway scores are genes x samples; transpose to samples x pathways
ps_t <- t(ps)
ps_df <- as.data.frame(ps_t)
ps_df$sample_id <- rownames(ps_t)

cox_df <- merge(clin, ps_df, by = "sample_id", all.x = TRUE)
cat(sprintf("  Merged data: %d samples\n", nrow(cox_df)))

# ============================================================
# 2. Prepare survival data
# ============================================================
cat("\n=== [2] Preparing survival data ===\n\n")

# Filter for valid RFS
cox_df <- cox_df[!is.na(cox_df$rfs_event) & !is.na(cox_df$rfs_delay) &
                   cox_df$rfs_delay > 0, ]
cox_df$rfs_event <- as.numeric(cox_df$rfs_event)
cox_df$rfs_delay <- as.numeric(cox_df$rfs_delay)

cat(sprintf("  Samples with valid RFS: %d\n", nrow(cox_df)))
cat(sprintf("  RFS events: %d / %d (%.1f%%)\n",
            sum(cox_df$rfs_event), nrow(cox_df),
            mean(cox_df$rfs_event) * 100))
cat(sprintf("  Median follow-up: %.1f months (range: %.0f-%.0f)\n",
            median(cox_df$rfs_delay), min(cox_df$rfs_delay), max(cox_df$rfs_delay)))

# --- 2b. Extract pathway matrix ---
pathway_names <- rownames(ps)
pw_matrix <- as.matrix(cox_df[, pathway_names, drop = FALSE])
cat(sprintf("  Pathway feature matrix: %d samples x %d pathways\n",
            nrow(pw_matrix), ncol(pw_matrix)))

# Check for missing values
if (any(is.na(pw_matrix))) {
  cat("  WARNING: Missing values in pathway scores - imputing with column medians\n")
  for (j in seq_len(ncol(pw_matrix))) {
    pw_matrix[is.na(pw_matrix[, j]), j] <- median(pw_matrix[, j], na.rm = TRUE)
  }
}

# Standardize pathways for LASSO
pw_scaled <- scale(pw_matrix)

# ============================================================
# 3. Step 1: Univariate Cox screening
# ============================================================
cat("\n=== [3] Univariate Cox screening ===\n\n")

uni_results <- data.frame(
  Pathway = pathway_names,
  HR = NA, HR_lower = NA, HR_upper = NA,
  p_value = NA, stringsAsFactors = FALSE
)

for (i in seq_along(pathway_names)) {
  pw <- pathway_names[i]
  tryCatch({
    fit <- coxph(Surv(rfs_delay, rfs_event) ~ pw_scaled[, i], data = cox_df)
    s <- summary(fit)
    uni_results$HR[i] <- exp(coef(fit))
    uni_results$HR_lower[i] <- exp(confint(fit)[1])
    uni_results$HR_upper[i] <- exp(confint(fit)[2])
    uni_results$p_value[i] <- s$coefficients[1, "Pr(>|z|)"]
  }, error = function(e) {
    uni_results$p_value[i] <- NA
  })
}

# Sort by p-value
uni_results <- uni_results[order(uni_results$p_value), ]

# Print top pathways
cat("  Top 10 pathways by univariate Cox p-value:\n")
cat(sprintf("  %-50s %8s %8s %8s %8s\n", "Pathway", "HR", "Lower", "Upper", "p"))
for (i in 1:min(10, nrow(uni_results))) {
  r <- uni_results[i, ]
  cat(sprintf("  %-50s %7.3f %7.3f %7.3f %7.4f\n",
              substr(r$Pathway, 1, 50), r$HR, r$HR_lower, r$HR_upper, r$p_value))
}

# NOTE: LASSO-Cox was attempted but selected 0 pathways at all lambda
# values, indicating the per-pathway survival signal is too weak for L1
# penalization. Switching to standard clinical modeling workflow:
#   p<0.05 univariate filter → backward AIC → final model

# Use p<0.05 as selection threshold (standard in clinical Cox modeling)
candidate_mask <- which(uni_results$p_value < 0.05 & !is.na(uni_results$p_value))
candidate_pathways <- uni_results$Pathway[candidate_mask]
cat(sprintf("\n  Pathways passing univariate filter (p<0.05): %d / %d\n",
            length(candidate_pathways), length(pathway_names)))

# List candidates
cat("  Candidate pathways:\n")
for (pw in candidate_pathways) {
  r <- uni_results[uni_results$Pathway == pw, ]
  cat(sprintf("    %s  (HR=%.3f, p=%.4f)\n", pw, r$HR, r$p_value))
}

if (length(candidate_pathways) < 2) {
  cat("  WARNING: Too few candidates. Relaxing to p<0.10...\n")
  candidate_mask <- which(uni_results$p_value < 0.1 & !is.na(uni_results$p_value))
  candidate_pathways <- uni_results$Pathway[candidate_mask]
  cat(sprintf("  Pathways passing univariate filter (p<0.10): %d / %d\n",
              length(candidate_pathways), length(pathway_names)))
}

# Save univariate results
write.csv(uni_results, file.path(NOMO_TAB_DIR, "univariate_cox_screening.csv"),
          row.names = FALSE)
cat(sprintf("  -> Saved: %s\n",
            file.path(NOMO_TAB_DIR, "univariate_cox_screening.csv")))

# ============================================================
# 4. Build risk score from top univariate pathways
# ============================================================
cat("\n=== [4] Risk score computation ===\n\n")

# Standardize pathway scores
for (pw in candidate_pathways) {
  cname <- make.names(pw)
  cox_df[[cname]] <- pw_scaled[, which(pathway_names == pw)]
}

# --- 4a. Pathway-only Cox with backward AIC ---
pw_formula_vars <- paste0("`", make.names(candidate_pathways), "`", collapse = " + ")
pw_formula <- as.formula(paste0("Surv(rfs_delay, rfs_event) ~ ", pw_formula_vars))

cox_pw <- coxph(pw_formula, data = cox_df)
cox_pw_step <- step(cox_pw, direction = "backward", trace = 0)
pw_selected <- setdiff(names(coef(cox_pw_step)), "Intercept")

cat("  Pathway-only model (backward AIC):\n")
cat(sprintf("    Pathways retained: %d\n", length(pw_selected)))
for (pw in pw_selected) {
  cat(sprintf("    %s\n", pw))
}

# --- 4b. Clinical-only model ---
cox_df$sex_bin <- ifelse(cox_df$sex == "Male", 1, 0)
cox_df$age_num <- as.numeric(cox_df$age)
cox_df$stage_III_IV <- ifelse(grepl("3|4|III|IV", cox_df$tnm_stage), 1, 0)
cox_df$mmr_dmmr <- ifelse(cox_df$mmr_status == "dMMR", 1, 0)
cox_df$location_distal <- ifelse(cox_df$tumor_location == "distal", 1, 0)

clin_formula <- Surv(rfs_delay, rfs_event) ~ age_num + sex_bin + stage_III_IV +
  mmr_dmmr + location_distal
cox_clin <- coxph(clin_formula, data = cox_df)
cox_clin_step <- step(cox_clin, direction = "backward", trace = 0)
clin_selected <- setdiff(names(coef(cox_clin_step)), "Intercept")

cat("\n  Clinical-only model (backward AIC):\n")
cat(sprintf("    Variables retained: %d\n", length(clin_selected)))
for (v in clin_selected) {
  cat(sprintf("    %s\n", v))
}

# --- 4c. Combined model: pathway + clinical ---
if (length(pw_selected) > 0 && length(clin_selected) > 0) {
  combined_vars <- unique(c(pw_selected, clin_selected))
  combined_formula <- as.formula(
    paste0("Surv(rfs_delay, rfs_event) ~ ",
           paste0("`", combined_vars, "`", collapse = " + "))
  )
  cox_combined <- coxph(combined_formula, data = cox_df)
  cox_final <- step(cox_combined, direction = "backward", trace = 0)
} else if (length(pw_selected) > 0) {
  cox_final <- cox_pw_step
} else {
  cox_final <- cox_clin_step
}

final_vars <- setdiff(names(coef(cox_final)), "Intercept")
final_formula <- formula(cox_final)
cat("\n  Final combined model (backward AIC):\n")
cat(sprintf("    Variables retained: %d\n", length(final_vars)))
for (v in final_vars) {
  cat(sprintf("    %s\n", v))
}

# --- 4d. Compute risk score ---
if (length(final_vars) > 0) {
  # Linear predictor from final model
  risk_score <- predict(cox_final, type = "lp")
  cox_df$risk_score <- risk_score
  cat(sprintf("\n  Risk score range: %.3f to %.3f\n",
              min(risk_score), max(risk_score)))
} else {
  cat("\n  WARNING: No variables retained. Setting risk score to 0.\n")
  cox_df$risk_score <- 0
}

# Stratify by median
cox_df$risk_group <- ifelse(cox_df$risk_score > median(cox_df$risk_score),
                             "High risk", "Low risk")
cox_df$risk_group <- factor(cox_df$risk_group,
                             levels = c("Low risk", "High risk"))
cat(sprintf("  Low risk: %d, High risk: %d\n",
            sum(cox_df$risk_group == "Low risk"),
            sum(cox_df$risk_group == "High risk")))

# --- 4e. Compare C-index ---
c_index_clin <- concordance(cox_clin_step)$concordance
c_index_pw <- if (length(pw_selected) > 0) concordance(cox_pw_step)$concordance else NA
c_index_full <- concordance(cox_final)$concordance

cat(sprintf("\n  C-index comparison:\n"))
cat(sprintf("    Clinical-only: %.3f\n", c_index_clin))
if (!is.na(c_index_pw)) cat(sprintf("    Pathway-only: %.3f\n", c_index_pw))
cat(sprintf("    Combined:      %.3f\n", c_index_full))

# ============================================================
# 5. Multivariable Cox results & forest plot
# ============================================================
cat("\n=== [5] Multivariable Cox Results ===\n\n")

# Save final model summary
sink(file.path(NOMO_TAB_DIR, "cox_final_model.txt"))
cat("=== Final Cox Regression Model ===\n\n")
cat(sprintf("Cohort: GSE39582 XELOX subgroup\n"))
cat(sprintf("Samples: %d, Events: %d (%.1f%%)\n\n",
            nrow(cox_df), sum(cox_df$rfs_event), mean(cox_df$rfs_event)*100))
cat("--- Model selection ---\n")
cat(sprintf("Pathway candidates: %d (p<0.05 univariate filter)\n", length(candidate_pathways)))
cat(sprintf("Clinical candidates: age, sex, stage, MMR, location\n"))
cat("Method: backward AIC (pathway-only → clinical-only → combined)\n\n")
cat("--- Final model ---\n")
cat(deparse(final_formula, width.cutoff = 80))
cat("\n\n")
cat("--- Final model summary ---\n")
print(summary(cox_final))
sink()

# Extract results table
coef_final <- summary(cox_final)$coefficients
ci_final <- summary(cox_final)$conf.int
multi_results <- data.frame(
  Variable = names(coef(cox_final)),
  HR = round(exp(coef(cox_final)), 3),
  HR_lower = round(exp(confint(cox_final)[, 1]), 3),
  HR_upper = round(exp(confint(cox_final)[, 2]), 3),
  p_value = round(coef_final[, "Pr(>|z|)"], 4),
  stringsAsFactors = FALSE
)
write.csv(multi_results, file.path(NOMO_TAB_DIR, "cox_final_results.csv"),
          row.names = FALSE)
cat(sprintf("  -> Saved: %s\n",
            file.path(NOMO_TAB_DIR, "cox_final_results.csv")))

# C-index from final model
c_index_val <- concordance(cox_final)$concordance
cat(sprintf("  Clinical-only C-index: %.3f\n", c_index_clin))
cat(sprintf("  Final model C-index:   %.3f\n\n", c_index_full))

# --- Forest plot ---
if (nrow(multi_results) > 0) {
  nv <- nrow(multi_results)
  pdf(file.path(NOMO_FIG_DIR, "cox_forest_plot.pdf"),
      width = 9, height = max(4, nv * 0.5 + 2))
  par(mar = c(4, 10, 3, 4))
  hr_max <- max(multi_results$HR_upper, na.rm = TRUE) * 1.4
  plot(NA, xlim = c(0, max(hr_max, 2)),
       ylim = c(0.5, nv + 1),
       xlab = "Hazard Ratio (95% CI)", ylab = "", yaxt = "n",
       main = "Final Cox Model: Pathway + Clinical Factors")
  axis(2, at = nv:1, labels = rev(multi_results$Variable), las = 2, cex.axis = 0.8)
  abline(v = 1, lty = 2, col = "gray")
  for (i in seq_len(nv)) {
    idx <- nv - i + 1
    col_use <- ifelse(multi_results$p_value[idx] < 0.05, "darkred", "gray40")
    points(multi_results$HR[idx], i, pch = 15, cex = 1.2, col = col_use)
    segments(multi_results$HR_lower[idx], i, multi_results$HR_upper[idx], i,
             lwd = 2, col = col_use)
    text(hr_max * 0.95, i,
         sprintf("%.2f (%.2f-%.2f)", multi_results$HR[idx],
                 multi_results$HR_lower[idx], multi_results$HR_upper[idx]),
         cex = 0.7, adj = 0)
  }
  dev.off()
  cat(sprintf("  -> Saved: %s\n",
              file.path(NOMO_FIG_DIR, "cox_forest_plot.pdf")))
}

# ============================================================
# 6. Kaplan-Meier curves (risk-stratified)
# ============================================================
cat("\n=== [6] KM Curves ===\n\n")

km_fit <- survfit(Surv(rfs_delay, rfs_event) ~ risk_group, data = cox_df)

pdf(file.path(NOMO_FIG_DIR, "km_risk_groups.pdf"), width = 8, height = 7)
survminer::ggsurvplot(
  km_fit,
  data = cox_df,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  xlab = "Relapse-Free Survival (months)",
  ylab = "RFS Probability",
  title = "Risk-stratified RFS in GSE39582 XELOX Cohort",
  palette = c("steelblue", "darkred"),
  legend.title = "Risk Group",
  ggtheme = theme_bw()
)
dev.off()
cat(sprintf("  -> Saved: %s\n",
            file.path(NOMO_FIG_DIR, "km_risk_groups.pdf")))

# Log-rank
lr_test <- survdiff(Surv(rfs_delay, rfs_event) ~ risk_group, data = cox_df)
lr_p <- 1 - pchisq(lr_test$chisq, df = 1)
cat(sprintf("  Log-rank p = %.4f\n", lr_p))

# ============================================================
# 7. Time-dependent ROC
# ============================================================
cat("\n=== [7] Time-dependent ROC ===\n\n")

surv_times <- c(12, 36, 60)

tryCatch({
  roc_td <- timeROC(
    T = cox_df$rfs_delay,
    delta = cox_df$rfs_event,
    marker = cox_df$risk_score,
    cause = 1,
    times = surv_times,
    iid = TRUE
  )

  cat("  Time-dependent AUC:\n")
  for (i in seq_along(surv_times)) {
    cat(sprintf("    %d-month AUC = %.3f (SE = %.4f)\n",
                surv_times[i], roc_td$AUC[i], roc_td$inference$vect_sd_1[i]))
  }

  pdf(file.path(NOMO_FIG_DIR, "time_dependent_roc.pdf"), width = 8, height = 8)
  plot(roc_td, time = surv_times[1], col = "orange", lwd = 2,
       title = "Time-dependent ROC for Risk Score")
  for (i in 2:length(surv_times)) {
    plot(roc_td, time = surv_times[i], col = c("forestgreen", "darkred")[i - 1],
         lwd = 2, add = TRUE)
  }
  legend("bottomright",
         legend = sprintf("%d-month AUC = %.3f", surv_times, roc_td$AUC),
         col = c("orange", "forestgreen", "darkred"),
         lwd = 2, cex = 0.9)
  dev.off()
  cat(sprintf("  -> Saved: %s\n",
              file.path(NOMO_FIG_DIR, "time_dependent_roc.pdf")))
}, error = function(e) {
  cat(sprintf("  WARNING: timeROC failed: %s\n", e$message))
})

# ============================================================
# 8. Nomogram
# ============================================================
cat("\n=== [8] Nomogram ===\n\n")

if (length(final_vars) >= 2) {
  tryCatch({
    # Set time unit to Months (fixes rms "30 Days" display bug)
    cox_df$rfs_delay <- as.numeric(cox_df$rfs_delay)
    units(cox_df$rfs_delay) <- "Month"

    # Prepare variables for datadist
    nomo_vars <- c(final_vars, "rfs_delay", "rfs_event")
    dd <- datadist(cox_df[, nomo_vars])
    options(datadist = "dd")

    # Build cph formula from final model
    cph_formula <- as.formula(
      paste0("Surv(rfs_delay, rfs_event) ~ ",
             paste0("`", final_vars, "`", collapse = " + "))
    )
    cox_cph <- cph(cph_formula, data = cox_df, x = TRUE, y = TRUE, surv = TRUE)

    # Survival functions
    surv_fn_12 <- function(lp) Survival(cox_cph)(times = 12, lp = lp)
    surv_fn_36 <- function(lp) Survival(cox_cph)(times = 36, lp = lp)
    surv_fn_60 <- function(lp) Survival(cox_cph)(times = 60, lp = lp)

    # Nomogram
    pdf(file.path(NOMO_FIG_DIR, "nomogram.pdf"), width = 14, height = 9)
    nom <- nomogram(cox_cph,
                    fun = list(surv_fn_12, surv_fn_36, surv_fn_60),
                    funlabel = c("12-month RFS", "36-month RFS", "60-month RFS"),
                    maxscale = 100, lp = TRUE)
    plot(nom, xfrac = 0.35, cex.var = 0.7, cex.axis = 0.65)
    dev.off()
    cat(sprintf("  -> Saved: %s\n",
                file.path(NOMO_FIG_DIR, "nomogram.pdf")))

    # ============================================================
    # 9. Calibration plots + quantitative metrics
    # ============================================================
    cat("\n=== [9] Calibration ===\n\n")

    cal_metrics <- data.frame(
      TimePoint = surv_times,
      CalibrationSlope = NA,
      CalibrationIntercept = NA,
      BrierScore = NA,
      stringsAsFactors = FALSE
    )

    for (i in seq_along(surv_times)) {
      st <- surv_times[i]
      tryCatch({
        # --- 9a. rms::calibrate() for slope & intercept ---
        cal_file <- file.path(NOMO_FIG_DIR, sprintf("calibration_%dmo.pdf", st))
        pdf(cal_file, width = 7, height = 7)
        cal <- calibrate(cox_cph, u = st, B = 200)
        plot(cal, subtitles = FALSE,
             xlab = paste0("Predicted ", st, "-month RFS"),
             ylab = paste0("Observed ", st, "-month RFS"),
             main = paste0("Bootstrap Calibration (", st, " months, B=200)"))
        dev.off()
        cat(sprintf("  -> Saved: %s\n", cal_file))

        # Extract calibration-in-the-large (intercept) and calibration slope
        # rms::calibrate returns: predicted x, calibrated x, calibrated y, ...
        # We fit a logistic recalibration model:
        #   observed ~ logit(predicted)
        # to get slope and intercept
        pred_probs <- cal[,"predicted"]
        obs_probs  <- cal[,"calibrated.corrected"]
        valid_idx   <- !is.na(pred_probs) & !is.na(obs_probs)
        if (sum(valid_idx) >= 5) {
          cal_model <- glm(obs_probs ~ logit(pred_probs),
                           data = data.frame(pred_probs = pred_probs[valid_idx],
                                             obs_probs = obs_probs[valid_idx]),
                           family = gaussian)
          cal_metrics$CalibrationSlope[i]     <- round(coef(cal_model)[2], 3)
          cal_metrics$CalibrationIntercept[i] <- round(coef(cal_model)[1], 3)
          cat(sprintf("    %d-month: slope=%.3f, intercept=%.3f\n",
                      st, cal_metrics$CalibrationSlope[i],
                      cal_metrics$CalibrationIntercept[i]))
        }

        # --- 9b. Brier score via pec ---
        # pec::pec computes prediction error curves (Brier score)
        if (requireNamespace("pec", quietly = TRUE)) {
          suppressPackageStartupMessages(library(pec))
          # Build a simple coxph from the final formula for pec compatibility
          cox_for_pec <- coxph(final_formula, data = cox_df, x = TRUE, y = TRUE)
          # Compute apparent Brier score at time st
          pe <- tryCatch({
            pec::pec(object = list("Model" = cox_for_pec),
                     formula = final_formula,
                     data = cox_df,
                     times = st,
                     start = 1,
                     exact = FALSE,
                     reference = FALSE)
          }, error = function(e) NULL)
          if (!is.null(pe) && !is.null(pe$AppErr$Model)) {
            cal_metrics$BrierScore[i] <- round(pe$AppErr$Model, 4)
            cat(sprintf("    %d-month Brier: %.4f\n", st, cal_metrics$BrierScore[i]))
          }
        } else {
          cat(sprintf("    %d-month Brier: pec package not available, skipped\n", st))
        }
      }, error = function(e) {
        cat(sprintf("  %d-month calibration metrics failed: %s\n", st, e$message))
      })
    }

    # Save calibration metrics
    write.csv(cal_metrics, file.path(NOMO_TAB_DIR, "calibration_metrics.csv"),
              row.names = FALSE)
    cat(sprintf("  -> Saved: %s\n",
                file.path(NOMO_TAB_DIR, "calibration_metrics.csv")))
  }, error = function(e) {
    cat(sprintf("  WARNING: Nomogram failed: %s\n", e$message))
  })
} else {
  cat("  Skipping Nomogram: <2 variables in final model\n")
}

# ============================================================
# 10. Bootstrap model stability
# ============================================================
cat("\n=== [10] Bootstrap model stability ===\n\n")

B_stab <- 500
stab_results <- data.frame(matrix(NA, nrow = B_stab, ncol = length(final_vars)))
colnames(stab_results) <- final_vars

cat(sprintf("  Running %d bootstrap Cox iterations...\n", B_stab))
for (b in seq_len(B_stab)) {
  if (b %% 100 == 0) cat(sprintf("    Iteration %d / %d\n", b, B_stab))
  set.seed(b)
  boot_idx <- sample(nrow(cox_df), replace = TRUE)
  boot_df <- cox_df[boot_idx, ]
  tryCatch({
    boot_cox <- coxph(final_formula, data = boot_df)
    boot_coef <- coef(boot_cox)
    stab_results[b, names(boot_coef)] <- boot_coef
  }, error = function(e) {})
}

# Coefficient stability summary
cat("\n  Bootstrap coefficient stability (B=500):\n")
cat(sprintf("  %-55s %8s %8s %8s %8s\n",
            "Variable", "Mean", "SD", "CI_low", "CI_high"))
for (v in final_vars) {
  valid_coef <- stab_results[[v]][!is.na(stab_results[[v]])]
  if (length(valid_coef) > 10) {
    ci <- quantile(valid_coef, c(0.025, 0.975))
    cat(sprintf("  %-55s %7.3f %7.3f %7.3f %7.3f\n",
                substr(v, 1, 55),
                mean(valid_coef), sd(valid_coef), ci[1], ci[2]))

    # Direction consistency
    pct_pos <- mean(valid_coef > 0) * 100
    pct_neg <- mean(valid_coef < 0) * 100
    direction_stab <- max(pct_pos, pct_neg)
    cat(sprintf("  %55s Direction stability: %.1f%%\n", "", direction_stab))
  }
}

# Save stability
write.csv(stab_results, file.path(NOMO_TAB_DIR, "bootstrap_cox_stability.csv"),
          row.names = FALSE)
cat(sprintf("\n  -> Saved: %s\n",
            file.path(NOMO_TAB_DIR, "bootstrap_cox_stability.csv")))

# ============================================================
# 11. Save risk scores
# ============================================================
cat("\n=== [11] Saving results ===\n\n")

risk_df <- cox_df[, c("sample_id", "title", "rfs_event", "rfs_delay",
                        "risk_score", "risk_group")]
write.csv(risk_df, file.path(NOMO_TAB_DIR, "risk_scores.csv"), row.names = FALSE)
cat(sprintf("  -> Saved: %s\n",
            file.path(NOMO_TAB_DIR, "risk_scores.csv")))

# ============================================================
# 12. Summary report
# ============================================================
cat("\n=== [12] Analysis Summary ===\n\n")

sink(file.path(NOMO_TAB_DIR, "nomogram_analysis_report.txt"))

cat("=== Nomogram & Cox Prognostic Model Report ===\n")
cat("==============================================\n\n")
cat(sprintf("Cohort: GSE39582 XELOX subgroup\n"))
cat(sprintf("Samples with valid RFS: %d\n", nrow(cox_df)))
cat(sprintf("RFS events: %d / %d (%.1f%%)\n",
            sum(cox_df$rfs_event), nrow(cox_df), mean(cox_df$rfs_event)*100))
cat(sprintf("Median follow-up: %.1f months\n\n", median(cox_df$rfs_delay)))

cat("--- Model strategy ---\n")
cat("LASSO-Cox selected 0 pathways (survival signal too weak for L1)\n")
cat("Alternative: Univariate p<0.05 filter -> backward AIC selection\n\n")

cat("--- Univariate Cox (top 10) ---\n")
for (i in 1:min(10, nrow(uni_results))) {
  r <- uni_results[i, ]
  cat(sprintf("  %s: HR=%.3f (%.3f-%.3f), p=%.4f\n",
              r$Pathway, r$HR, r$HR_lower, r$HR_upper, r$p_value))
}

cat(sprintf("\n--- Pathway candidates (p<0.05): %d ---\n", length(candidate_pathways)))
for (pw in candidate_pathways) {
  cat(sprintf("  %s\n", pw))
}

cat(sprintf("\n--- Final model variables ---\n"))
cat(deparse(final_formula, width.cutoff = 60))
cat("\n")
cat(sprintf("\nClinical-only pathway candidates (p<0.05): %d\n", length(candidate_pathways)))
cat(sprintf("Pathway-only C-index: %.3f\n", c_index_pw))
cat(sprintf("Clinical-only C-index: %.3f\n", c_index_clin))
cat(sprintf("Combined C-index:      %.3f\n", c_index_full))

cat(sprintf("\n--- Performance ---\n"))
cat(sprintf("Final C-index: %.3f\n", c_index_val))
cat(sprintf("Log-rank p: %.4f\n", lr_p))

cat("\n--- Bootstrap stability (B=500) ---\n")
for (v in final_vars) {
  valid_coef <- stab_results[[v]][!is.na(stab_results[[v]])]
  if (length(valid_coef) > 10) {
    ci <- quantile(valid_coef, c(0.025, 0.975))
    pct_pos <- mean(valid_coef > 0) * 100
    direction_stab <- max(pct_pos, 100 - pct_pos)
    cat(sprintf("  %s: mean=%.3f, 95%%CI=[%.3f,%.3f], dir_stab=%.0f%%\n",
                v, mean(valid_coef), ci[1], ci[2], direction_stab))
  }
}

sink()
cat(sprintf("  -> Saved: %s\n",
            file.path(NOMO_TAB_DIR, "nomogram_analysis_report.txt")))

cat("\n=== Script 29 complete ===\n\n")
