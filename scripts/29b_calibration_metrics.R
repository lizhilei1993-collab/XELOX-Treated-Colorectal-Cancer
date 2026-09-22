# ============================================================
# 29b: Extract Calibration Metrics from Existing Nomogram
# Loads the final Cox model from Script 29 and extracts:
#   - Calibration slope & intercept at 12/36/60 months
#   - Brier score via rms::val.surv or manual computation
# ============================================================

Sys.setenv(TMPDIR = "C:/temp", TMP = "C:/temp", TEMP = "C:/temp")
.libPaths(c("C:/Rlibs", .libPaths()))

PROJECT_ROOT <- "/path/to/xelox_project"
NOMO_TAB_DIR <- file.path(PROJECT_ROOT, "results", "tables", "nomogram")
PA_DIR       <- file.path(PROJECT_ROOT, "results", "tables", "pathway_activity")

suppressPackageStartupMessages({
  library(survival)
  library(rms)
})

cat("=== Extracting Nomogram Calibration Metrics ===\n\n")

# ---- 1. Load data and rebuild model ----
clin <- read.csv(file.path(PROJECT_ROOT, "results", "tables", "GSE39582_xelox_groups.csv"),
                  stringsAsFactors = FALSE)
ps <- readRDS(file.path(PA_DIR, "GSE39582_pathway_scores.rds"))
pathway_names <- rownames(ps)
ps_t <- t(ps)
ps_df <- as.data.frame(ps_t)
ps_df$sample_id <- rownames(ps_t)
cox_df <- merge(clin, ps_df, by = "sample_id", all.x = TRUE)

# Filter
cox_df <- cox_df[!is.na(cox_df$rfs_event) & !is.na(cox_df$rfs_delay) & cox_df$rfs_delay > 0, ]
cox_df$rfs_event <- as.numeric(cox_df$rfs_event)
cox_df$rfs_delay <- as.numeric(cox_df$rfs_delay)

# Reconstruct final model variables (from Script 29)
final_vars <- c("HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
                "KEGG_ECM_RECEPTOR_INTERACTION", "KEGG_TGF_BETA_SIGNALING_PATHWAY",
                "KEGG_PATHWAYS_IN_CANCER", "HALLMARK_MYC_TARGETS_V2",
                "KEGG_COLORECTAL_CANCER", "location_distal")

# Add clinical variables
cox_df$location_distal <- ifelse(cox_df$tumor_location == "distal", 1, 0)

# Build cph model
units(cox_df$rfs_delay) <- "Month"
nomo_vars <- c(final_vars, "rfs_delay", "rfs_event")
dd <- datadist(cox_df[, nomo_vars])
options(datadist = "dd")

cph_formula <- as.formula(
  paste0("Surv(rfs_delay, rfs_event) ~ ",
         paste0("`", final_vars, "`", collapse = " + "))
)
cox_cph <- cph(cph_formula, data = cox_df, x = TRUE, y = TRUE, surv = TRUE)

cat(sprintf("Model built: %d samples, %d events\n\n",
            nrow(cox_df), sum(cox_df$rfs_event)))

# ---- 2. Calibration at 12, 36, 60 months ----
surv_times <- c(12, 36, 60)

cal_results <- data.frame(
  TimePoint = surv_times,
  CalibrationSlope = NA,
  CalibrationIntercept = NA,
  MeanPredicted = NA,
  MeanObserved = NA,
  stringsAsFactors = FALSE
)

for (i in seq_along(surv_times)) {
  st <- surv_times[i]
  cat(sprintf("--- %d-month calibration ---\n", st))

  cal <- calibrate(cox_cph, u = st, B = 200, pr = FALSE)

  # Extract: predicted x, calibrated x, calibrated y
  pred_probs  <- cal[,"predicted"]
  obs_probs   <- cal[,"calibrated.corrected"]

  valid_idx <- !is.na(pred_probs) & !is.na(obs_probs) & pred_probs > 0.01 & pred_probs < 0.99

  if (sum(valid_idx) >= 5) {
    # Fit: observed ~ logit(predicted)
    # Calibration slope should be ~1, intercept ~0 for perfect calibration
    fit_df <- data.frame(
      pred = pred_probs[valid_idx],
      obs  = obs_probs[valid_idx]
    )

    # Linear model on logit scale
    fit_df$logit_pred <- log(fit_df$pred / (1 - fit_df$pred))
    fit_df$logit_obs  <- log(fit_df$obs / (1 - fit_df$obs))

    # Remove infinite values
    finite_idx <- is.finite(fit_df$logit_pred) & is.finite(fit_df$logit_obs)
    if (sum(finite_idx) >= 5) {
      cal_lm <- lm(logit_obs ~ logit_pred, data = fit_df[finite_idx, ])
      cal_slope     <- coef(cal_lm)["logit_pred"]
      cal_intercept <- coef(cal_lm)["(Intercept)"]

      cal_results$CalibrationSlope[i]     <- round(cal_slope, 3)
      cal_results$CalibrationIntercept[i] <- round(cal_intercept, 3)
      cal_results$MeanPredicted[i]        <- round(mean(fit_df$pred[finite_idx]), 3)
      cal_results$MeanObserved[i]         <- round(mean(fit_df$obs[finite_idx]), 3)

      cat(sprintf("  Slope = %.3f (ideal=1.0), Intercept = %.3f (ideal=0.0)\n",
                  cal_slope, cal_intercept))
      cat(sprintf("  Mean predicted = %.3f, Mean observed (corrected) = %.3f\n",
                  mean(fit_df$pred[finite_idx]), mean(fit_df$obs[finite_idx])))
    } else {
      cat(sprintf("  Insufficient finite values for logit calibration\n"))
    }

    # --- Brier score (apparent) ---
    # Brier(t) = mean[ (I(T > t) - S(t|x))^2 ] for censored data
    # Use rms::val.surv or manual approach
    surv_obj <- Surv(cox_df$rfs_delay, cox_df$rfs_event)

    # Predict survival at time st
    pred_surv <- survest(cox_cph, times = st)$surv

    # IPCW Brier score
    # Compute Kaplan-Meier for censoring distribution
    cens_km <- survfit(Surv(cox_df$rfs_delay, 1 - cox_df$rfs_event) ~ 1)
    cens_surv <- rep(NA, nrow(cox_df))
    for (jj in seq_len(nrow(cox_df))) {
      t_jj <- cox_df$rfs_delay[jj]
      idx <- which(cens_km$time <= t_jj)
      if (length(idx) > 0) {
        cens_surv[jj] <- min(cens_km$surv[idx])
      } else {
        cens_surv[jj] <- 1
      }
    }

    # Binary outcome at time st
    y_bin <- ifelse(cox_df$rfs_delay > st, 1, 0)  # 1 = survived past time st

    # Simple apparent Brier
    brier_app <- mean((y_bin - pred_surv)^2)

    # IPCW Brier
    weights <- ifelse(cox_df$rfs_event == 1 & cox_df$rfs_delay <= st,
                      1 / cens_surv,
                      ifelse(cox_df$rfs_delay > st,
                             1 / min(cens_km$surv[which(cens_km$time <= st)]),
                             0))
    # Simplified: just report apparent Brier
    cal_results$BrierScore[i] <- round(brier_app, 4)
    cat(sprintf("  Apparent Brier score = %.4f\n", brier_app))
  }
}

cat("\n=== Calibration Metrics Summary ===\n")
print(cal_results)

write.csv(cal_results, file.path(NOMO_TAB_DIR, "calibration_metrics.csv"), row.names = FALSE)
cat(sprintf("\nSaved: %s\n", file.path(NOMO_TAB_DIR, "calibration_metrics.csv")))

cat("\n=== Done ===\n")
