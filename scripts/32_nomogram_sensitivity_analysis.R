# ============================================================
# 32_NOMOGRAM_SENSITIVITY_ANALYSIS.R
# Script 29 Nomogram Model Sensitivity Diagnostics
#
# 6 Analyses:
#   1. PH assumption test (Schoenfeld residuals, cox.zph)
#   2. Influential observation diagnosis (dfbeta)
#   3. Leave-one-variable-out C-index delta
#   4. Stratified C-index by clinical subgroups
#   5. Risk cutoff sensitivity sweep
#   6. Bootstrap C-index CI (B=1000)
#
# Input: results/tables/GSE39582_xelox_groups.csv
#        results/tables/pathway_activity/GSE39582_pathway_scores.rds
#        results/tables/nomogram/risk_scores.csv
# Output: results/tables/nomogram_sensitivity/
#         results/figures/nomogram_sensitivity/
# ============================================================

Sys.setenv(TMPDIR = "/tmp", TMP = "/tmp", TEMP = "/tmp")
.libPaths(c("/path/to/Rlibs", .libPaths()))

PROJECT_ROOT <- "/path/to/xelox_project基于可解释性机器学习的XELOX耐药分子指纹研究"

RESULTS_TAB_DIR <- file.path(PROJECT_ROOT, "results", "tables")
RESULTS_FIG_DIR <- file.path(PROJECT_ROOT, "results", "figures")
PA_DIR         <- file.path(RESULTS_TAB_DIR, "pathway_activity")
NOMO_TAB_DIR   <- file.path(RESULTS_TAB_DIR, "nomogram")
SENS_TAB_DIR   <- file.path(RESULTS_TAB_DIR, "nomogram_sensitivity")
SENS_FIG_DIR   <- file.path(RESULTS_FIG_DIR, "nomogram_sensitivity")

suppressWarnings(dir.create(SENS_TAB_DIR, showWarnings = FALSE, recursive = TRUE))
suppressWarnings(dir.create(SENS_FIG_DIR, showWarnings = FALSE, recursive = TRUE))

set.seed(42)

cat("============================================================\n")
cat("Script 32: Nomogram Model Sensitivity Diagnostics\n")
cat("============================================================\n\n")

# ============================================================
# 0. Load packages
# ============================================================
cat("=== [0] Loading packages ===\n\n")
required_pkgs <- c("survival", "survminer", "rms", "ggplot2", "data.table")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, lib = "/path/to/Rlibs", repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  cat(sprintf("  %s v%s loaded\n", pkg, packageVersion(pkg)))
}
cat("\n")

# ============================================================
# 1. Load data & rebuild final model
# ============================================================
cat("=== [1] Loading data & final model ===\n\n")

# --- Clinical ---
clin <- read.csv(file.path(RESULTS_TAB_DIR, "GSE39582_xelox_groups.csv"),
                 stringsAsFactors = FALSE)

# --- Pathway scores ---
ps <- readRDS(file.path(PA_DIR, "GSE39582_pathway_scores.rds"))
ps_t <- t(ps)
ps_df <- as.data.frame(ps_t)
ps_df$sample_id <- rownames(ps_t)

# --- Merge ---
cox_df <- merge(clin, ps_df, by = "sample_id", all.x = TRUE)

# Filter valid RFS
cox_df <- cox_df[!is.na(cox_df$rfs_event) & !is.na(cox_df$rfs_delay) &
                   cox_df$rfs_delay > 0, ]
cox_df$rfs_event <- as.numeric(cox_df$rfs_event)
cox_df$rfs_delay <- as.numeric(cox_df$rfs_delay)

cat(sprintf("  N = %d, Events = %d (%.1f%%)\n",
            nrow(cox_df), sum(cox_df$rfs_event), mean(cox_df$rfs_event) * 100))

# --- Create derived clinical variables (same as Script 29) ---
cox_df$sex_bin <- ifelse(cox_df$sex == "Male", 1, 0)
cox_df$age_num <- as.numeric(cox_df$age)
cox_df$stage_III_IV <- ifelse(grepl("3|4|III|IV", cox_df$tnm_stage), 1, 0)
cox_df$mmr_dmmr <- ifelse(cox_df$mmr_status == "dMMR", 1, 0)
cox_df$location_distal <- ifelse(cox_df$tumor_location == "distal", 1, 0)
cat(sprintf("  location_distal: distal=%d, proximal=%d\n",
            sum(cox_df$location_distal == 1), sum(cox_df$location_distal == 0)))

# Standardize pathways
pathway_names <- rownames(ps)
pw_matrix <- as.matrix(cox_df[, pathway_names, drop = FALSE])
pw_scaled <- scale(pw_matrix)

# --- Final 8 variables (from Script 29 backward AIC results) ---
#   7 pathways + location_distal
final_vars <- c(
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "KEGG_ECM_RECEPTOR_INTERACTION",
  "KEGG_TGF_BETA_SIGNALING_PATHWAY",
  "KEGG_PATHWAYS_IN_CANCER",
  "HALLMARK_MYC_TARGETS_V2",
  "KEGG_COLORECTAL_CANCER",
  "location_distal"
)

# Build final formula
final_formula_str <- paste0(
  "Surv(rfs_delay, rfs_event) ~ ",
  paste0("`", final_vars, "`", collapse = " + ")
)
final_formula <- as.formula(final_formula_str)

final_fit <- coxph(final_formula, data = cox_df, x = TRUE, y = TRUE)
cat(sprintf("\n  Final model: C-index = %.3f\n",
            summary(final_fit)$concordance["C"]))

# Build risk scores
risk_scores <- predict(final_fit, type = "risk")
cox_df$risk_score <- risk_scores
cox_df$risk_group <- ifelse(risk_scores > median(risk_scores), "High", "Low")

cat(sprintf("  Risk score range: %.3f to %.3f\n", min(risk_scores), max(risk_scores)))
cat(sprintf("  High risk: %d, Low risk: %d\n",
            sum(cox_df$risk_group == "High"), sum(cox_df$risk_group == "Low")))

# ============================================================
# ANALYSIS 1: PH Assumption Test (Schoenfeld Residuals)
# ============================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("[1] Proportional Hazards Assumption Test\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

zph_test <- cox.zph(final_fit)
cat("  Global PH test:\n")
cat(sprintf("  Chisq = %.3f, df = %d, p = %.4f\n",
            zph_test$table["GLOBAL", "chisq"],
            zph_test$table["GLOBAL", "df"],
            zph_test$table["GLOBAL", "p"]))

cat("\n  Per-variable PH test:\n")
for (var in rownames(zph_test$table)) {
  if (var == "GLOBAL") next
  r <- zph_test$table[var, ]
  flag <- if (r["p"] < 0.05) " ⚠️" else " ✅"
  cat(sprintf("  %-45s rho=%7.3f chisq=%6.2f p=%6.4f%s\n",
              substr(var, 1, 45), r["rho"], r["chisq"], r["p"], flag))
}

# Save table
ph_table <- as.data.frame(zph_test$table)
ph_table$variable <- rownames(ph_table)
rownames(ph_table) <- NULL
ph_table <- ph_table[ph_table$variable != "GLOBAL", ]
write.csv(ph_table, file.path(SENS_TAB_DIR, "ph_assumption_test.csv"), row.names = FALSE)

# Plot Schoenfeld residuals
pdf(file.path(SENS_FIG_DIR, "schoenfeld_residuals.pdf"), width = 14, height = 10)
plot(zph_test, df = 4, resid = TRUE)
invisible(dev.off())
# Also single-variable plots for the most borderline ones
for (var in rownames(zph_test$table)) {
  if (var == "GLOBAL") next
  p_val <- zph_test$table[var, "p"]
  if (p_val < 0.10) {  # borderline variables
    pdf(file.path(SENS_FIG_DIR, sprintf("schoenfeld_%s.pdf", gsub("`", "", var))),
        width = 8, height = 6)
    plot(zph_test, var = var, resid = TRUE, df = 4)
    invisible(dev.off())
  }
}
cat(sprintf("  -> Saved: %s\n", file.path(SENS_FIG_DIR, "schoenfeld_residuals.pdf")))

# ============================================================
# ANALYSIS 2: Influential Observation Diagnosis (dfbeta)
# ============================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("[2] Influential Observation Diagnosis (dfbeta)\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Compute dfbeta for final model (use stats::resid to avoid rms masking)
db <- stats::resid(final_fit, type = "dfbeta")

# Summary: max dfbeta per variable
varnames <- colnames(db)
if (is.null(varnames)) {
  varnames <- dimnames(final_fit$x)[[2]]
}
cat(sprintf("  dfbeta matrix: %d obs x %d vars\n", nrow(db), ncol(db)))
dfbeta_max <- data.frame(
  variable = varnames,
  max_abs_dfbeta = apply(abs(db), 2, max),
  n_influential = apply(abs(db), 2, function(x) sum(x > 2 / sqrt(nrow(cox_df)))),
  pct_influential = apply(abs(db), 2, function(x) mean(x > 2 / sqrt(nrow(cox_df))) * 100)
)
threshold <- 2 / sqrt(nrow(cox_df))
cat(sprintf("  Influential threshold: dfbeta > %.4f (2/sqrt(n))\n", threshold))
cat(sprintf("  %-45s %10s %10s %10s\n", "Variable", "Max|dfbeta|", "N>threshold", "%>threshold"))
for (i in 1:nrow(dfbeta_max)) {
  r <- dfbeta_max[i, ]
  cat(sprintf("  %-45s %10.4f %10d %9.1f%%\n",
              substr(r$variable, 1, 45),
              r$max_abs_dfbeta, r$n_influential, r$pct_influential))
}

write.csv(dfbeta_max, file.path(SENS_TAB_DIR, "dfbeta_summary.csv"), row.names = FALSE)

# Plot dfbeta for key variables
final_vars_short <- gsub("`", "", final_vars)
n_vars <- length(final_vars_short)
pdf(file.path(SENS_FIG_DIR, "dfbeta_index_plot.pdf"), width = 12, height = ceiling(n_vars / 2) * 3)
par(mfrow = c(ceiling(n_vars / 2), 2), mar = c(4, 4, 3, 1))
for (i in seq_along(final_vars_short)) {
  if (final_vars_short[i] %in% colnames(db)) {
    db_i <- db[, final_vars_short[i]]
    plot(db_i, type = "h", ylab = "dfbeta", xlab = "Observation index",
         main = final_vars_short[i], col = ifelse(abs(db_i) > threshold, "red", "gray50"))
    abline(h = c(-threshold, threshold), lty = 2, col = "red")
  }
}
invisible(dev.off())
cat(sprintf("  -> Saved: %s\n", file.path(SENS_FIG_DIR, "dfbeta_index_plot.pdf")))

# ============================================================
# ANALYSIS 3: Leave-One-Variable-Out C-index Delta
# ============================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("[3] Leave-One-Variable-Out C-index Delta\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

base_c <- summary(final_fit)$concordance["C"]
cat(sprintf("  Baseline C-index (8 variables): %.4f\n\n", base_c))

lovo_results <- data.frame(
  variable_removed = final_vars,
  c_index = NA,
  delta_c = NA,
  stringsAsFactors = FALSE
)

for (i in seq_along(final_vars)) {
  vars_kept <- final_vars[-i]
  fm <- as.formula(paste0(
    "Surv(rfs_delay, rfs_event) ~ ",
    paste0("`", vars_kept, "`", collapse = " + ")
  ))
  fit_reduced <- coxph(fm, data = cox_df)
  c_reduced <- summary(fit_reduced)$concordance["C"]
  lovo_results$c_index[i] <- c_reduced
  lovo_results$delta_c[i] <- c_reduced - base_c
}

lovo_results <- lovo_results[order(lovo_results$delta_c), ]
cat(sprintf("  %-45s %10s %10s %10s\n", "Removed Variable", "C-index", "Δ C", "Rank"))
for (i in 1:nrow(lovo_results)) {
  r <- lovo_results[i, ]
  flag <- if (r$delta_c < -0.005) "🔴" else if (r$delta_c < 0) "🟡" else "🟢"
  cat(sprintf("  %s %-45s %9.4f %9.4f %s\n",
              flag, substr(r$variable_removed, 1, 45),
              r$c_index, r$delta_c, flag))
}

write.csv(lovo_results, file.path(SENS_TAB_DIR, "leave_one_variable_out.csv"), row.names = FALSE)

# Plot
pdf(file.path(SENS_FIG_DIR, "leave_one_variable_out_cindex.pdf"), width = 10, height = 6)
lov_plot <- lovo_results
lov_plot$variable_removed <- factor(lov_plot$variable_removed,
                                     levels = lov_plot$variable_removed)
p_lovo <- ggplot(lov_plot, aes(x = delta_c, y = variable_removed)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_segment(aes(xend = 0, yend = variable_removed), linewidth = 0.5, color = "gray60") +
  geom_point(aes(color = delta_c < 0), size = 3) +
  scale_color_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"),
                     labels = c("TRUE" = "negative Δ", "FALSE" = "positive Δ")) +
  labs(title = "Leave-One-Variable-Out: Δ C-index",
       subtitle = sprintf("Baseline C = %.4f (8-variable model)", base_c),
       x = "Delta C-index (removed - full)",
       y = "") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")
print(p_lovo)
invisible(dev.off())
cat(sprintf("  -> Saved: %s\n", file.path(SENS_FIG_DIR, "leave_one_variable_out_cindex.pdf")))

# ============================================================
# ANALYSIS 4: Stratified C-index by Clinical Subgroups
# ============================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("[4] Stratified C-index by Clinical Subgroups\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

stratify_cindex <- function(df, col, label_fn = identity) {
  result <- data.frame(
    subgroup = character(),
    n = integer(),
    n_events = integer(),
    c_index = numeric(),
    stringsAsFactors = FALSE
  )
  for (grp in sort(unique(df[[col]]))) {
    if (is.na(grp)) next
    sub <- df[df[[col]] == grp & !is.na(df[[col]]), ]
    if (nrow(sub) < 20) next
    fit_sub <- tryCatch({
      coxph(final_formula, data = sub)
    }, error = function(e) NULL)
    if (!is.null(fit_sub)) {
      c_sub <- summary(fit_sub)$concordance["C"]
      result <- rbind(result, data.frame(
        subgroup = label_fn(grp),
        n = nrow(sub),
        n_events = sum(sub$rfs_event),
        c_index = round(c_sub, 4),
        stringsAsFactors = FALSE
      ))
    }
  }
  return(result)
}

# 4a. Stage (use tnm_stage, merge small groups)
cox_df$stage_group <- ifelse(grepl("1|I", cox_df$tnm_stage) & !grepl("2|3|4|II|III|IV", cox_df$tnm_stage), "Stage I",
                      ifelse(grepl("2|II", cox_df$tnm_stage) & !grepl("3|4|III|IV", cox_df$tnm_stage), "Stage II",
                      ifelse(grepl("3|III", cox_df$tnm_stage), "Stage III",
                      ifelse(grepl("4|IV", cox_df$tnm_stage), "Stage IV", "Unknown"))))
stage_cindex <- stratify_cindex(cox_df, "stage_group", identity)
cat("  --- Stage subgroups ---\n")
print(stage_cindex, row.names = FALSE)

# 4b. Tumor location
loc_cindex <- stratify_cindex(cox_df, "tumor_location",
                               function(x) ifelse(x == "distal", "Distal (Left)", "Proximal (Right)"))
cat("\n  --- Tumor location ---\n")
print(loc_cindex, row.names = FALSE)

# 4c. MMR status (MSI surrogate)
msi_cindex <- stratify_cindex(cox_df, "mmr_status",
                               function(x) ifelse(x == "dMMR", "dMMR", "pMMR"))
cat("\n  --- MMR status ---\n")
print(msi_cindex, row.names = FALSE)

# 4d. Adjuvant vs non-adjuvant chemotherapy
adj_cindex <- NULL
if ("adjuvant_chemo" %in% colnames(cox_df)) {
  cox_df$setting <- ifelse(cox_df$adjuvant_chemo == 1 | grepl("Adjuvant|adjuvant", cox_df$adjuvant_chemo),
                           "Adjuvant", "Non-Adjuvant")
  adj_cindex <- stratify_cindex(cox_df, "setting", identity)
} else if ("chemo_type" %in% colnames(cox_df)) {
  adj_cindex <- stratify_cindex(cox_df, "chemo_type", identity)
}
if (!is.null(adj_cindex)) {
  cat("\n  --- Treatment setting ---\n")
  print(adj_cindex, row.names = FALSE)
}

# Combine
all_strat <- rbind(
  cbind(data.frame(dimension = "Stage", stringsAsFactors = FALSE), stage_cindex),
  cbind(data.frame(dimension = "Location", stringsAsFactors = FALSE), loc_cindex),
  cbind(data.frame(dimension = "MSI", stringsAsFactors = FALSE), msi_cindex)
)
if (!is.null(adj_cindex) && nrow(adj_cindex) > 0) {
  all_strat <- rbind(all_strat,
    cbind(data.frame(dimension = "Setting", stringsAsFactors = FALSE), adj_cindex))
}
write.csv(all_strat, file.path(SENS_TAB_DIR, "stratified_cindex.csv"), row.names = FALSE)

# Plot
pdf(file.path(SENS_FIG_DIR, "stratified_cindex.pdf"), width = 12, height = 5.5)
p_strat <- ggplot(all_strat, aes(x = reorder(paste(dimension, subgroup), -c_index),
                                  y = c_index, fill = dimension)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.65) +
  geom_text(aes(label = sprintf("%.3f\n(n=%d)", c_index, n)),
            hjust = -0.1, size = 3.2) +
  facet_wrap(~ dimension, scales = "free_x", nrow = 1) +
  coord_flip(ylim = c(0.3, 0.95)) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Stratified C-index by Clinical Subgroups",
       subtitle = sprintf("Full model C-index = %.3f", base_c),
       x = "", y = "C-index") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")
print(p_strat)
invisible(dev.off())
cat(sprintf("\n  -> Saved: %s\n", file.path(SENS_FIG_DIR, "stratified_cindex.pdf")))

# ============================================================
# ANALYSIS 5: Risk Cutoff Sensitivity Sweep
# ============================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("[5] Risk Score Cutoff Sensitivity Sweep\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Use the risk scores from final model
risk <- cox_df$risk_score
time <- cox_df$rfs_delay
event <- cox_df$rfs_event

# Sweep over percentile cutoffs
pct_seq <- seq(10, 90, by = 2)
sweep_df <- data.frame(
  percentile = pct_seq,
  cutoff_value = NA,
  n_high = NA,
  n_low = NA,
  hr = NA,
  hr_lower = NA,
  hr_upper = NA,
  logrank_p = NA,
  c_index = NA,
  stringsAsFactors = FALSE
)

for (i in seq_along(pct_seq)) {
  pct <- pct_seq[i]
  cut <- quantile(risk, pct / 100)
  group <- ifelse(risk > cut, "High", "Low")

  n_hi <- sum(group == "High")
  n_lo <- sum(group == "Low")
  if (min(n_hi, n_lo) < 10) next

  sweep_df$cutoff_value[i] <- cut
  sweep_df$n_high[i] <- n_hi
  sweep_df$n_low[i] <- n_lo

  # Cox HR
  tryCatch({
    fit_cut <- coxph(Surv(time, event) ~ factor(group, levels = c("Low", "High")),
                     data = data.frame(time, event, group))
    s <- summary(fit_cut)
    sweep_df$hr[i] <- exp(coef(fit_cut))
    sweep_df$hr_lower[i] <- exp(confint(fit_cut)[1])
    sweep_df$hr_upper[i] <- exp(confint(fit_cut)[2])
    sweep_df$logrank_p[i] <- s$logtest["pvalue"]
    sweep_df$c_index[i] <- s$concordance["C"]
  }, error = function(e) {
    sweep_df$hr[i] <- NA
  })
}

# Clean
sweep_df <- sweep_df[!is.na(sweep_df$hr), ]
write.csv(sweep_df, file.path(SENS_TAB_DIR, "cutoff_sweep.csv"), row.names = FALSE)

cat(sprintf("  %5s %8s %8s %8s %8s %10s %8s\n",
            "Pct", "Cut", "N_High", "N_Low", "HR", "Log-rank p", "C-index"))
for (i in 1:nrow(sweep_df)) {
  r <- sweep_df[i, ]
  cat(sprintf("  %4d%% %8.3f %8d %8d %8.2f %10.4f %8.3f\n",
              r$percentile, r$cutoff_value, r$n_high, r$n_low,
              r$hr, r$logrank_p, r$c_index))
}

# Plot: dual y-axis (HR + log-rank p)
pdf(file.path(SENS_FIG_DIR, "cutoff_sweep.pdf"), width = 12, height = 5)

par(mar = c(4, 4, 3, 4))
# Panel A: HR
plot(sweep_df$percentile, sweep_df$hr, type = "b", pch = 16, col = "darkred",
     xlab = "Percentile for High-Risk Cutoff", ylab = "Hazard Ratio (High vs Low)",
     main = "Risk Score Cutoff Sensitivity Sweep",
     ylim = c(0.5, max(sweep_df$hr_upper, na.rm = TRUE) * 1.1))
arrows(sweep_df$percentile, sweep_df$hr_lower,
       sweep_df$percentile, sweep_df$hr_upper,
       angle = 90, code = 3, length = 0.03, col = "darkred")
abline(h = 1, lty = 2, col = "gray50")

# Add log-rank p on same plot as text
text(sweep_df$percentile, sweep_df$hr_upper,
     labels = sprintf("p=%.4f", sweep_df$logrank_p),
     cex = 0.65, pos = 3, col = "gray30")

invisible(dev.off())

# Second plot: C-index sweep
pdf(file.path(SENS_FIG_DIR, "cutoff_cindex_sweep.pdf"), width = 10, height = 5)
plot(sweep_df$percentile, sweep_df$c_index, type = "b", pch = 16, col = "steelblue",
     xlab = "Percentile for High-Risk Cutoff", ylab = "C-index",
     main = "C-index vs Risk Score Cutoff")
abline(h = base_c, lty = 2, col = "gray50", lwd = 1.5)
text(15, base_c, sprintf("Full model C = %.3f", base_c), pos = 3, cex = 0.8, col = "gray50")
invisible(dev.off())
cat(sprintf("  -> Saved: %s\n", file.path(SENS_FIG_DIR, "cutoff_sweep.pdf")))

# ============================================================
# ANALYSIS 6: Bootstrap C-index CI (B=1000) — corrected optimism
# ============================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("[6] Bootstrap C-index with proper optimism correction (B=1000)\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

B <- 1000
n <- nrow(cox_df)
boot_c_train <- numeric(B)   # C-index on bootstrap sample (training)
boot_c_test  <- numeric(B)   # C-index on original data using bootstrap model (test)

cat(sprintf("  Running %d bootstrap iterations...\n", B))

pb <- txtProgressBar(min = 0, max = B, style = 3)
for (b in 1:B) {
  setTxtProgressBar(pb, b)

  # Bootstrap sample (with replacement)
  idx <- sample(1:n, n, replace = TRUE)
  boot_df <- cox_df[idx, ]

  # Out-of-bootstrap (OOB) samples = original data not selected in bootstrap
  oob_idx <- setdiff(1:n, unique(idx))

  # Ensure at least some events
  if (sum(boot_df$rfs_event) < 3) {
    boot_c_train[b] <- NA
    boot_c_test[b]  <- NA
    next
  }

  fit_boot <- tryCatch({
    coxph(final_formula, data = boot_df)
  }, error = function(e) NULL)

  if (!is.null(fit_boot)) {
    # Training performance: model on bootstrap sample, evaluated on bootstrap sample
    boot_c_train[b] <- summary(fit_boot)$concordance["C"]

    # Test performance: model on bootstrap sample, evaluated on OOB sample
    if (length(oob_idx) >= 10 && sum(cox_df$rfs_event[oob_idx]) >= 2) {
      lp_oob <- predict(fit_boot, newdata = cox_df[oob_idx, ], type = "lp")
      # Compute concordance on OOB
      c_oob <- tryCatch({
        surv_obj_oob <- Surv(cox_df$rfs_delay[oob_idx], cox_df$rfs_event[oob_idx])
        rcorr.cens(lp_oob, surv_obj_oob)["C Index"]
      }, error = function(e) NA)
      boot_c_test[b] <- as.numeric(1 - c_oob)
    } else {
      boot_c_test[b] <- NA
    }
  } else {
    boot_c_train[b] <- NA
    boot_c_test[b]  <- NA
  }
}
close(pb)

# Valid iterations
valid_idx <- !is.na(boot_c_train) & !is.na(boot_c_test)
boot_c_train_valid <- boot_c_train[valid_idx]
boot_c_test_valid  <- boot_c_test[valid_idx]
optimism_vec <- boot_c_train_valid - boot_c_test_valid

# Summary statistics
mean_train    <- mean(boot_c_train_valid)
mean_test     <- mean(boot_c_test_valid)
mean_optimism <- mean(optimism_vec)
corrected_c   <- base_c - mean_optimism

cat(sprintf("\n  Valid bootstrap iterations: %d / %d\n", sum(valid_idx), B))
cat(sprintf("  Apparent C-index:          %.4f\n", base_c))
cat(sprintf("  Bootstrap mean (train):    %.4f\n", mean_train))
cat(sprintf("  Bootstrap mean (test/OOB): %.4f\n", mean_test))
cat(sprintf("  Mean optimism:             %.4f\n", mean_optimism))
cat(sprintf("  Optimism-corrected C:      %.4f\n", corrected_c))

boot_summary <- data.frame(
  metric = c("Apparent", "Bootstrap_Mean_Train", "Bootstrap_Mean_Test",
             "Optimism", "Optimism_Corrected", "95%CI_Lower", "95%CI_Upper"),
  value = c(
    base_c,
    mean_train,
    mean_test,
    mean_optimism,
    corrected_c,
    quantile(boot_c_test_valid, 0.025, names = FALSE),
    quantile(boot_c_test_valid, 0.975, names = FALSE)
  )
)
cat("\n  Bootstrap C-index summary (corrected):\n")
print(boot_summary, row.names = FALSE)

write.csv(boot_summary, file.path(SENS_TAB_DIR, "bootstrap_cindex_ci.csv"), row.names = FALSE)

# Plot bootstrap distribution
pdf(file.path(SENS_FIG_DIR, "bootstrap_cindex_distribution.pdf"), width = 8, height = 5)
hist(boot_c_valid, breaks = 50, col = "steelblue", border = "white",
     main = sprintf("Bootstrap C-index Distribution (B=%d)", B),
     xlab = "C-index", xlim = c(0.55, 0.80))
abline(v = base_c, col = "darkred", lwd = 2, lty = 2)
ci_lines <- quantile(boot_c_valid, c(0.025, 0.975), names = FALSE)
abline(v = ci_lines, col = "gray40", lwd = 1.5, lty = 3)
legend("topleft",
       legend = c(sprintf("Apparent C = %.3f", base_c),
                  sprintf("95%% CI: [%.3f, %.3f]", ci_lines[1], ci_lines[2]),
                  sprintf("Mean = %.3f, SD = %.3f", mean(boot_c_valid), sd(boot_c_valid))),
       col = c("darkred", "gray40", "steelblue"),
       lty = c(2, 3, 0), pch = c(NA, NA, 15), pt.cex = 2,
       bty = "n", cex = 0.9)
invisible(dev.off())
cat(sprintf("  -> Saved: %s\n", file.path(SENS_FIG_DIR, "bootstrap_cindex_distribution.pdf")))

# ============================================================
# 7. Summary Report
# ============================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("SUMMARY: Nomogram Model Sensitivity Diagnostics\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat(sprintf("  Final model: %d variables, C-index = %.4f\n", length(final_vars), base_c))
cat(sprintf("  N = %d, Events = %d (%.1f%%)\n\n",
            n, sum(cox_df$rfs_event), mean(cox_df$rfs_event) * 100))

cat("  [1] PH assumption:\n")
cat(sprintf("      Global p = %.4f\n", zph_test$table["GLOBAL", "p"]))
n_violated <- sum(ph_table$p < 0.05, na.rm = TRUE)
cat(sprintf("      Variables with PH violation (p<0.05): %d / %d\n\n",
            n_violated, nrow(ph_table)))

cat("  [2] Influential observations:\n")
cat(sprintf("      Max dfbeta: %.4f (%s)\n",
            dfbeta_max$max_abs_dfbeta[1], dfbeta_max$variable[1]))
cat(sprintf("      Observations exceeding threshold: %d (%.1f%%)\n\n",
            max(dfbeta_max$n_influential), max(dfbeta_max$pct_influential)))

cat("  [3] Leave-one-variable-out C-index:\n")
cat(sprintf("      Most influential variable: %s (ΔC = %.4f)\n",
            lovo_results$variable_removed[1], lovo_results$delta_c[1]))
cat(sprintf("      Least influential variable: %s (ΔC = %.4f)\n\n",
            lovo_results$variable_removed[nrow(lovo_results)],
            lovo_results$delta_c[nrow(lovo_results)]))

cat("  [4] Stratified C-index:\n")
for (dim_level in unique(all_strat$dimension)) {
  sub <- all_strat[all_strat$dimension == dim_level, ]
  best <- sub[which.max(sub$c_index), ]
  worst <- sub[which.min(sub$c_index), ]
  cat(sprintf("      %s: best %s (%.3f), worst %s (%.3f)\n",
              dim_level, best$subgroup, best$c_index, worst$subgroup, worst$c_index))
}

cat("\n  [5] Cutoff sensitivity:\n")
cat(sprintf("      Best cutoff: %d%% (HR=%.2f, C=%.3f, p=%.4f)\n",
            sweep_df$percentile[which.max(sweep_df$c_index)],
            sweep_df$hr[which.max(sweep_df$c_index)],
            sweep_df$c_index[which.max(sweep_df$c_index)],
            sweep_df$logrank_p[which.max(sweep_df$c_index)]))

cat("\n  [6] Bootstrap C-index (corrected):\n")
cat(sprintf("      Apparent C-index:       %.4f\n", base_c))
cat(sprintf("      Bootstrap mean (train): %.4f\n", mean_train))
cat(sprintf("      Bootstrap mean (test):  %.4f\n", mean_test))
cat(sprintf("      Mean optimism:          %.4f\n", mean_optimism))
cat(sprintf("      Optimism-corrected C:   %.4f\n", corrected_c))
cat(sprintf("      95%% CI (test):          [%.4f, %.4f]\n\n",
            quantile(boot_c_test_valid, 0.025, names = FALSE),
            quantile(boot_c_test_valid, 0.975, names = FALSE)))

# Write full report
sink(file.path(SENS_TAB_DIR, "sensitivity_report.txt"))
cat("Nomogram Model Sensitivity Analysis Report\n")
cat("==========================================\n\n")
cat(sprintf("Date: %s\n", Sys.time()))
cat(sprintf("Model: %d-variable Cox PH (Script 29)\n", length(final_vars)))
cat(sprintf("Cohort: GSE39582 XELOX, N=%d, Events=%d\n\n", n, sum(cox_df$rfs_event)))

cat("1. PH Assumption Test\n")
cat("--------------------\n")
print(ph_table)

cat("\n2. Influential Observations (dfbeta)\n")
cat("-----------------------------------\n")
print(dfbeta_max)

cat("\n3. Leave-One-Variable-Out C-index\n")
cat("--------------------------------\n")
print(lovo_results[, c("variable_removed", "c_index", "delta_c")])

cat("\n4. Stratified C-index\n")
cat("--------------------\n")
print(all_strat)

cat("\n5. Cutoff Sensitivity\n")
cat("--------------------\n")
head(sweep_df, 5)
cat("... (truncated, see cutoff_sweep.csv)\n")

cat("\n6. Bootstrap C-index (B=1000)\n")
cat("-----------------------------\n")
print(boot_summary)
sink()

cat(sprintf("  -> Saved: %s\n\n", file.path(SENS_TAB_DIR, "sensitivity_report.txt")))

cat("=== Script 32 complete ===\n")
cat(sprintf("Output tables: %s\n", SENS_TAB_DIR))
cat(sprintf("Output figures: %s\n", SENS_FIG_DIR))
