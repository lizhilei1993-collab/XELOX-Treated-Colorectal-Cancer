# ============================================================
# 27_GENE_LASSO_COX_FAIR_COMPARISON.R
# Gene-level LASSO-Cox vs. Pathway-level LASSO-Cox Fair Comparison
#
# Purpose:
#   Isolate the unique contribution of pathway-level feature
#   engineering by applying IDENTICAL LASSO-Cox methodology to both
#   the 2,537-gene candidate pool and the 44-pathway ssGSEA scores
#   within the same GSE39582 XELOX cohort.
#
#   The original analysis confounded three factors:
#     Gene-level: 10-algorithm ensemble, 10 features, EPV=3.4
#     Pathway-level: backward AIC Cox, 8 features, EPV=9.9
#
#   This script holds modeling method constant (LASSO-Cox) and
#   isolates the feature-type effect (genes vs. pathways).
#
# Method:
#   1. Gene-level: LASSO-Cox (cv.glmnet) on 2,537 candidate genes
#   2. Pathway-level: LASSO-Cox on 44 ssGSEA pathway scores
#   3. Compute C-index (apparent + bootstrap-corrected) for both
#   4. Report EPV, selected features, coefficient stability
#   5. Save comparison table + forest plots
#
# Input:
#   results/tables/pathway_activity/GSE39582_gene_expression.rds
#   results/tables/pathway_activity/GSE39582_pathway_scores.rds
#   results/tables/GSE39582_xelox_groups.csv
#
# Output:
#   results/tables/fair_comparison/gene_lasso_cox_results.csv
#   results/tables/fair_comparison/pathway_lasso_cox_results.csv
#   results/tables/fair_comparison/fair_comparison_summary.csv
#   results/figures/fair_comparison/
# ============================================================

Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.936")
.libPaths(c("/path/to/Rlibs", .libPaths()))

PROJECT_ROOT <- "X:/"
RESULTS_TAB_DIR <- file.path(PROJECT_ROOT, "results", "tables")
RESULTS_FIG_DIR <- file.path(PROJECT_ROOT, "results", "figures")
PA_DIR          <- file.path(RESULTS_TAB_DIR, "pathway_activity")

# Output to sandbox-safe /tmp/
OUT_TAB_DIR <- "/tmp/fc_output/tables"
OUT_FIG_DIR <- "/tmp/fc_output/figures"
FC_TAB_DIR   <- file.path(OUT_TAB_DIR, "fair_comparison")
FC_FIG_DIR   <- file.path(OUT_FIG_DIR, "fair_comparison")

suppressWarnings(dir.create(FC_TAB_DIR, showWarnings = FALSE, recursive = TRUE))
suppressWarnings(dir.create(FC_FIG_DIR, showWarnings = FALSE, recursive = TRUE))

set.seed(42)

cat("============================================================\n")
cat("Script 27: LASSO-Cox Fair Comparison (Gene vs. Pathway)\n")
cat("============================================================\n\n")

# ============================================================
# 0. Load packages
# ============================================================
cat("=== [0] Loading packages ===\n\n")
required_pkgs <- c("survival", "glmnet", "rms", "ggplot2", "dplyr", "tidyr")
for (pkg in required_pkgs) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  cat(sprintf("  %s v%s loaded\n", pkg, packageVersion(pkg)))
}
cat("\n")


# ============================================================
# 1. Load and prepare data
# ============================================================
cat("=== [1] Loading data ===\n\n")

# Clinical data
clin <- read.csv(file.path(RESULTS_TAB_DIR, "GSE39582_xelox_groups.csv"),
                  stringsAsFactors = FALSE)
rownames(clin) <- clin$sample_id
cat(sprintf("  Clinical: %d samples\n", nrow(clin)))

# Gene expression (pre-filtered top 500 by univariate Cox p-value)
gene_expr <- readRDS("/tmp/fc_output/GSE39582_top500_genes.rds")
cat(sprintf("  Gene expr (top 500 pre-filtered): %d genes x %d samples\n", nrow(gene_expr), ncol(gene_expr)))

# Pathway scores (44 pathways)
pw_scores <- readRDS(file.path(PA_DIR, "GSE39582_pathway_scores.rds"))
cat(sprintf("  Pathway scores: %d pathways x %d samples\n", nrow(pw_scores), ncol(pw_scores)))

# Common samples
common_samples <- intersect(clin$sample_id, colnames(gene_expr))
cat(sprintf("  Common samples: %d\n", length(common_samples)))

# Build analysis dataframe
cox_df <- clin[common_samples, c("sample_id", "rfs_event", "rfs_delay",
                                  "tumor_location")]
cox_df$rfs_event <- as.numeric(cox_df$rfs_event)
cox_df$rfs_delay <- as.numeric(cox_df$rfs_delay)
cox_df <- cox_df[!is.na(cox_df$rfs_event) & !is.na(cox_df$rfs_delay) & cox_df$rfs_delay > 0, ]
valid_samples <- cox_df$sample_id

cat(sprintf("  Valid samples: %d\n", nrow(cox_df)))
cat(sprintf("  RFS events: %d (%.1f%%)\n", sum(cox_df$rfs_event), mean(cox_df$rfs_event) * 100))

# Prepare matrices
gene_mat <- t(gene_expr[, valid_samples])      # samples x genes
pw_mat   <- t(pw_scores[, valid_samples])       # samples x pathways
surv_obj <- Surv(cox_df$rfs_delay, cox_df$rfs_event)

cat(sprintf("  Gene matrix: %d x %d\n", nrow(gene_mat), ncol(gene_mat)))
cat(sprintf("  Pathway matrix: %d x %d\n", nrow(pw_mat), ncol(pw_mat)))


# ============================================================
# 2. Helper: LASSO-Cox with CV
# ============================================================
run_lasso_cox <- function(X, y, label, n_bootstrap = 1000) {
  # X = samples x features matrix
  # y = Surv object
  # Returns: list(selected, coefs, c_index, c_boot, c_corrected, epv)

  cat(sprintf("\n  [%s] Running LASSO-Cox...\n", label))

  # 2a. LASSO with cross-validation
  cv_fit <- cv.glmnet(X, y, family = "cox", alpha = 1, nfolds = 10)
  lambda_1se <- cv_fit$lambda.1se

  # Extract coefficients at lambda.1se
  coefs_raw <- coef(cv_fit, s = "lambda.1se")
  selected_idx <- which(coefs_raw[, 1] != 0)
  n_selected <- length(selected_idx)

  cat(sprintf("    lambda.1se = %.4f, selected = %d features\n", lambda_1se, n_selected))

  if (n_selected == 0) {
    cat("    WARNING: No features selected at lambda.1se! Using lambda.min.\n")
    coefs_raw <- coef(cv_fit, s = "lambda.min")
    selected_idx <- which(coefs_raw[, 1] != 0)
    n_selected <- length(selected_idx)
    cat(sprintf("    lambda.min selected = %d features\n", n_selected))
  }

  selected_names <- colnames(X)[selected_idx]
  selected_coefs <- coefs_raw[selected_idx, 1]

  coef_df <- data.frame(
    Feature = selected_names,
    Coefficient = selected_coefs,
    stringsAsFactors = FALSE
  )
  coef_df <- coef_df[order(abs(coef_df$Coefficient), decreasing = TRUE), ]

  # 2b. Apparent C-index (linear predictor from LASSO)
  lp <- predict(cv_fit, newx = X, s = "lambda.1se")[, 1]
  c_app <- as.numeric(1 - rcorr.cens(lp, y)["C Index"])
  epv <- sum(y[, 2]) / n_selected

  cat(sprintf("    Apparent C-index = %.4f, EPV = %.1f\n", c_app, epv))

  # 2b2. Early exit if 0 features — no bootstrap possible
  if (n_selected == 0) {
    cat("    SKIPPING bootstrap: 0 features selected.\n")
    return(list(
      label = label, n_features = 0, selected_features = coef_df,
      c_apparent = c_app, c_bootstrap = c(NA), optimism = NA,
      c_corrected = c_app, epv = epv, lambda_1se = lambda_1se,
      linear_predictor = lp, coef_stability = data.frame(),
      n_boot_effective = 0
    ))
  }

  # 2c. Bootstrap optimism correction
  c_boot <- numeric(n_bootstrap)
  for (b in seq_len(n_bootstrap)) {
    boot_idx <- sample(nrow(X), replace = TRUE)
    X_boot <- X[boot_idx, ]
    y_boot <- y[boot_idx]

    cv_boot <- tryCatch({
      cv.glmnet(X_boot, y_boot, family = "cox", alpha = 1, nfolds = 10)
    }, error = function(e) NULL)

    if (is.null(cv_boot)) {
      c_boot[b] <- NA
      next
    }

    # Performance on bootstrap sample (apparent)
    lp_boot <- tryCatch({
      predict(cv_boot, newx = X_boot, s = "lambda.1se")[, 1]
    }, error = function(e) rep(0, nrow(X_boot)))

    c_boot_app <- tryCatch({
      as.numeric(1 - rcorr.cens(lp_boot, y_boot)["C Index"])
    }, error = function(e) NA)

    # Performance on original sample (test)
    lp_orig <- tryCatch({
      predict(cv_boot, newx = X, s = "lambda.1se")[, 1]
    }, error = function(e) rep(0, nrow(X)))

    c_boot_test <- tryCatch({
      as.numeric(1 - rcorr.cens(lp_orig, y)["C Index"])
    }, error = function(e) NA)

    if (!is.na(c_boot_app) && !is.na(c_boot_test)) {
      c_boot[b] <- c_boot_app - c_boot_test  # optimism
    } else {
      c_boot[b] <- NA
    }
  }

  optimism <- mean(c_boot, na.rm = TRUE)
  c_corrected <- c_app - optimism

  # 2d. Bootstrap coefficient stability
  boot_coefs <- matrix(NA, nrow = n_bootstrap, ncol = n_selected)
  colnames(boot_coefs) <- selected_names

  for (b in seq_len(n_bootstrap)) {
    boot_idx <- sample(nrow(X), replace = TRUE)
    X_boot <- X[boot_idx, ]
    y_boot <- y[boot_idx]

    cv_boot <- tryCatch({
      cv.glmnet(X_boot, y_boot, family = "cox", alpha = 1, nfolds = 10)
    }, error = function(e) NULL)

    if (!is.null(cv_boot)) {
      coefs_b <- coef(cv_boot, s = "lambda.1se")
      for (j in seq_along(selected_names)) {
        idx <- which(rownames(coefs_b) == selected_names[j])
        if (length(idx) > 0) boot_coefs[b, j] <- coefs_b[idx, 1]
      }
    }
  }

  # Summary statistics for each coefficient
  coef_stability <- data.frame(
    Feature = selected_names,
    Median = apply(boot_coefs, 2, median, na.rm = TRUE),
    Q025 = apply(boot_coefs, 2, quantile, 0.025, na.rm = TRUE),
    Q975 = apply(boot_coefs, 2, quantile, 0.975, na.rm = TRUE),
    PropNonZero = colMeans(!is.na(boot_coefs) & boot_coefs != 0),
    stringsAsFactors = FALSE
  )

  cat(sprintf("    Bootstrap-corrected C-index = %.4f (optimism = %.4f)\n",
              c_corrected, optimism))
  cat(sprintf("    Effective bootstrap reps: %d / %d\n",
              sum(!is.na(c_boot)), n_bootstrap))

  return(list(
    label = label,
    n_features = n_selected,
    selected_features = coef_df,
    c_apparent = c_app,
    c_bootstrap = c_boot,
    optimism = optimism,
    c_corrected = c_corrected,
    epv = epv,
    lambda_1se = lambda_1se,
    linear_predictor = lp,
    coef_stability = coef_stability,
    n_boot_effective = sum(!is.na(c_boot))
  ))
}


# ============================================================
# 3. Run LASSO-Cox on gene-level data
# ============================================================
cat("\n=== [2] Gene-level LASSO-Cox ===\n")
gene_result <- run_lasso_cox(gene_mat, surv_obj, "Gene-level", n_bootstrap = 200)


# ============================================================
# 4. Run LASSO-Cox on pathway-level data
# ============================================================
cat("\n=== [3] Pathway-level LASSO-Cox ===\n")
pw_result <- run_lasso_cox(pw_mat, surv_obj, "Pathway-level", n_bootstrap = 1000)


# ============================================================
# 5. Build comparison summary
# ============================================================
cat("\n=== [4] Comparison Summary ===\n\n")

comparison <- data.frame(
  Metric = c(
    "Feature type", "Modeling method",
    "Total candidate features", "Features selected",
    "Apparent C-index", "Bootstrap optimism",
    "Bootstrap-corrected C-index", "Events", "EPV"
  ),
  Gene_Level = c(
    "Individual genes", "LASSO-Cox (cv.glmnet)",
    ncol(gene_mat), gene_result$n_features,
    sprintf("%.4f", gene_result$c_apparent),
    sprintf("%.4f", gene_result$optimism),
    sprintf("%.4f", gene_result$c_corrected),
    sum(cox_df$rfs_event),
    sprintf("%.1f", gene_result$epv)
  ),
  Pathway_Level = c(
    "ssGSEA pathway scores", "LASSO-Cox (cv.glmnet)",
    ncol(pw_mat), pw_result$n_features,
    sprintf("%.4f", pw_result$c_apparent),
    sprintf("%.4f", pw_result$optimism),
    sprintf("%.4f", pw_result$c_corrected),
    sum(cox_df$rfs_event),
    sprintf("%.1f", pw_result$epv)
  ),
  stringsAsFactors = FALSE
)

cat("  Comparison Summary:\n")
for (i in seq_len(nrow(comparison))) {
  cat(sprintf("  %-30s | %-20s | %-20s\n",
              comparison$Metric[i], comparison$Gene_Level[i], comparison$Pathway_Level[i]))
}


# ============================================================
# 6. Save results
# ============================================================
cat("\n=== [5] Saving results ===\n\n")

# 6a. Gene-level results
write.csv(gene_result$selected_features,
          file.path(FC_TAB_DIR, "gene_lasso_cox_results.csv"),
          row.names = FALSE)
write.csv(gene_result$coef_stability,
          file.path(FC_TAB_DIR, "gene_lasso_coef_stability.csv"),
          row.names = FALSE)

# 6b. Pathway-level results
write.csv(pw_result$selected_features,
          file.path(FC_TAB_DIR, "pathway_lasso_cox_results.csv"),
          row.names = FALSE)
write.csv(pw_result$coef_stability,
          file.path(FC_TAB_DIR, "pathway_lasco_coef_stability.csv"),
          row.names = FALSE)

# 6c. Comparison summary
write.csv(comparison,
          file.path(FC_TAB_DIR, "fair_comparison_summary.csv"),
          row.names = FALSE)

cat("  Saved:\n")
cat(sprintf("    %s/gene_lasso_cox_results.csv\n", FC_TAB_DIR))
cat(sprintf("    %s/pathway_lasso_cox_results.csv\n", FC_TAB_DIR))
cat(sprintf("    %s/fair_comparison_summary.csv\n", FC_TAB_DIR))


# ============================================================
# 7. Visualizations
# ============================================================
cat("\n=== [6] Generating visualizations ===\n\n")

# 7a. Side-by-side C-index bar chart
c_idx_df <- data.frame(
  Model = factor(c("Gene-level LASSO-Cox", "Gene-level LASSO-Cox",
                   "Pathway-level LASSO-Cox", "Pathway-level LASSO-Cox"),
                 levels = c("Gene-level LASSO-Cox", "Pathway-level LASSO-Cox")),
  Type = c("Apparent", "Bootstrap-corrected", "Apparent", "Bootstrap-corrected"),
  C_Index = c(gene_result$c_apparent, gene_result$c_corrected,
              pw_result$c_apparent, pw_result$c_corrected)
)

pdf(file.path(FC_FIG_DIR, "fair_comparison_c_index.pdf"), width = 6, height = 5)
p <- ggplot(c_idx_df, aes(x = Model, y = C_Index, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", C_Index)),
            position = position_dodge(width = 0.8), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Apparent" = "#69b3a2", "Bootstrap-corrected" = "#404080")) +
  labs(title = "Fair Comparison: LASSO-Cox C-index (Gene vs. Pathway)",
       subtitle = paste0("GSE39582 XELOX cohort, n=", nrow(cox_df),
                         ", events=", sum(cox_df$rfs_event)),
       y = "C-index", x = "") +
  ylim(0, max(c_idx_df$C_Index) * 1.15) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p)
dev.off()
cat(sprintf("  -> %s\n", file.path(FC_FIG_DIR, "fair_comparison_c_index.pdf")))

# 7b. Gene-level coefficient forest plot
if (gene_result$n_features > 1) {
  plot_df_gene <- gene_result$coef_stability
  plot_df_gene$Feature <- factor(plot_df_gene$Feature,
                                  levels = plot_df_gene$Feature[order(plot_df_gene$Median)])

  pdf(file.path(FC_FIG_DIR, "gene_lasso_coef_forest.pdf"), width = 8,
      height = max(3, gene_result$n_features * 0.4))
  p <- ggplot(plot_df_gene, aes(x = Median, y = Feature)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(size = 2.5, color = "#2196F3") +
    geom_errorbarh(aes(xmin = Q025, xmax = Q975), height = 0.2, color = "#2196F3") +
    labs(title = paste0("Gene-level LASSO-Cox Coefficients (", gene_result$n_features, " genes)"),
         subtitle = "Bootstrap median, 2.5th, and 97.5th percentiles (B=1000)",
         x = "Log Hazard Ratio", y = "") +
    theme_minimal()
  print(p)
  dev.off()
  cat(sprintf("  -> %s\n", file.path(FC_FIG_DIR, "gene_lasso_coef_forest.pdf")))
}

# 7c. Pathway-level coefficient forest plot
if (pw_result$n_features > 1) {
  plot_df_pw <- pw_result$coef_stability
  plot_df_pw$Feature <- factor(plot_df_pw$Feature,
                                levels = plot_df_pw$Feature[order(plot_df_pw$Median)])

  pdf(file.path(FC_FIG_DIR, "pathway_lasso_coef_forest.pdf"), width = 8,
      height = max(3, pw_result$n_features * 0.5))
  p <- ggplot(plot_df_pw, aes(x = Median, y = Feature)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(size = 2.5, color = "#4CAF50") +
    geom_errorbarh(aes(xmin = Q025, xmax = Q975), height = 0.2, color = "#4CAF50") +
    labs(title = paste0("Pathway-level LASSO-Cox Coefficients (", pw_result$n_features, " pathways)"),
         subtitle = "Bootstrap median, 2.5th, and 97.5th percentiles (B=1000)",
         x = "Log Hazard Ratio", y = "") +
    theme_minimal()
  print(p)
  dev.off()
  cat(sprintf("  -> %s\n", file.path(FC_FIG_DIR, "pathway_lasso_coef_forest.pdf")))
}

# 7d. Cross-validated lambda plot (gene)
pdf(file.path(FC_FIG_DIR, "gene_lasso_lambda.pdf"), width = 7, height = 5)
cv_gene <- cv.glmnet(gene_mat, surv_obj, family = "cox", alpha = 1, nfolds = 10)
plot(cv_gene)
title("Gene-level LASSO: Cross-validated Lambda")
dev.off()
cat(sprintf("  -> %s\n", file.path(FC_FIG_DIR, "gene_lasso_lambda.pdf")))

# 7e. Cross-validated lambda plot (pathway)
pdf(file.path(FC_FIG_DIR, "pathway_lasso_lambda.pdf"), width = 7, height = 5)
cv_pw <- cv.glmnet(pw_mat, surv_obj, family = "cox", alpha = 1, nfolds = 10)
plot(cv_pw)
title("Pathway-level LASSO: Cross-validated Lambda")
dev.off()
cat(sprintf("  -> %s\n", file.path(FC_FIG_DIR, "pathway_lasso_lambda.pdf")))


# ============================================================
# 8. Final summary
# ============================================================
cat("\n============================================================\n")
cat("Script 27 completed successfully.\n")
cat("============================================================\n")
cat(sprintf("  Gene-level:    C = %.4f (corrected %.4f), EPV = %.1f, %d genes\n",
            gene_result$c_apparent, gene_result$c_corrected,
            gene_result$epv, gene_result$n_features))
cat(sprintf("  Pathway-level: C = %.4f (corrected %.4f), EPV = %.1f, %d pathways\n",
            pw_result$c_apparent, pw_result$c_corrected,
            pw_result$epv, pw_result$n_features))
cat(sprintf("  Delta C (pathway - gene) = %.4f\n",
            pw_result$c_corrected - gene_result$c_corrected))
