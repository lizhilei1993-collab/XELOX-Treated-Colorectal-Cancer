# ============================================================
# 33_GENE_LEVEL_AIC_COX.R
# Gene-level AIC Cox Model as Ablation Control
#
# Purpose:
#   Build a gene-level prognostic model using backward AIC selection
#   on the same GSE39582 XELOX cohort (n=227, 79 RFS events).
#   This serves as an ablation control to demonstrate that pathway-level
#   aggregation provides better prognostic value than individual genes.
#
# Method:
#   1. Load gene expression matrix for GSE39582 (2,537 genes)
#   2. Univariate Cox screening (p < 0.05)
#   3. Backward AIC selection on significant genes
#   4. Compare C-index with pathway-level 8-variable Nomogram
#
# Input:
#   results/tables/pathway_activity/GSE39582_gene_expression.rds
#   results/tables/GSE39582_xelox_groups.csv
#
# Output:
#   results/tables/gene_level/gene_cox_univariate.csv
#   results/tables/gene_level/gene_cox_aic_model.csv
#   results/tables/gene_level/ablation_comparison.csv
#   results/figures/gene_level/forest_plot_gene_cox.pdf
# ============================================================

Sys.setenv(TMPDIR = "/tmp", TMP = "/tmp", TEMP = "/tmp")
.libPaths(c("/path/to/Rlibs", .libPaths()))

PROJECT_ROOT <- Sys.getenv("XELOX_ROOT",
  "/path/to/xelox_project基于可解释性机器学习的XELOX耐药分子指纹研究")

RESULTS_TAB_DIR <- file.path(PROJECT_ROOT, "results", "tables")
RESULTS_FIG_DIR <- file.path(PROJECT_ROOT, "results", "figures")
PA_DIR          <- file.path(RESULTS_TAB_DIR, "pathway_activity")
GENE_TAB_DIR    <- file.path(RESULTS_TAB_DIR, "gene_level")
GENE_FIG_DIR    <- file.path(RESULTS_FIG_DIR, "gene_level")
NOMO_TAB_DIR    <- file.path(RESULTS_TAB_DIR, "nomogram")

suppressWarnings(dir.create(GENE_TAB_DIR, showWarnings = FALSE, recursive = TRUE))
suppressWarnings(dir.create(GENE_FIG_DIR, showWarnings = FALSE, recursive = TRUE))

set.seed(42)

cat("============================================================\n")
cat("Script 33: Gene-level AIC Cox Model (Ablation Control)\n")
cat("============================================================\n\n")

# ============================================================
# 0. Load required packages
# ============================================================
cat("=== [0] Loading packages ===\n\n")
required_pkgs <- c("survival", "survminer", "rms", "MASS", "ggplot2")
for (pkg in required_pkgs) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  cat(sprintf("  %s v%s loaded\n", pkg, packageVersion(pkg)))
}
cat("\n")

# ============================================================
# 1. Load data
# ============================================================
cat("=== [1] Loading data ===\n\n")

# Load clinical data
clin <- read.csv(file.path(RESULTS_TAB_DIR, "GSE39582_xelox_groups.csv"),
                  stringsAsFactors = FALSE)
rownames(clin) <- clin$sample_id
cat(sprintf("  Clinical data: %d samples\n", nrow(clin)))

# Load gene expression matrix
gene_expr <- readRDS(file.path(PA_DIR, "GSE39582_gene_expression.rds"))
cat(sprintf("  Gene expression: %d genes x %d samples\n", nrow(gene_expr), ncol(gene_expr)))

# Common samples
common_samples <- intersect(clin$sample_id, colnames(gene_expr))
cat(sprintf("  Common samples: %d\n", length(common_samples)))

# Build analysis dataframe
cox_df <- clin[common_samples, c("sample_id", "rfs_event", "rfs_delay",
                                  "sex", "age", "tnm_stage", "tumor_location",
                                  "mmr_status", "kras", "braf")]
cox_df$rfs_event <- as.numeric(cox_df$rfs_event)
cox_df$rfs_delay <- as.numeric(cox_df$rfs_delay)

# Remove missing survival data
cox_df <- cox_df[!is.na(cox_df$rfs_event) & !is.na(cox_df$rfs_delay) & cox_df$rfs_delay > 0, ]
valid_samples <- cox_df$sample_id

cat(sprintf("  Valid samples for analysis: %d\n", nrow(cox_df)))
cat(sprintf("  RFS events: %d (%.1f%%)\n", sum(cox_df$rfs_event), mean(cox_df$rfs_event) * 100))

# Subset gene expression to valid samples
gene_expr_sub <- t(gene_expr[, valid_samples])  # samples x genes
cat(sprintf("  Gene matrix for analysis: %d samples x %d genes\n", nrow(gene_expr_sub), ncol(gene_expr_sub)))

# ============================================================
# 2. Univariate Cox screening
# ============================================================
cat("\n=== [2] Univariate Cox screening ===\n\n")

surv_obj <- Surv(cox_df$rfs_delay, cox_df$rfs_event)

# Function to compute univariate Cox
compute_univariate_cox <- function(gene_name) {
  gene_vals <- gene_expr_sub[, gene_name]
  if (sd(gene_vals, na.rm = TRUE) < 1e-10) {
    return(data.frame(Gene = gene_name, HR = NA, HR_lower = NA, HR_upper = NA,
                      p_value = NA, stringsAsFactors = FALSE))
  }
  
  tryCatch({
    cox_fit <- coxph(surv_obj ~ gene_vals)
    coefs <- summary(cox_fit)$coefficients
    ci <- exp(confint(cox_fit))
    
    data.frame(
      Gene = gene_name,
      HR = round(exp(coefs[1, "coef"]), 3),
      HR_lower = round(ci[1, 1], 3),
      HR_upper = round(ci[1, 2], 3),
      p_value = round(coefs[1, "Pr(>|z|)"], 6),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(Gene = gene_name, HR = NA, HR_lower = NA, HR_upper = NA,
               p_value = NA, stringsAsFactors = FALSE)
  })
}

# Run univariate Cox for all genes
cat("  Running univariate Cox for all genes...\n")
all_genes <- colnames(gene_expr_sub)
cox_results <- do.call(rbind, lapply(all_genes, compute_univariate_cox))

# Count significant genes
sig_genes <- cox_results[!is.na(cox_results$p_value) & cox_results$p_value < 0.05, ]
cat(sprintf("  Total genes tested: %d\n", nrow(cox_results)))
cat(sprintf("  Genes with p < 0.05: %d (%.1f%%)\n", nrow(sig_genes), nrow(sig_genes)/nrow(cox_results)*100))

# Save univariate results
write.csv(cox_results, file.path(GENE_TAB_DIR, "gene_cox_univariate.csv"), row.names = FALSE)
cat(sprintf("  -> Saved: %s\n", file.path(GENE_TAB_DIR, "gene_cox_univariate.csv")))

# ============================================================
# 3. Backward AIC selection
# ============================================================
cat("\n=== [3] Backward AIC selection ===\n\n")

# Limit to top 50 genes to make AIC selection feasible
max_genes_for_aic <- 50
if (nrow(sig_genes) > max_genes_for_aic) {
  cat(sprintf("  Limiting to top %d genes by p-value (from %d significant)\n",
              max_genes_for_aic, nrow(sig_genes)))
  sig_genes <- head(sig_genes[order(sig_genes$p_value), ], max_genes_for_aic)
} else if (nrow(sig_genes) < 2) {
  cat("  WARNING: Too few significant genes for AIC selection\n")
  cat("  Using top 20 genes by p-value instead\n")
  sig_genes <- head(cox_results[order(cox_results$p_value), ], 20)
}

# Prepare data for multivariable model
sig_gene_names <- sig_genes$Gene
gene_matrix_sig <- gene_expr_sub[, sig_gene_names, drop = FALSE]

# Scale gene values
gene_matrix_scaled <- scale(gene_matrix_sig)

# Build initial multivariable Cox model
# Create a combined dataframe with survival + scaled gene data
cox_model_df <- as.data.frame(gene_matrix_scaled)
cox_model_df$rfs_delay <- cox_df$rfs_delay
cox_model_df$rfs_event <- cox_df$rfs_event

formula_str <- paste0("Surv(rfs_delay, rfs_event) ~ ",
                      paste0("`", sig_gene_names, "`", collapse = " + "))
cox_full <- coxph(as.formula(formula_str), data = cox_model_df)

cat(sprintf("  Full model: %d genes\n", length(sig_gene_names)))
cat(sprintf("  Full model AIC: %.2f\n", AIC(cox_full)))

# Backward AIC selection
cat("  Running backward AIC selection...\n")
cox_aic <- stepAIC(cox_full, direction = "backward", trace = 0, data = cox_model_df)

# Extract selected genes
selected_genes <- names(coef(cox_aic))
cat(sprintf("  Selected genes after AIC: %d\n", length(selected_genes)))

# Build results table
aic_coefs <- summary(cox_aic)$coefficients
aic_ci <- exp(confint(cox_aic))

aic_results <- data.frame(
  Gene = selected_genes,
  HR = round(exp(aic_coefs[, "coef"]), 3),
  HR_lower = round(aic_ci[, 1], 3),
  HR_upper = round(aic_ci[, 2], 3),
  p_value = round(aic_coefs[, "Pr(>|z|)"], 4),
  Coefficient = round(aic_coefs[, "coef"], 4),
  stringsAsFactors = FALSE
)

# Sort by p-value
aic_results <- aic_results[order(aic_results$p_value), ]

cat("\n  Selected genes:\n")
for (i in seq_len(nrow(aic_results))) {
  cat(sprintf("    %2d. %-20s HR=%.3f (%.3f-%.3f), p=%.4f\n",
              i, aic_results$Gene[i], aic_results$HR[i],
              aic_results$HR_lower[i], aic_results$HR_upper[i],
              aic_results$p_value[i]))
}

# Save AIC model results
write.csv(aic_results, file.path(GENE_TAB_DIR, "gene_cox_aic_model.csv"), row.names = FALSE)
cat(sprintf("\n  -> Saved: %s\n", file.path(GENE_TAB_DIR, "gene_cox_aic_model.csv")))

# ============================================================
# 4. Compute C-index for gene-level model
# ============================================================
cat("\n=== [4] Computing C-index ===\n\n")

# C-index for gene-level model
tryCatch({
  cindex_gene <- concordance(cox_aic)$concordance
  cat(sprintf("  Gene-level AIC model C-index: %.3f\n", cindex_gene))
}, error = function(e) {
  # Manual C-index calculation
  pred_lp <- predict(cox_aic, type = "lp")
  concordant <- 0
  total <- 0
  for (i in seq_len(nrow(cox_df))) {
    for (j in seq_len(nrow(cox_df))) {
      if (i < j && cox_df$rfs_delay[i] < cox_df$rfs_delay[j] && cox_df$rfs_event[i] == 1) {
        total <- total + 1
        if (pred_lp[i] > pred_lp[j]) concordant <- concordant + 1
      }
    }
  }
  cindex_gene <<- ifelse(total > 0, concordant / total, NA)
  cat(sprintf("  Gene-level AIC model C-index (manual): %.3f\n", cindex_gene))
})

# Load pathway-level C-index from existing nomogram
nomo_model_file <- file.path(NOMO_TAB_DIR, "cox_final_model.txt")
if (file.exists(nomo_model_file)) {
  nomo_text <- readLines(nomo_model_file)
  concordance_line <- grep("Concordance", nomo_text, value = TRUE)
  if (length(concordance_line) > 0) {
    cindex_pathway <- as.numeric(sub(".*Concordance=\\s*([0-9.]+).*", "\\1", concordance_line[1]))
  } else {
    cindex_pathway <- NA
  }
  cat(sprintf("  Pathway-level Nomogram C-index: %.3f\n", cindex_pathway))
} else {
  cat("  WARNING: Pathway-level model file not found\n")
  cindex_pathway <- NA
}

# ============================================================
# 5. Ablation comparison
# ============================================================
cat("\n=== [5] Ablation comparison ===\n\n")

comparison_df <- data.frame(
  Model = c("Pathway-level Nomogram (8 variables)", "Gene-level AIC Cox"),
  Features = c(8, length(selected_genes)),
  C_index = c(cindex_pathway, cindex_gene),
  stringsAsFactors = FALSE
)

cat("  Model Comparison:\n")
cat(sprintf("  %-40s C-index = %.3f\n", comparison_df$Model[1], comparison_df$C_index[1]))
cat(sprintf("  %-40s C-index = %.3f\n", comparison_df$Model[2], comparison_df$C_index[2]))

if (!is.na(cindex_pathway) && !is.na(cindex_gene)) {
  cindex_diff <- cindex_gene - cindex_pathway
  cat(sprintf("\n  C-index difference (Gene - Pathway): %+.3f\n", cindex_diff))
  if (cindex_diff < 0) {
    cat("  => Pathway-level model outperforms gene-level model\n")
  } else {
    cat("  => Gene-level model outperforms pathway-level model\n")
  }
}

# Save comparison
write.csv(comparison_df, file.path(GENE_TAB_DIR, "ablation_comparison.csv"), row.names = FALSE)
cat(sprintf("\n  -> Saved: %s\n", file.path(GENE_TAB_DIR, "ablation_comparison.csv")))

# ============================================================
# 6. Forest plot for gene-level model
# ============================================================
cat("\n=== [6] Generating forest plot ===\n\n")

pdf(file.path(GENE_FIG_DIR, "forest_plot_gene_cox.pdf"), width = 10, height = 6)
par(mar = c(4, 12, 3, 4))

n_genes <- nrow(aic_results)
plot(NA, xlim = c(0, max(aic_results$HR_upper, na.rm = TRUE) * 1.3),
     ylim = c(0.5, n_genes + 1),
     xlab = "Hazard Ratio (95% CI)", ylab = "", yaxt = "n",
     main = "Gene-level AIC Cox Model: Selected Genes")
axis(2, at = n_genes:1, labels = rev(aic_results$Gene), las = 2)
abline(v = 1, lty = 2, col = "gray")

for (i in seq_len(n_genes)) {
  idx <- n_genes - i + 1
  points(aic_results$HR[i], i, pch = 15, cex = 1.2, col = "darkred")
  segments(aic_results$HR_lower[i], i, aic_results$HR_upper[i], i, lwd = 2, col = "darkred")
  text(max(aic_results$HR_upper, na.rm = TRUE) * 1.1, i,
       sprintf("%.2f (%.2f-%.2f)", aic_results$HR[i], aic_results$HR_lower[i], aic_results$HR_upper[i]),
       cex = 0.8, adj = 0)
}

dev.off()
cat(sprintf("  -> Saved: %s\n", file.path(GENE_FIG_DIR, "forest_plot_gene_cox.pdf")))

# ============================================================
# 7. Summary
# ============================================================
cat("\n=== [7] Analysis Summary ===\n\n")

cat("Gene-level AIC Cox Model:\n")
cat(sprintf("  - Total genes tested: %d\n", nrow(cox_results)))
cat(sprintf("  - Genes with p < 0.05: %d\n", nrow(sig_genes)))
cat(sprintf("  - Genes selected by AIC: %d\n", length(selected_genes)))
cat(sprintf("  - C-index: %.3f\n", cindex_gene))
tryCatch(cat(sprintf("  - AIC: %.2f\n", AIC(cox_aic))), error = function(e) cat("  - AIC: N/A\n"))

cat("\nPathway-level Nomogram (from Script 29):\n")
cat(sprintf("  - Features: 8 (7 pathways + 1 clinical)\n"))
cat(sprintf("  - C-index: %.3f\n", cindex_pathway))

cat("\nAblation Conclusion:\n")
if (!is.na(cindex_pathway) && !is.na(cindex_gene)) {
  if (cindex_pathway > cindex_gene) {
    cat("  Pathway-level aggregation provides superior prognostic value\n")
    cat("  compared to individual gene selection.\n")
  } else {
    cat("  Gene-level model shows comparable or better performance.\n")
    cat("  Consider further investigation of pathway aggregation benefit.\n")
  }
}

cat("\n=== Script 33 complete ===\n\n")
