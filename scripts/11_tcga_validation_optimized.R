# ============================================================
# 11_TCGA_VALIDATION_OPTIMIZED.R
# XELOX Resistance - TCGA Validation Optimization
# Multi-strategy validation for imbalanced external data
# ============================================================
# Purpose:
#   Strategy A: Weighted Fingerprint Risk Score (SHAP-weighted z-score)
#   Strategy B: Multi-level gene set scoring (8/27/144 genes)
#   Strategy C: Directional concordance analysis (binomial test)
#   Strategy D: Permutation test for AUC significance (1000 iters)
#
# Root cause: current direct XGBoost prediction gives AUC=0.5375
#   due to platform effects (microarray -> RNA-seq) + extreme imbalance (5R/48S)
#
# Output:
#   tables/ml_phase2/TCGA_validation_optimized_report.txt
#   tables/ml_phase2/risk_scores_comparison.csv
#   figures/ml_phase2/Validation_risk_score_boxplot.pdf
#   figures/ml_phase2/Validation_multi_level_auc.pdf
#   figures/ml_phase2/Validation_permutation_test.pdf
# ============================================================

Sys.setenv(TMPDIR = "/tmp", TMP = "/tmp", TEMP = "/tmp")
Sys.setlocale("LC_ALL", "C")
.libPaths(c("/path/to/Rlibs", .libPaths()))

PROJECT_ROOT     <- "/path/to/xelox_project"
DATA_TCGA_DIR    <- file.path(PROJECT_ROOT, "data", "tcga")
RESULTS_TAB_DIR  <- file.path(PROJECT_ROOT, "results", "tables")
RESULTS_FIG_DIR  <- file.path(PROJECT_ROOT, "results", "figures")
ML2_TAB_DIR      <- file.path(RESULTS_TAB_DIR, "ml_phase2")
ML2_FIG_DIR      <- file.path(RESULTS_FIG_DIR, "ml_phase2")

# Create output subdirectory
OUT_DIR <- file.path(ML2_TAB_DIR, "validation_optimized")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

cat("============================================================\n")
cat("XELOX Resistance - TCGA Validation Optimization\n")
cat("Multi-Strategy External Validation\n")
cat("============================================================\n\n")

# ============================================================
# 0. Load required packages
# ============================================================
cat("=== [0] Loading packages ===\n\n")
required_pkgs <- c("pROC", "ggplot2", "data.table")
for (pkg in required_pkgs) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  cat(sprintf("  %s v%s loaded\n", pkg, packageVersion(pkg)))
}
cat("\n")

# ============================================================
# 1. Load Phase 2 outputs
# ============================================================
cat("=== [1] Loading Phase 2 outputs ===\n\n")

# 1a. SHAP summary (top fingerprint genes with importance weights)
shap_file <- file.path(ML2_TAB_DIR, "SHAP_summary_phase2.csv")
shap_df <- read.csv(shap_file, stringsAsFactors = FALSE)
rownames(shap_df) <- shap_df$Gene
cat(sprintf("  SHAP summary: %d genes\n", nrow(shap_df)))

# 1b. Multi-algorithm voting data (all 173 voting genes with vote counts)
votes_file <- file.path(ML2_TAB_DIR, "MultiAlgo_feature_votes.csv")
votes_df <- read.csv(votes_file, stringsAsFactors = FALSE)
cat(sprintf("  Voting data: %d genes\n", nrow(votes_df)))

# 1c. Fingerprint genes
fingerprint_file <- file.path(ML2_TAB_DIR, "XELOX_resistance_fingerprint.csv")
fingerprint_df <- read.csv(fingerprint_file, stringsAsFactors = FALSE)
fp_genes <- fingerprint_df$Gene
cat(sprintf("  Fingerprint genes: %d\n", length(fp_genes)))

# ============================================================
# 2. Load GSE39582 direction info from DEG data
# ============================================================
cat("=== [2] Loading GSE39582 direction information ===\n\n")

deg_file <- file.path(RESULTS_TAB_DIR, "DEG_GSE39582_mapped.csv")
deg_df <- read.csv(deg_file, stringsAsFactors = FALSE)

# Build gene -> logFC and direction map
# Take the most significant probe per gene (smallest P.Value)
deg_clean <- deg_df[!is.na(deg_df$gene_symbol) & deg_df$gene_symbol != "", ]
gene_order <- order(deg_clean$P.Value)
deg_clean <- deg_clean[gene_order, ]
deg_dedup <- deg_clean[!duplicated(deg_clean$gene_symbol), ]
gene_direction <- setNames(deg_dedup$logFC, deg_dedup$gene_symbol)

# Helper: resolve a multi-gene entry (e.g. "GGT1 /// GGT2 /// GGTLC1")
# Returns the first gene in the entry that has a known direction
resolve_direction <- function(gene_entry, direction_map) {
  parts <- trimws(strsplit(gene_entry, " /// ")[[1]])
  for (p in parts) {
    p_upper <- toupper(p)
    # Try exact match
    if (p_upper %in% names(direction_map)) return(direction_map[p_upper])
    # Try mixed case match
    name_match <- names(direction_map)[toupper(names(direction_map)) == p_upper]
    if (length(name_match) > 0) return(direction_map[name_match[1]])
  }
  return(NA)
}

# Apply to fingerprint genes
cat("  Fingerprint gene directions (logFC in GSE39582, resistant vs sensitive):\n")
fp_directions <- sapply(fp_genes, function(g) resolve_direction(g, gene_direction))
names(fp_directions) <- fp_genes

for (i in seq_along(fp_genes)) {
  cat(sprintf("    %s: logFC = %s\n", fp_genes[i], 
              ifelse(is.na(fp_directions[i]), "NA", sprintf("%+.4f", fp_directions[i]))))
}
cat("\n")

# ============================================================
# 3. Load TCGA expression data
# ============================================================
cat("=== [3] Loading TCGA expression data ===\n\n")

tcga_file <- file.path(DATA_TCGA_DIR, "TCGA_COAD_READ_gene_pool_expression.csv")
tcga_full <- read.csv(tcga_file, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("  TCGA data: %d samples x %d columns\n", nrow(tcga_full), ncol(tcga_full)))

# Filter to samples with response labels
tcga_resp <- tcga_full[tcga_full$response %in% c("resistant", "sensitive"), ]
cat(sprintf("  Labeled samples: %d\n", nrow(tcga_resp)))
cat(sprintf("    Resistant: %d, Sensitive: %d\n", 
            sum(tcga_resp$response == "resistant"), 
            sum(tcga_resp$response == "sensitive")))

# Binary response: 1 = resistant, 0 = sensitive
y_tcga <- ifelse(tcga_resp$response == "resistant", 1, 0)

# Metadata columns
meta_cols <- c("patient_barcode", "sample_barcode", "response", "cancer_type")
expr_cols <- setdiff(colnames(tcga_full), meta_cols)
cat(sprintf("  Expression columns: %d\n", length(expr_cols)))

# ============================================================
# 4. Strategy A: Weighted Fingerprint Risk Score
# ============================================================
cat("\n=== [4] Strategy A: Weighted Fingerprint Risk Score ===\n\n")

# Helper: match a multi-gene entry to TCGA columns (returns first matching column)
match_gene_to_tcga <- function(gene_entry, tcga_cols) {
  parts <- trimws(strsplit(gene_entry, " /// ")[[1]])
  for (p in parts) {
    p_upper <- toupper(p)
    match_idx <- which(toupper(tcga_cols) == p_upper)
    if (length(match_idx) > 0) return(tcga_cols[match_idx[1]])
  }
  return(NA)
}

# Map each fingerprint gene to TCGA columns
fp_to_tcga <- sapply(fp_genes, function(g) match_gene_to_tcga(g, expr_cols))
fp_found <- fp_to_tcga[!is.na(fp_to_tcga)]
cat(sprintf("  Fingerprint genes found in TCGA: %d / %d\n", length(fp_found), length(fp_genes)))

# Build risk score for each available gene
# risk_score = Σ( SHAP_weight * sign(GSE39582_direction) * z_score )
risk_gene_list <- names(fp_found)
risk_weights <- numeric(length(risk_gene_list))
risk_signs <- numeric(length(risk_gene_list))
gene_z_scores <- matrix(NA, nrow(tcga_resp), length(risk_gene_list))

for (i in seq_along(risk_gene_list)) {
  gene <- risk_gene_list[i]
  tcga_col <- fp_found[gene]
  
  # SHAP weight
  if (gene %in% rownames(shap_df)) {
    risk_weights[i] <- shap_df[gene, "MeanAbsSHAP"]
  } else {
    # For multi-gene entries not directly in SHAP df, use average of available
    risk_weights[i] <- mean(shap_df$MeanAbsSHAP, na.rm = TRUE)
  }
  
  # Direction sign from GSE39582
  dir_val <- fp_directions[gene]
  if (!is.na(dir_val)) {
    risk_signs[i] <- sign(dir_val)
  } else {
    risk_signs[i] <- 0  # neutral if unknown
  }
  
  # Z-score normalize TCGA expression (skip zero-variance genes)
  expr_vals <- as.numeric(tcga_resp[[tcga_col]])
  sd_val <- sd(expr_vals, na.rm = TRUE)
  if (is.na(sd_val) || sd_val == 0) {
    cat(sprintf("  %s: SKIPPED (zero variance in TCGA)\n", gene))
    risk_weights[i] <- 0  # effectively exclude
    gene_z_scores[, i] <- 0
    next
  }
  z_vals <- (expr_vals - mean(expr_vals, na.rm = TRUE)) / sd_val
  gene_z_scores[, i] <- z_vals
  
  cat(sprintf("  %s: weight=%.4f, sign=%+d, z-range=[%.2f, %.2f]\n",
              gene, risk_weights[i], risk_signs[i], min(z_vals), max(z_vals)))
}

# Compute composite risk score
risk_score_A <- rowSums(gene_z_scores * matrix(risk_weights * risk_signs, 
                                                nrow = nrow(tcga_resp), 
                                                ncol = length(risk_gene_list), 
                                                byrow = TRUE), na.rm = TRUE)
# Final NaN/Inf check
risk_score_A[is.nan(risk_score_A) | is.infinite(risk_score_A)] <- 0
cat(sprintf("\n  Risk Score A statistics:\n"))
cat(sprintf("    Range: [%.4f, %.4f]\n", min(risk_score_A), max(risk_score_A)))
cat(sprintf("    Mean (Resistant): %.4f\n", mean(risk_score_A[y_tcga == 1])))
cat(sprintf("    Mean (Sensitive): %.4f\n", mean(risk_score_A[y_tcga == 0])))

# ROC analysis for risk score A
if (length(unique(y_tcga)) < 2 || sd(risk_score_A) == 0) {
  cat(sprintf("\n  SKIPPING ROC: no variation in risk score or only one class\n"))
  auc_A <- NA; ci_A <- c(NA, NA, NA)
  wt_A <- list(p.value = NA)
} else {
  roc_A <- roc(y_tcga, risk_score_A, quiet = TRUE)
  auc_A <- as.numeric(roc_A$auc)
  ci_A <- ci.auc(roc_A)
  cat(sprintf("\n  ROC Analysis:\n"))
  cat(sprintf("    AUC: %.4f (95%% CI: %.4f - %.4f)\n", auc_A, ci_A[1], ci_A[3]))

  # Wilcoxon rank-sum test
  wt_A <- wilcox.test(risk_score_A ~ y_tcga, alternative = "two.sided")
  cat(sprintf("    Wilcoxon p-value: %.6f\n", wt_A$p.value))
}

# Store results (outside else so risk_results always exists)
risk_results <- data.frame(
  Strategy = "A_SHAPweighted_Top8",
  GeneSet = "8 fingerprint genes (SHAP-weighted)",
  N_Genes = length(risk_gene_list),
  AUC = ifelse(is.na(auc_A), NA, auc_A),
  AUC_CI_lower = ifelse(is.na(ci_A[1]), NA, ci_A[1]),
  AUC_CI_upper = ifelse(is.na(ci_A[3]), NA, ci_A[3]),
  Wilcoxon_P = ifelse(is.null(wt_A$p.value), NA, wt_A$p.value),
  Mean_Risk_Resistant = mean(risk_score_A[y_tcga == 1]),
  Mean_Risk_Sensitive = mean(risk_score_A[y_tcga == 0]),
  stringsAsFactors = FALSE
)

# ============================================================
# 5. Strategy B: Multi-level gene set scoring
# ============================================================
cat("\n=== [5] Strategy B: Multi-level gene set scoring ===\n\n")

# Define gene sets with vote-count weights
# Set 1: Top 8 fingerprint (same as Strategy A, for comparison)
# Set 2: ≥3 voting genes (27 found in TCGA)
# Set 3: All voting genes (144 found in TCGA)

build_gene_set <- function(vote_df, min_votes, tcga_cols, direction_map, shap_weights) {
  # Filter genes by minimum vote
  candidate_genes <- vote_df$Gene[vote_df$Votes >= min_votes]
  
  # For each gene, find TCGA column
  result_list <- list()
  for (g in candidate_genes) {
    tcga_col <- match_gene_to_tcga(g, tcga_cols)
    if (is.na(tcga_col)) next
    
    # Direction from GSE39582 DEG data
    dir_val <- resolve_direction(g, direction_map)
    if (is.na(dir_val)) next  # skip genes without direction info
    
    # Weight: vote count
    weight <- vote_df$Votes[vote_df$Gene == g]
    
    # Sign: direction of effect
    sign_val <- sign(dir_val)
    
    result_list[[length(result_list) + 1]] <- list(
      gene = g,
      tcga_col = tcga_col,
      weight = weight,
      sign = sign_val
    )
  }
  return(result_list)
}

# Process each gene set
gene_sets <- list(
  list(name = "B_Vote3plus", label = "≥3 Votes (27 TCGA genes)", min_votes = 3),
  list(name = "C_AllVoting", label = "All 173 Voting (144 TCGA genes)", min_votes = 1)
)

for (gs in gene_sets) {
  cat(sprintf("  Processing: %s (min_votes >= %d)\n", gs$name, gs$min_votes))
  
  # Build gene set
  gene_set <- build_gene_set(votes_df, gs$min_votes, expr_cols, gene_direction, shap_df)
  
  if (length(gene_set) < 5) {
    cat(sprintf("    Too few genes (%d), skipping\n", length(gene_set)))
    next
  }
  cat(sprintf("    Found %d genes in TCGA with direction info\n", length(gene_set)))
  
  # Compute weighted risk score (skip genes with zero variance -> NaN z-scores)
  risk_score <- rep(0, nrow(tcga_resp))
  n_used <- 0
  effective_weight <- 0
  for (gs_gene in gene_set) {
    expr_vals <- as.numeric(tcga_resp[[gs_gene$tcga_col]])
    sd_val <- sd(expr_vals, na.rm = TRUE)
    if (is.na(sd_val) || sd_val == 0) {
      cat(sprintf("    Skipping %s (zero variance in TCGA)\n", gs_gene$gene))
      next
    }
    z_vals <- (expr_vals - mean(expr_vals, na.rm = TRUE)) / sd_val
    risk_score <- risk_score + gs_gene$weight * gs_gene$sign * z_vals
    n_used <- n_used + 1
    effective_weight <- effective_weight + gs_gene$weight
  }
  # Normalize by effective weight
  if (effective_weight > 0) {
    risk_score <- risk_score / effective_weight
  }
  cat(sprintf("    Genes effectively used: %d / %d\n", n_used, length(gene_set)))
  
  cat(sprintf("    Normalized risk score range: [%.4f, %.4f]\n", min(risk_score), max(risk_score)))
  cat(sprintf("    Mean (Resistant): %.4f, Mean (Sensitive): %.4f\n",
              mean(risk_score[y_tcga == 1]), mean(risk_score[y_tcga == 0])))
  
  # Clean NaN/Inf values before ROC
  risk_score[is.nan(risk_score) | is.infinite(risk_score)] <- 0
  
  # Check that we have both classes
  if (length(unique(y_tcga)) < 2 || all(risk_score == risk_score[1])) {
    cat(sprintf("    SKIPPING: no variation in risk score or only one class present\n"))
    next
  }
  
  # ROC
  roc_obj <- roc(y_tcga, risk_score, quiet = TRUE)
  auc_val <- as.numeric(roc_obj$auc)
  ci_obj <- ci.auc(roc_obj)
  cat(sprintf("    AUC: %.4f (95%% CI: %.4f - %.4f)\n", auc_val, ci_obj[1], ci_obj[3]))
  
  # Wilcoxon
  wt <- wilcox.test(risk_score ~ y_tcga, alternative = "two.sided")
  cat(sprintf("    Wilcoxon p-value: %.6f\n", wt$p.value))
  
  # Store
  risk_results <- rbind(risk_results, data.frame(
    Strategy = gs$name,
    GeneSet = gs$label,
    N_Genes = length(gene_set),
    AUC = auc_val,
    AUC_CI_lower = ci_obj[1],
    AUC_CI_upper = ci_obj[3],
    Wilcoxon_P = wt$p.value,
    Mean_Risk_Resistant = mean(risk_score[y_tcga == 1]),
    Mean_Risk_Sensitive = mean(risk_score[y_tcga == 0]),
    stringsAsFactors = FALSE
  ))
}

cat("\n")

# ============================================================
# 6. Strategy C: Directional Concordance
# ============================================================
cat("=== [6] Strategy C: Directional Concordance ===\n\n")

# For each available fingerprint gene, compute the mean expression
# in resistant vs sensitive TCGA samples and check if the direction
# matches GSE39582

concordance_results <- data.frame(
  Gene = character(),
  GSE39582_logFC = numeric(),
  GSE39582_Direction = character(),
  TCGA_Mean_Resistant = numeric(),
  TCGA_Mean_Sensitive = numeric(),
  TCGA_logFC = numeric(),
  TCGA_Direction = character(),
  Concordant = logical(),
  stringsAsFactors = FALSE
)

n_concordant <- 0
n_total_concordance <- 0

for (gene in names(fp_found)) {
  tcga_col <- fp_found[gene]
  
  # GSE39582 direction
  dir_val <- fp_directions[gene]
  if (is.na(dir_val)) next
  gse_dir <- ifelse(dir_val > 0, "up", "down")
  
  # TCGA: mean expression in each group
  resistant_vals <- as.numeric(tcga_resp[tcga_resp$response == "resistant", tcga_col])
  sensitive_vals <- as.numeric(tcga_resp[tcga_resp$response == "sensitive", tcga_col])
  
  mean_res <- mean(resistant_vals, na.rm = TRUE)
  mean_sens <- mean(sensitive_vals, na.rm = TRUE)
  tcga_logFC <- mean_res - mean_sens
  tcga_dir <- ifelse(tcga_logFC > 0, "up", "down")
  
  # Check concordance
  concordant <- (gse_dir == tcga_dir)
  if (concordant) n_concordant <- n_concordant + 1
  n_total_concordance <- n_total_concordance + 1
  
  concordance_results <- rbind(concordance_results, data.frame(
    Gene = gene,
    GSE39582_logFC = round(dir_val, 4),
    GSE39582_Direction = gse_dir,
    TCGA_Mean_Resistant = round(mean_res, 4),
    TCGA_Mean_Sensitive = round(mean_sens, 4),
    TCGA_logFC = round(tcga_logFC, 4),
    TCGA_Direction = tcga_dir,
    Concordant = concordant,
    stringsAsFactors = FALSE
  ))
}

# Display concordance
cat("  Directional concordance per gene:\n")
print(concordance_results, row.names = FALSE)

# Binomial test
binom_p <- binom.test(n_concordant, n_total_concordance, p = 0.5, alternative = "greater")$p.value
cat(sprintf("\n  Concordance: %d / %d (%.1f%%)\n", n_concordant, n_total_concordance,
            100 * n_concordant / n_total_concordance))
cat(sprintf("  Binomial test (one-sided): p = %.4f\n", binom_p))

# ============================================================
# 7. Strategy D: Permutation Test (1000 iterations)
# ============================================================
cat("\n=== [7] Strategy D: Permutation Test ===\n\n")

if (is.na(auc_A) || sd(risk_score_A) == 0) {
  cat("  SKIPPING: risk score A is not valid for permutation test\n")
  n_perm <- 0
  perm_aucs <- numeric(0)
  empirical_p <- NA
} else {
  n_perm <- 1000
  perm_aucs <- numeric(n_perm)
  perm_scores <- risk_score_A  # Use Strategy A risk score

cat(sprintf("  Running %d permutations...\n", n_perm))
pb <- txtProgressBar(min = 0, max = n_perm, style = 3)

for (i in 1:n_perm) {
  # Shuffle labels
  y_perm <- sample(y_tcga)
  # Compute AUC
  perm_roc <- roc(y_perm, perm_scores, quiet = TRUE)
  perm_aucs[i] <- as.numeric(perm_roc$auc)
  setTxtProgressBar(pb, i)
}
close(pb)

# Empirical p-value (proportion of permuted AUCs >= observed AUC)
empirical_p <- mean(perm_aucs >= auc_A)
cat(sprintf("\n  Observed AUC: %.4f\n", auc_A))
cat(sprintf("  Permutation mean AUC: %.4f\n", mean(perm_aucs)))
cat(sprintf("  Permutation SD: %.4f\n", sd(perm_aucs)))
cat(sprintf("  Empirical p-value: %.4f (%d/%d >= %.4f)\n", 
            empirical_p, sum(perm_aucs >= auc_A), n_perm, auc_A))
}  # end else for permutation test

# ============================================================
# 8. Figures
# ============================================================
cat("\n=== [8] Generating figures ===\n\n")

# 8a. Risk score boxplot (Strategy A)
tcga_plot_df <- data.frame(
  risk_score = risk_score_A,
  response = tcga_resp$response,
  stringsAsFactors = FALSE
)

p1 <- ggplot(tcga_plot_df, aes(x = response, y = risk_score, fill = response)) +
  geom_boxplot(alpha = 0.7, outlier.size = 2, outlier.alpha = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2.5) +
  scale_fill_manual(values = c("resistant" = "#E41A1C", "sensitive" = "#377EB8")) +
  labs(
    title = "TCGA Validation: SHAP-Weighted Risk Score",
    subtitle = sprintf("AUC = %.4f | Wilcoxon p = %.4f | %d samples (R=%d, S=%d)",
                       auc_A, wt_A$p.value, nrow(tcga_resp),
                       sum(y_tcga == 1), sum(y_tcga == 0)),
    x = "XELOX Response", y = "Fingerprint Risk Score"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(ML2_FIG_DIR, "Validation_risk_score_boxplot.pdf"), p1, 
       width = 6, height = 7)
cat("  Risk score boxplot saved\n")

# 8b. Multi-level AUC comparison (bar chart)
auc_plot_df <- risk_results
auc_plot_df$AUC_Display <- sprintf("%.4f", auc_plot_df$AUC)
auc_plot_df$Label <- paste0("n=", auc_plot_df$N_Genes)

p2 <- ggplot(auc_plot_df, aes(x = GeneSet, y = AUC, fill = GeneSet)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.6) +
  geom_errorbar(aes(ymin = AUC_CI_lower, ymax = AUC_CI_upper), width = 0.15, size = 0.8) +
  geom_text(aes(label = AUC_Display), vjust = -0.8, size = 3.8, fontface = "bold") +
  geom_text(aes(label = Label), vjust = -2.5, size = 3.2, color = "grey40") +
  scale_fill_brewer(palette = "Set2") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "TCGA Validation: Multi-Level Gene Set Comparison",
    subtitle = sprintf("AUC with 95%% DeLong CI | Permutation p = %.4f", empirical_p),
    x = "Gene Set Strategy", y = "AUC"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 20, hjust = 0.8)
  )

ggsave(file.path(ML2_FIG_DIR, "Validation_multi_level_auc.pdf"), p2, 
       width = 7, height = 6)
cat("  Multi-level AUC comparison saved\n")

# 8c. Permutation test distribution
perm_plot_df <- data.frame(auc = perm_aucs)

p3 <- ggplot(perm_plot_df, aes(x = auc)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.7) +
  geom_vline(xintercept = auc_A, color = "#E41A1C", linetype = "dashed", size = 1.2) +
  annotate("text", x = auc_A + 0.03, y = max(table(cut(perm_aucs, 30))) * 0.9,
           label = sprintf("Observed\nAUC = %.4f\np = %.4f", auc_A, empirical_p),
           color = "#E41A1C", size = 4, hjust = 0, fontface = "bold") +
  labs(
    title = "Permutation Test: AUC Null Distribution",
    subtitle = sprintf("%d permutations | Mean AUC = %.4f | SD = %.4f",
                       n_perm, mean(perm_aucs), sd(perm_aucs)),
    x = "AUC under Random Labels", y = "Frequency"
  ) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(ML2_FIG_DIR, "Validation_permutation_test.pdf"), p3, 
       width = 7, height = 6)
cat("  Permutation test distribution saved\n")

# ============================================================
# 9. Save results
# ============================================================
cat("\n=== [9] Saving results ===\n\n")

# Save risk scores comparison
write.csv(risk_results, file.path(OUT_DIR, "risk_scores_comparison.csv"), row.names = FALSE)
cat("  Risk scores comparison saved\n")

# Save risk scores for individual samples
sample_risk_df <- data.frame(
  sample_barcode = tcga_resp$sample_barcode,
  response = tcga_resp$response,
  cancer_type = tcga_resp$cancer_type,
  fingerprint_risk_score = risk_score_A,
  stringsAsFactors = FALSE
)
write.csv(sample_risk_df, file.path(OUT_DIR, "TCGA_risk_scores_per_sample.csv"), row.names = FALSE)
cat("  Per-sample risk scores saved\n")

# Save concordance results
write.csv(concordance_results, file.path(OUT_DIR, "directional_concordance.csv"), row.names = FALSE)
cat("  Directional concordance saved\n")

# Save permutation results
write.csv(data.frame(permuted_AUC = perm_aucs), 
          file.path(OUT_DIR, "permutation_test_results.csv"), row.names = FALSE)
cat("  Permutation results saved\n")

# ============================================================
# 10. Comprehensive Report
# ============================================================
cat("\n=== [10] Writing comprehensive report ===\n\n")

report <- c(
  "==============================================================================",
  "XELOX Resistance - TCGA Validation Optimization Report",
  paste("Date:", Sys.Date()),
  "==============================================================================",
  "",
  "--- Background ---",
  sprintf("  Original direct XGBoost prediction gave AUC = 0.5375"),
  sprintf("  Root cause: platform effects (microarray -> RNA-seq) + extreme imbalance"),
  "",
  "--- TCGA Validation Data ---",
  sprintf("  Total samples: %d (Resistant=%d, Sensitive=%d)",
          nrow(tcga_resp), sum(y_tcga == 1), sum(y_tcga == 0)),
  sprintf("  Cancer types: COAD=%d, READ=%d",
          sum(tcga_resp$cancer_type == "COAD"),
          sum(tcga_resp$cancer_type == "READ")),
  sprintf("  Gene set coverage: %d/%d fingerprint genes found in TCGA",
          length(fp_found), length(fp_genes)),
  "",
  "--- Strategy A: SHAP-Weighted Risk Score (Primary) ---",
  sprintf("  Genes: %d fingerprint genes with SHAP importance weighting",
          length(risk_gene_list)),
  sprintf("  Method: z-score normalization + SHAP weight * direction sign"),
  sprintf("  AUC: %.4f (95%% CI: %.4f - %.4f)", auc_A, ci_A[1], ci_A[3]),
  sprintf("  Wilcoxon p-value: %.6f", wt_A$p.value),
  sprintf("  Mean risk score (Resistant): %.4f", mean(risk_score_A[y_tcga == 1])),
  sprintf("  Mean risk score (Sensitive): %.4f", mean(risk_score_A[y_tcga == 0]))
)

# Add Strategy B results
report <- c(report, "",
  "--- Strategy B: Multi-Level Gene Set Scoring ---")
for (i in seq_len(nrow(risk_results))) {
  report <- c(report,
    sprintf("  %s: AUC=%.4f (CI: %.4f-%.4f), p=%.6f, n_genes=%d",
            risk_results$Strategy[i], risk_results$AUC[i],
            risk_results$AUC_CI_lower[i], risk_results$AUC_CI_upper[i],
            risk_results$Wilcoxon_P[i], risk_results$N_Genes[i]))
}

report <- c(report, "",
  "--- Strategy C: Directional Concordance ---",
  sprintf("  Concordant genes: %d / %d (%.1f%%)", 
          n_concordant, n_total_concordance, 
          100 * n_concordant / n_total_concordance),
  sprintf("  Binomial test p-value: %.4f", binom_p),
  "")

for (i in seq_len(nrow(concordance_results))) {
  report <- c(report,
    sprintf("  %s: GSE39582=%s, TCGA=%s -> %s",
            concordance_results$Gene[i],
            concordance_results$GSE39582_Direction[i],
            concordance_results$TCGA_Direction[i],
            ifelse(concordance_results$Concordant[i], "CONCORDANT", "DISCORDANT")))
}

report <- c(report, "",
  "--- Strategy D: Permutation Test (1000 iterations) ---",
  sprintf("  Observed AUC: %.4f", auc_A),
  sprintf("  Permutation mean AUC: %.4f (SD = %.4f)", mean(perm_aucs), sd(perm_aucs)),
  sprintf("  Empirical p-value: %.4f", empirical_p),
  sprintf("  Interpretation: %s",
          ifelse(empirical_p < 0.05, 
                 "AUC is statistically significant (p < 0.05)",
                 "AUC is NOT statistically significant (p >= 0.05)")),
  "",
  "--- Summary ---",
  sprintf("  Best AUC: %.4f (%s)", 
          max(risk_results$AUC), 
          risk_results$Strategy[which.max(risk_results$AUC)]),
  sprintf("  Directional concordance: %.1f%% of genes match between datasets",
          100 * n_concordant / n_total_concordance),
  sprintf("  Permutation significance: %s",
          ifelse(empirical_p < 0.05, "SIGNIFICANT", "NOT SIGNIFICANT")),
  "",
  "--- Output Files ---",
  sprintf("  Risk scores: %s", file.path(OUT_DIR, "risk_scores_comparison.csv")),
  sprintf("  Per-sample scores: %s", file.path(OUT_DIR, "TCGA_risk_scores_per_sample.csv")),
  sprintf("  Concordance: %s", file.path(OUT_DIR, "directional_concordance.csv")),
  sprintf("  Permutation: %s", file.path(OUT_DIR, "permutation_test_results.csv")),
  sprintf("  Figures: %s", ML2_FIG_DIR),
  "",
  "=============================================================================="
)

writeLines(report, file.path(OUT_DIR, "TCGA_validation_optimized_report.txt"))
cat(paste(report, collapse = "\n"), "\n")

cat("\n=== TCGA Validation Optimization Complete! ===\n")
cat(sprintf("  Report: %s\n", file.path(OUT_DIR, "TCGA_validation_optimized_report.txt")))
cat(sprintf("  Figures: %s\n", ML2_FIG_DIR))
