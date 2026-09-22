# ============================================================
# 27b_BOOTSTRAP_FIX.R
# Re-run bootstrap stability analysis for merged training set
#
# Load existing merged_expression_combat.rds and redo bootstrap
# with B=500 (more stable than 1000) and robust error handling.
# ============================================================

Sys.setenv(TMPDIR = "C:/temp", TMP = "C:/temp", TEMP = "C:/temp")
.libPaths(c("C:/Rlibs", .libPaths()))

PROJECT_ROOT <- "C:/xelox_research"
DATA_GEO_DIR <- file.path(PROJECT_ROOT, "data", "geo")
RESULTS_TAB_DIR <- file.path(PROJECT_ROOT, "results", "tables")
PA_DIR <- file.path(RESULTS_TAB_DIR, "pathway_activity")
OUT_TAB_DIR <- file.path(RESULTS_TAB_DIR, "prs_merged")
OUT_FIG_DIR <- file.path(PROJECT_ROOT, "results", "figures", "prs_merged")

set.seed(42)

cat("============================================================\n")
cat("27b: Bootstrap Stability Fix for Merged Training\n")
cat("============================================================\n\n")

suppressPackageStartupMessages({
  library(sva)
  library(GSVA)
  library(GSEABase)
  library(glmnet)
  library(pROC)
  library(ggplot2)
})
cat("  Packages loaded\n\n")

# ============================================================
# 1. Load pre-computed merged data
# ============================================================
cat("=== [1] Loading merged data ===\n\n")

# Load merged ComBat expression
merged_combat <- readRDS(file.path(OUT_TAB_DIR, "merged_expression_combat.rds"))
cat(sprintf("  Merged expression: %d genes x %d samples\n",
            nrow(merged_combat), ncol(merged_combat)))

# Load sample info
sample_info <- read.csv(file.path(OUT_TAB_DIR, "merged_sample_info.csv"),
                         stringsAsFactors = FALSE)
merged_y <- sample_info$Response
names(merged_y) <- sample_info$Sample
batch <- sample_info$Batch

cat(sprintf("  Samples: %d, Resistant=%d, Sensitive=%d\n",
            length(merged_y), sum(merged_y == 1), sum(merged_y == 0)))
cat(sprintf("  Batches: GSE39582=%d, GSE19860=%d\n",
            sum(batch == "GSE39582"), sum(batch == "GSE19860")))

# ============================================================
# 2. Recompute ssGSEA (or load pre-computed)
# ============================================================
cat("=== [2] Loading/recomputing ssGSEA ===\n\n")

gmt_h <- "C:/temp/h.all.v2024.1.symbols.gmt"
gmt_kegg <- "C:/temp/c2.cp.kegg_legacy.v2024.1.symbols.gmt"

read_gmt <- function(gmt_file) {
  lines <- readLines(gmt_file, warn = FALSE)
  gs_list <- list()
  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    gs_name <- parts[1]
    genes <- parts[-(1:2)]
    genes <- genes[nchar(genes) > 0]
    gs_list[[gs_name]] <- genes
  }
  return(gs_list)
}

hallmark_all <- read_gmt(gmt_h)
kegg_all <- read_gmt(gmt_kegg)

target_pathways <- list()
kegg_targets <- c(
  "KEGG_NUCLEOTIDE_EXCISION_REPAIR", "KEGG_BASE_EXCISION_REPAIR",
  "KEGG_MISMATCH_REPAIR", "KEGG_HOMOLOGOUS_RECOMBINATION",
  "KEGG_P53_SIGNALING_PATHWAY", "KEGG_APOPTOSIS", "KEGG_CELL_CYCLE",
  "KEGG_ABC_TRANSPORTERS", "KEGG_GLUTATHIONE_METABOLISM",
  "KEGG_DRUG_METABOLISM_CYTOCHROME_P450", "KEGG_DRUG_METABOLISM_OTHER_ENZYMES",
  "KEGG_METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450",
  "KEGG_WNT_SIGNALING_PATHWAY", "KEGG_TGF_BETA_SIGNALING_PATHWAY",
  "KEGG_NOTCH_SIGNALING_PATHWAY", "KEGG_FOCAL_ADHESION",
  "KEGG_ECM_RECEPTOR_INTERACTION",
  "KEGG_MAPK_SIGNALING_PATHWAY", "KEGG_MTOR_SIGNALING_PATHWAY",
  "KEGG_JAK_STAT_SIGNALING_PATHWAY",
  "KEGG_CHEMOKINE_SIGNALING_PATHWAY", "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY", "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_COLORECTAL_CANCER", "KEGG_PATHWAYS_IN_CANCER"
)
hallmark_targets <- c(
  "HALLMARK_DNA_REPAIR", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_APOPTOSIS", "HALLMARK_ANGIOGENESIS", "HALLMARK_HYPOXIA",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_GLYCOLYSIS",
  "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_P53_PATHWAY",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING", "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2", "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MTORC1_SIGNALING", "HALLMARK_TGF_BETA_SIGNALING"
)

for (id in kegg_targets) {
  if (id %in% names(kegg_all)) target_pathways[[id]] <- kegg_all[[id]]
}
for (id in hallmark_targets) {
  if (id %in% names(hallmark_all)) target_pathways[[id]] <- hallmark_all[[id]]
}
cat(sprintf("  Target pathways: %d\n\n", length(target_pathways)))

# Try loading pre-computed pathway scores, otherwise recompute
ps_file <- file.path(OUT_TAB_DIR, "merged_pathway_scores.rds")
if (file.exists(ps_file)) {
  merged_es <- readRDS(ps_file)
  cat(sprintf("  Loaded pre-computed ssGSEA: %d x %d\n",
              nrow(merged_es), ncol(merged_es)))
} else {
  cat("  Computing ssGSEA...\n")
  gs_list <- lapply(names(target_pathways), function(nm) {
    GeneSet(setName = nm, geneIds = target_pathways[[nm]])
  })
  gs_collection <- GeneSetCollection(gs_list)
  merged_es <- gsva(ssgseaParam(merged_combat, gs_collection,
                                 minSize = 10, maxSize = 500, verbose = FALSE))
  cat(sprintf("  ssGSEA complete: %d x %d\n", nrow(merged_es), ncol(merged_es)))
}

# ============================================================
# 3. Prepare for bootstrap (same as Script 27 LASSO section)
# ============================================================
cat("=== [3] Preparing LASSO input ===\n\n")

train_ps <- merged_es
train_y <- merged_y

common_samples <- intersect(colnames(train_ps), names(train_y))
train_ps <- train_ps[, common_samples, drop = FALSE]
train_y <- train_y[common_samples]

train_x <- t(train_ps)
pathway_names <- colnames(train_x)

cat(sprintf("  Training: %d pathways x %d samples\n",
            nrow(train_ps), ncol(train_ps)))
cat(sprintf("  Outcome: Resistant=%d, Sensitive=%d\n\n",
            sum(train_y == 1), sum(train_y == 0)))

# Univariate AUC pre-filter
compute_auc <- function(x, y) {
  tryCatch({
    roc_obj <- roc(y, x, direction = "<", quiet = TRUE)
    as.numeric(auc(roc_obj))
  }, error = function(e) NA)
}

auc_vec <- apply(train_ps, 1, function(pw) compute_auc(pw, train_y))
names(auc_vec) <- pathway_names
auc_filter <- abs(auc_vec - 0.5) > 0.05
auc_filter[is.na(auc_filter)] <- FALSE

cat(sprintf("  AUC filter passed: %d / %d pathways\n",
            sum(auc_filter), length(auc_filter)))

train_x_filtered <- train_x[, auc_filter, drop = FALSE]
train_x_scaled <- scale(train_x_filtered)
n_features <- ncol(train_x_filtered)
cat(sprintf("  Filtered X: %d samples x %d features\n\n",
            nrow(train_x_filtered), n_features))

# ============================================================
# 4. Bootstrap stability analysis (B=500, robust)
# ============================================================
cat("=== [4] Bootstrap stability (B=500) ===\n\n")

B <- 500
boot_coef_matrix <- matrix(0, nrow = B, ncol = n_features)
colnames(boot_coef_matrix) <- colnames(train_x_filtered)
boot_selected_count <- 0
boot_lambda_min_count <- 0

set.seed(42)
for (b in 1:B) {
  # Stratified bootstrap
  idx_res <- which(train_y == 1)
  idx_sens <- which(train_y == 0)
  boot_idx <- c(
    sample(idx_res, size = length(idx_res), replace = TRUE),
    sample(idx_sens, size = length(idx_sens), replace = TRUE)
  )
  
  x_boot <- train_x_scaled[boot_idx, , drop = FALSE]
  y_boot <- train_y[boot_idx]
  
  cv_boot <- tryCatch(
    cv.glmnet(x = x_boot, y = y_boot, family = "binomial",
               alpha = 1, nfolds = 5, type.measure = "deviance"),
    error = function(e) NULL
  )
  
  if (is.null(cv_boot)) next
  
  # Try lambda.1se first
  coef_boot <- as.matrix(coef(cv_boot, s = "lambda.1se"))
  nz_1se <- sum(coef_boot[-1, 1] != 0)
  
  if (nz_1se == 0) {
    coef_boot <- as.matrix(coef(cv_boot, s = "lambda.min"))
    nz_min <- sum(coef_boot[-1, 1] != 0)
    if (nz_min > 0) boot_lambda_min_count <- boot_lambda_min_count + 1
  }
  
  boot_coef_matrix[b, ] <- coef_boot[-1, 1]
  if (sum(coef_boot[-1, 1] != 0) > 0) boot_selected_count <- boot_selected_count + 1
  
  if (b %% 100 == 0) cat(sprintf("  Bootstrap %d / %d\n", b, B))
}

cat(sprintf("\n  Bootstraps with >=1 selected: %d / %d\n",
            boot_selected_count, B))
cat(sprintf("  lambda.min fallbacks: %d / %d\n", boot_lambda_min_count, B))

# ============================================================
# 5. Selection frequency
# ============================================================
cat("\n=== [5] Computing stability metrics ===\n\n")

selection_freq <- colMeans(boot_coef_matrix != 0)
freq_df <- data.frame(
  Pathway = names(selection_freq),
  SelectionFreq = round(selection_freq, 4),
  stringsAsFactors = FALSE
)
freq_df <- freq_df[order(freq_df$SelectionFreq, decreasing = TRUE), ]

cat("  Top 15 most selected pathways:\n")
for (i in 1:min(15, nrow(freq_df))) {
  cat(sprintf("    %2d. %s: %.1f%%\n", i, freq_df$Pathway[i],
              freq_df$SelectionFreq[i] * 100))
}

cat("\n  Pathways with >10% selection frequency:\n")
freq_10 <- freq_df[freq_df$SelectionFreq > 0.1, ]
for (i in seq_len(nrow(freq_10))) {
  cat(sprintf("    %s: %.1f%%\n", freq_10$Pathway[i], freq_10$SelectionFreq[i] * 100))
}

# ============================================================
# 6. Pairwise Jaccard stability
# ============================================================
cat("\n  Computing pairwise Jaccard...\n")

compute_pairwise_jaccard <- function(coef_mat, n_sample = 200) {
  selected <- coef_mat != 0
  n_boot <- nrow(selected)
  if (n_boot < 2) return(NA)
  
  if (n_boot > n_sample) {
    idx <- sample(n_boot, n_sample)
    selected <- selected[idx, ]
    n_boot <- n_sample
  }
  
  jaccard_vals <- numeric(n_boot * (n_boot - 1) / 2)
  idx <- 1
  for (i in 1:(n_boot - 1)) {
    set_i <- which(selected[i, ])
    for (j in (i + 1):n_boot) {
      set_j <- which(selected[j, ])
      intersection <- length(intersect(set_i, set_j))
      union_len <- length(union(set_i, set_j))
      jaccard_vals[idx] <- ifelse(union_len > 0, intersection / union_len, 0)
      idx <- idx + 1
    }
  }
  return(mean(jaccard_vals, na.rm = TRUE))
}

valid_boot <- rowSums(boot_coef_matrix != 0) > 0
n_valid <- sum(valid_boot)

if (n_valid >= 10) {
  jaccard_stability <- compute_pairwise_jaccard(boot_coef_matrix[valid_boot, ])
  cat(sprintf("  Pairwise Jaccard stability: %.4f (n_valid=%d)\n",
              jaccard_stability, n_valid))
} else {
  jaccard_stability <- NA
  cat(sprintf("  Jaccard: NA (only %d valid bootstraps)\n", n_valid))
}

# ============================================================
# 7. Sign consistency and mean coefficients
# ============================================================
cat("\n  Computing sign consistency...\n")

mean_coef <- colMeans(boot_coef_matrix, na.rm = TRUE)
mean_abs_coef_nonzero <- apply(boot_coef_matrix, 2, function(col) {
  nz <- col[col != 0]
  if (length(nz) == 0) return(NA)
  return(mean(abs(nz), na.rm = TRUE))
})

sign_consistency <- apply(boot_coef_matrix, 2, function(col) {
  non_zero <- col != 0
  if (sum(non_zero) == 0) return(NA)
  pos_frac <- mean(col[non_zero] > 0, na.rm = TRUE)
  return(max(pos_frac, 1 - pos_frac))
})

cv_vec <- apply(boot_coef_matrix, 2, function(col) {
  nz <- col[col != 0]
  if (length(nz) < 2) return(NA)
  abs_mean <- abs(mean(nz))
  if (abs_mean < 1e-10) return(NA)
  sd(nz) / abs_mean
})

stability_df <- data.frame(
  Pathway = colnames(boot_coef_matrix),
  SelectionFreq = round(selection_freq, 4),
  MeanCoef = round(mean_coef, 6),
  MeanAbsCoefNonzero = round(mean_abs_coef_nonzero, 6),
  SignConsistency = round(sign_consistency, 4),
  CV = round(cv_vec, 4),
  stringsAsFactors = FALSE
)
stability_df <- stability_df[order(stability_df$SelectionFreq, decreasing = TRUE), ]

# Save
write.csv(stability_df, file.path(OUT_TAB_DIR, "merged_bootstrap_stability.csv"), row.names = FALSE)
cat("  Saved stability table\n")

# ============================================================
# 8. Write bootstrap report
# ============================================================
cat("\n=== [6] Writing bootstrap report ===\n\n")

sink(file.path(OUT_TAB_DIR, "merged_bootstrap_report.txt"))
cat("Bootstrap Stability Report (GSE19860 Merged Training, N=204)\n")
cat("==============================================================\n\n")
cat(sprintf("Bootstraps: B = %d\n", B))
cat(sprintf("Valid bootstraps (>=1 nonzero): %d / %d (%.1f%%)\n",
            n_valid, B, n_valid / B * 100))
cat(sprintf("lambda.min fallbacks: %d / %d (%.1f%%)\n\n",
            boot_lambda_min_count, B, boot_lambda_min_count / B * 100))
cat(sprintf("Pairwise Jaccard stability: %.4f\n\n", jaccard_stability))

cat(sprintf("Pathways with >50%% selection frequency: %d\n",
            sum(stability_df$SelectionFreq > 0.5, na.rm = TRUE)))
cat(sprintf("Pathways with >30%% selection frequency: %d\n",
            sum(stability_df$SelectionFreq > 0.3, na.rm = TRUE)))
cat(sprintf("Pathways with >10%% selection frequency: %d\n\n",
            sum(stability_df$SelectionFreq > 0.1, na.rm = TRUE)))

cat("Full pathway stability table:\n")
cat(sprintf("%-50s %8s %12s %16s %16s %10s\n",
            "Pathway", "Freq%", "MeanCoef", "Mean|Coef|NZ", "SignCons%", "CV"))
cat(paste(rep("-", 120), collapse = ""), "\n")
for (i in seq_len(nrow(stability_df))) {
  cat(sprintf("%-50s %7.1f%% %11.6f %15.6f %15.1f%% %9.2f\n",
              stability_df$Pathway[i],
              stability_df$SelectionFreq[i] * 100,
              stability_df$MeanCoef[i],
              stability_df$MeanAbsCoefNonzero[i],
              stability_df$SignConsistency[i] * 100,
              stability_df$CV[i]))
}

cat("\n\n--- Final Model (from Script 27 LASSO) ---\n")
final_coefs <- read.csv(file.path(OUT_TAB_DIR, "merged_prs_coefficients.csv"),
                          stringsAsFactors = FALSE)
if (nrow(final_coefs) > 0) {
  for (i in seq_len(nrow(final_coefs))) {
    pw <- final_coefs$Pathway[i]
    freq <- if (pw %in% stability_df$Pathway) {
      stability_df$SelectionFreq[stability_df$Pathway == pw] * 100
    } else NA
    cat(sprintf("  %s (coef=%.4f, bootstrap_freq=%.1f%%)\n",
                pw, final_coefs$Coefficient[i], freq))
  }
}

cat("\n\n--- Comparison with original training (Script 25) ---\n")
cat("  Original (GSE39582+GSE104645, N=277): Jaccard = 0.1002\n")
cat(sprintf("  New (GSE39582+GSE19860, N=204): Jaccard = %.4f\n", jaccard_stability))
if (!is.na(jaccard_stability)) {
  delta <- jaccard_stability - 0.1002
  cat(sprintf("  Delta = %+.4f\n", delta))
  if (jaccard_stability > 0.15) {
    cat("  ✓ MODERATE IMPROVEMENT\n")
  } else if (jaccard_stability > 0.10) {
    cat("  ✓ MINOR IMPROVEMENT\n")
  } else {
    cat("  ✗ NO IMPROVEMENT (still unstable)\n")
  }
}

cat("\n--- Key Pathways Check ---\n")
key_pws <- c("HALLMARK_NOTCH_SIGNALING", "KEGG_NOTCH_SIGNALING_PATHWAY",
             "HALLMARK_APOPTOSIS", "KEGG_APOPTOSIS",
             "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "KEGG_ECM_RECEPTOR_INTERACTION")
for (pw in key_pws) {
  if (pw %in% stability_df$Pathway) {
    freq <- stability_df$SelectionFreq[stability_df$Pathway == pw] * 100
    cat(sprintf("  %s: bootstrap freq = %.1f%%\n", pw, freq))
  }
}

sink()
cat("  Bootstrap report saved\n\n")

# ============================================================
# 9. Stability curve
# ============================================================
if (n_valid >= 20) {
  cat("  Generating stability curve...\n")
  selected_mat <- boot_coef_matrix[valid_boot, ] != 0
  n_valid_use <- min(n_valid, 500)
  
  cumulative_jaccard <- numeric(n_valid_use - 1)
  for (k in 2:n_valid_use) {
    set_a <- which(selected_mat[1, ])
    set_b <- which(colMeans(selected_mat[1:k, , drop = FALSE]) > 0.5)
    intersection <- length(intersect(set_a, set_b))
    union_len <- length(union(set_a, set_b))
    cumulative_jaccard[k - 1] <- ifelse(union_len > 0, intersection / union_len, 0)
  }
  
  tryCatch({
    pdf(file.path(OUT_FIG_DIR, "bootstrap_stability_curve_merged.pdf"),
        width = 7, height = 5.5)
    plot(2:n_valid_use, cumulative_jaccard, type = "l", lwd = 2,
         col = "#2166AC", xlab = "Bootstrap Replicates",
         ylab = "Cumulative Jaccard Index",
         main = "Stability Curve — Merged Training (N=204)",
         ylim = c(0, 1))
    abline(h = 0.5, lty = 2, col = "gray50")
    text(n_valid_use * 0.8, 0.52, "Jaccard=0.5", col = "gray50", cex = 0.8)
    dev.off()
    cat("  Stability curve saved\n")
  }, error = function(e) cat("  Stability curve failed:", e$message, "\n"))
}

# ============================================================
# Final summary
# ============================================================
cat("\n============================================================\n")
cat("27b Bootstrap Fix Complete\n")
cat("============================================================\n")
cat(sprintf("  Valid bootstraps: %d / %d\n", n_valid, B))
cat(sprintf("  Pairwise Jaccard: %.4f\n", jaccard_stability))
cat(sprintf("  Pathways >10%%: %d\n",
            sum(stability_df$SelectionFreq > 0.1, na.rm = TRUE)))
cat(sprintf("  Pathways >50%%: %d\n",
            sum(stability_df$SelectionFreq > 0.5, na.rm = TRUE)))
cat("============================================================\n")
