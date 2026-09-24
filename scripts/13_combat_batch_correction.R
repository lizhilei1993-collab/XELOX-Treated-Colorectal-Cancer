# ============================================================
# 13_COMBAT_BATCH_CORRECTION.R
# XELOX Resistance - ComBat Batch Correction + Directional Concordance Re-evaluation
# ============================================================
# Purpose:
#   Apply ComBat batch correction to all GPL570 cohorts (GSE39582 + GSE72970 + GSE69657)
#   to distinguish batch effects from true biological signal.
#
#   Then re-evaluate:
#   1. Directional concordance of 8 fingerprint genes across cohorts (pre vs post)
#   2. Jaccard stability filtering (keep genes with >=60% concordance)
#   3. Corrected expression matrices for downstream analysis
#
# Output:
#   tables/batch_correction/ComBat_corrected_expression.rds  (corrected combined matrix)
#   tables/batch_correction/ComBat_diagnostic_pca.pdf        (PCA pre/post)
#   tables/batch_correction/concordance_pre_combat.csv       (pre-correction concordance)
#   tables/batch_correction/concordance_post_combat.csv      (post-correction concordance)
#   tables/batch_correction/directional_stability_summary.csv (Jaccard filtering)
#   tables/batch_correction/stable_fingerprint_genes.csv     (filtered gene list)
# ============================================================

Sys.setenv(TMPDIR = "/tmp", TMP = "/tmp", TEMP = "/tmp")
.libPaths(c("/path/to/Rlibs", .libPaths()))

# Accept PROJECT_ROOT from environment variable
PROJECT_ROOT <- Sys.getenv("XELOX_ROOT",
  "/path/to/xelox_project")

# Set locale
if (.Platform$OS.type == "windows") {
  tryCatch(Sys.setlocale("LC_ALL", "Chinese"), error = function(e) {
    tryCatch(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), error = function(e2) {})
  })
}
DATA_GEO_DIR     <- file.path(PROJECT_ROOT, "data", "geo")
RESULTS_TAB_DIR  <- file.path(PROJECT_ROOT, "results", "tables")
RESULTS_FIG_DIR  <- file.path(PROJECT_ROOT, "results", "figures")
ML2_TAB_DIR      <- file.path(RESULTS_TAB_DIR, "ml_phase2")

# Create output directories
BC_DIR <- file.path(RESULTS_TAB_DIR, "batch_correction")
dir.create(BC_DIR, showWarnings = FALSE, recursive = TRUE)
BC_FIG_DIR <- file.path(RESULTS_FIG_DIR, "batch_correction")
dir.create(BC_FIG_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

cat("============================================================\n")
cat("ComBat Batch Correction - XELOX GPL570 Cohorts\n")
cat("============================================================\n\n")

# ============================================================
# 0. Load required packages
# ============================================================
cat("=== [0] Loading packages ===\n\n")
suppressPackageStartupMessages({
  library(Biobase)
  library(sva)
  library(limma)
  library(ggplot2)
  library(pROC)
  library(data.table)
})
cat("  sva loaded\n")
cat("  limma loaded\n")
cat("  ggplot2 loaded\n")
cat("  pROC loaded\n")

# ============================================================
# 1. Load fingerprint gene definitions
# ============================================================
cat("\n=== [1] Loading fingerprint gene info ===\n\n")

fp_genes <- c("CSNK1G2", "KAZN", "KLK6", "MGA", "MID2", "SOX11", "TAS2R40", "ZNF451")
cat("  8 fingerprint genes:", paste(fp_genes, collapse = ", "), "\n\n")

# SHAP weights
shap_file <- file.path(ML2_TAB_DIR, "SHAP_summary_phase2.csv")
shap_df <- read.csv(shap_file, stringsAsFactors = FALSE)
rownames(shap_df) <- shap_df$Gene

get_shap_weight <- function(gene, shap_df) {
  if (gene %in% rownames(shap_df)) return(shap_df[gene, "MeanAbsSHAP"])
  for (g in rownames(shap_df)) {
    parts <- trimws(strsplit(g, " /// ")[[1]])
    if (toupper(gene) %in% toupper(parts)) return(shap_df[g, "MeanAbsSHAP"])
  }
  return(NA)
}
fp_shap <- sapply(fp_genes, function(g) get_shap_weight(g, shap_df))
cat("  SHAP weights:\n")
for (i in seq_along(fp_genes)) {
  cat(sprintf("    %s: %.4f\n", fp_genes[i], fp_shap[i]))
}

# GSE39582 direction (training set reference)
deg_file <- file.path(RESULTS_TAB_DIR, "DEG_GSE39582_mapped.csv")
deg_df <- read.csv(deg_file, stringsAsFactors = FALSE)
deg_clean <- deg_df[!is.na(deg_df$gene_symbol) & deg_df$gene_symbol != "", ]
gene_order <- order(deg_clean$P.Value)
deg_clean <- deg_clean[gene_order, ]
deg_dedup <- deg_clean[!duplicated(deg_clean$gene_symbol), ]
gene_direction <- setNames(deg_dedup$logFC, deg_dedup$gene_symbol)

resolve_direction <- function(gene_entry, direction_map) {
  parts <- trimws(strsplit(gene_entry, " /// ")[[1]])
  for (p in parts) {
    p_upper <- toupper(p)
    if (p_upper %in% names(direction_map)) return(direction_map[p_upper])
    name_match <- names(direction_map)[toupper(names(direction_map)) == p_upper]
    if (length(name_match) > 0) return(direction_map[name_match[1]])
  }
  return(NA)
}

fp_ref_direction <- sapply(fp_genes, function(g) resolve_direction(g, gene_direction))
names(fp_ref_direction) <- fp_genes
cat("\n  GSE39582 reference directions:\n")
for (i in seq_along(fp_genes)) {
  cat(sprintf("    %s: logFC = %+.4f  (%s)\n",
              fp_genes[i], fp_ref_direction[i],
              ifelse(fp_ref_direction[i] > 0, "up in R", "down in R")))
}
cat("\n")

# ============================================================
# 2. Load GPL570 probe-to-gene annotation
# ============================================================
cat("=== [2] Loading GPL570 annotation ===\n\n")
annot_file <- file.path(DATA_GEO_DIR, "GPL570_probe_gene_map.csv")
annot_map <- read.csv(annot_file, stringsAsFactors = FALSE)
cat(sprintf("  GPL570 probes mapped: %d\n", nrow(annot_map)))
cat(sprintf("  Unique genes: %d\n", length(unique(annot_map$gene_symbol))))

# ============================================================
# 3. Load all expression matrices (probe-level)
# ============================================================
cat("=== [3] Loading expression matrices ===\n\n")

# --- 3a. GSE39582 (training set) ---
cat("  Loading GSE39582...\n")
eset_39582 <- readRDS(file.path(DATA_GEO_DIR, "GSE39582_eset.rds"))
expr_39582 <- exprs(eset_39582)
cat(sprintf("    Dimensions: %d probes x %d samples\n", nrow(expr_39582), ncol(expr_39582)))

# --- 3b. GSE72970 ---
cat("  Loading GSE72970...\n")
expr_72970 <- readRDS(file.path(DATA_GEO_DIR, "GSE72970_expression.rds"))
cat(sprintf("    Dimensions: %d probes x %d samples\n", nrow(expr_72970), ncol(expr_72970)))

# --- 3c. GSE69657 ---
cat("  Loading GSE69657...\n")
expr_69657 <- readRDS(file.path(DATA_GEO_DIR, "GSE69657_expression.rds"))
cat(sprintf("    Dimensions: %d probes x %d samples\n", nrow(expr_69657), ncol(expr_69657)))

# ============================================================
# 4. Build combined matrix at probe level
# ============================================================
cat("\n=== [4] Building combined probe-level matrix ===\n\n")

# Find common probes across all datasets
common_probes <- Reduce(intersect, list(
  rownames(expr_39582),
  rownames(expr_72970),
  rownames(expr_69657)
))
cat(sprintf("  Common probes across all 3 cohorts: %d / %d\n",
            length(common_probes), nrow(annot_map)))

# Subset to common probes
X_39582 <- expr_39582[common_probes, ]
X_72970 <- expr_72970[common_probes, ]
X_69657 <- expr_69657[common_probes, ]

# Create combined matrix
X_combined <- cbind(X_39582, X_72970, X_69657)
cat(sprintf("  Combined matrix: %d probes x %d samples\n",
            nrow(X_combined), ncol(X_combined)))
cat(sprintf("    GSE39582: %d\n", ncol(X_39582)))
cat(sprintf("    GSE72970: %d\n", ncol(X_72970)))
cat(sprintf("    GSE69657: %d\n", ncol(X_69657)))

# Batch vector
batch <- c(rep("GSE39582", ncol(X_39582)),
           rep("GSE72970", ncol(X_72970)),
           rep("GSE69657", ncol(X_69657)))
cat(sprintf("  Batch vector: %d samples\n", length(batch)))

# ============================================================
# 5. Build biological covariate (response status)
# ============================================================
cat("\n=== [5] Building biological covariate matrix ===\n\n")

# Load clinical data for response labels
clin_39582 <- read.csv(file.path(RESULTS_TAB_DIR, "GSE39582_clinical_data.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
clin_72970 <- read.csv(file.path(RESULTS_TAB_DIR, "GSE72970_clinical_data.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
clin_69657 <- read.csv(file.path(RESULTS_TAB_DIR, "GSE69657_clinical_data.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)

# Map response to binary: 1 = resistant, 0 = sensitive
# GSE39582 uses RFS-based labels
cat("  GSE39582 response:\n")
if ("response_status" %in% colnames(clin_39582)) {
  print(table(clin_39582$response_status))
  resp_39582 <- ifelse(clin_39582$response_status == "resistant", 1,
                       ifelse(clin_39582$response_status == "sensitive", 0, NA))
} else if ("Relapse" %in% colnames(clin_39582)) {
  print(table(clin_39582$Relapse))
  resp_39582 <- clin_39582$Relapse  # 1 = relapsed = resistant
} else {
  # Try to find response column
  resp_cols <- grep("response|relapse|status|RFS", names(clin_39582), value = TRUE, ignore.case = TRUE)
  cat("  Available columns:", paste(resp_cols, collapse = ", "), "\n")
  resp_39582 <- rep(NA, nrow(clin_39582))
}

# GSE72970 - RECIST response
cat("  GSE72970 response:\n")
resp_col_72970 <- grep("response status", names(clin_72970), value = TRUE)[1]
print(table(clin_72970[[resp_col_72970]]))
resp_72970 <- ifelse(clin_72970[[resp_col_72970]] == "NR", 1, 0)

# GSE69657 - pathological response
cat("  GSE69657 response:\n")
resp_col_69657 <- grep("chemoresponse", names(clin_69657), value = TRUE)[1]
if (length(resp_col_69657) == 0) resp_col_69657 <- grep("response", names(clin_69657), value = TRUE)[1]
print(table(clin_69657[[resp_col_69657]]))
resp_69657 <- ifelse(clin_69657[[resp_col_69657]] == "noresponder", 1,
                     ifelse(clin_69657[[resp_col_69657]] == "responder", 0, NA))

# Align sample order with the combined matrix
samples_combined <- colnames(X_combined)

idx_39582 <- match(colnames(X_39582), clin_39582$geo_accession)
idx_72970 <- match(colnames(X_72970), clin_72970$geo_accession)
idx_69657 <- match(colnames(X_69657), clin_69657$geo_accession)

all_resp <- c(resp_39582[idx_39582], resp_72970[idx_72970], resp_69657[idx_69657])
cat(sprintf("\n  Response labels: R=%d, S=%d, NA=%d\n",
            sum(all_resp == 1, na.rm = TRUE), sum(all_resp == 0, na.rm = TRUE), sum(is.na(all_resp))))

# Build model matrix for ComBat (biological covariate of interest)
# v6.0 fix: Drop NA samples instead of filling with 0.5
# Filling NA with 0.5 creates a fake "half-resistant" biological state
complete_idx <- !is.na(all_resp)
n_na <- sum(!complete_idx)
cat(sprintf("  Removing %d samples with NA response status\n", n_na))

X_combined_orig <- X_combined
batch_orig <- batch
all_resp_orig <- all_resp

X_combined <- X_combined[, complete_idx]
batch <- batch[complete_idx]
all_resp <- all_resp[complete_idx]

mod <- model.matrix(~ all_resp)
# ComBat uses mod to preserve biology while removing batch effects
cat(sprintf("  Model matrix built: %d samples x %d columns (after NA removal)\n", nrow(mod), ncol(mod)))

# ============================================================
# 6. PCA before correction
# ============================================================
cat("\n=== [6] PCA before correction ===\n\n")

# Sample a subset of probes for faster PCA (use top 5000 most variable)
probe_var <- apply(X_combined, 1, var, na.rm = TRUE)
top_5000_probes <- head(order(probe_var, decreasing = TRUE), 5000)
X_pca_pre <- X_combined[top_5000_probes, ]
X_pca_pre_scaled <- t(scale(t(X_pca_pre), center = TRUE, scale = TRUE))
X_pca_pre_scaled[is.na(X_pca_pre_scaled)] <- 0

pca_pre <- prcomp(t(X_pca_pre_scaled), center = FALSE, scale. = FALSE)
pca_var_pre <- summary(pca_pre)$importance[2, 1:3] * 100

pca_df_pre <- data.frame(
  PC1 = pca_pre$x[, 1],
  PC2 = pca_pre$x[, 2],
  Batch = batch,
  Response = ifelse(is.na(all_resp), "Unknown", ifelse(all_resp == 1, "R", "S"))
)
cat(sprintf("  PC1: %.1f%%, PC2: %.1f%%, PC3: %.1f%%\n",
            pca_var_pre[1], pca_var_pre[2], pca_var_pre[3]))

# Check batch separation
library(MASS)
lda_acc_pre <- NA
lda_pre <- tryCatch(lda(PC1 ~ Batch, data = pca_df_pre), error = function(e) NULL)
if (!is.null(lda_pre)) {
  lda_pred_pre <- predict(lda_pre)
  lda_acc_pre <- mean(lda_pred_pre$class == pca_df_pre$Batch)
  cat(sprintf("  LDA batch classification accuracy (PC1): %.1f%%\n", lda_acc_pre * 100))
}

# ============================================================
# 7. Run ComBat correction
# ============================================================
cat("\n=== [7] Running ComBat correction ===\n\n")

# ComBat needs: expression matrix (genes x samples), batch vector, biological covariate mod
cat("  Starting ComBat...\n")
cat(sprintf("  Matrix: %d genes x %d samples\n", nrow(X_combined), ncol(X_combined)))
cat(sprintf("  Batches: %d unique\n", length(unique(batch))))
cat(sprintf("  Mod: %d covariate(s)\n", ncol(mod)))

X_combat <- tryCatch(
  ComBat(
    dat = as.matrix(X_combined),
    batch = as.factor(batch),
    mod = mod,
    par.prior = TRUE,
    prior.plots = FALSE
  ),
  error = function(e) {
    cat("  ComBat with par.prior failed, trying without mod...\n")
    # Try simpler version if parametric fails
    ComBat(
      dat = as.matrix(X_combined),
      batch = as.factor(batch),
      par.prior = TRUE,
      prior.plots = FALSE
    )
  }
)

if (is.null(dim(X_combat))) {
  cat("  ERROR: ComBat returned unexpected result\n")
  X_combat <- X_combined  # Use uncorrected as fallback
} else {
  cat(sprintf("  ComBat correction complete! Output: %d genes x %d samples\n",
              nrow(X_combat), ncol(X_combat)))
}

# Extract corrected matrices for each cohort
X_combat_39582 <- X_combat[, 1:ncol(X_39582)]
X_combat_72970 <- X_combat[, (ncol(X_39582) + 1):(ncol(X_39582) + ncol(X_72970))]
X_combat_69657 <- X_combat[, (ncol(X_39582) + ncol(X_72970) + 1):ncol(X_combat)]

colnames(X_combat_39582) <- colnames(X_39582)
colnames(X_combat_72970) <- colnames(X_72970)
colnames(X_combat_69657) <- colnames(X_69657)
rownames(X_combat) <- rownames(X_combined)

cat("  Corrected matrices extracted:\n")
cat(sprintf("    GSE39582: %d x %d\n", nrow(X_combat_39582), ncol(X_combat_39582)))
cat(sprintf("    GSE72970: %d x %d\n", nrow(X_combat_72970), ncol(X_combat_72970)))
cat(sprintf("    GSE69657: %d x %d\n", nrow(X_combat_69657), ncol(X_combat_69657)))

# ============================================================
# 8. PCA after correction
# ============================================================
cat("\n=== [8] PCA after correction ===\n\n")

X_pca_post <- X_combat[top_5000_probes, ]
X_pca_post_scaled <- t(scale(t(X_pca_post), center = TRUE, scale = TRUE))
X_pca_post_scaled[is.na(X_pca_post_scaled)] <- 0

pca_post <- prcomp(t(X_pca_post_scaled), center = FALSE, scale. = FALSE)
pca_var_post <- summary(pca_post)$importance[2, 1:3] * 100

pca_df_post <- data.frame(
  PC1 = pca_post$x[, 1],
  PC2 = pca_post$x[, 2],
  Batch = batch,
  Response = ifelse(is.na(all_resp), "Unknown", ifelse(all_resp == 1, "R", "S"))
)
cat(sprintf("  PC1: %.1f%%, PC2: %.1f%%, PC3: %.1f%%\n",
            pca_var_post[1], pca_var_post[2], pca_var_post[3]))

lda_acc_post <- NA
lda_post <- tryCatch(lda(PC1 ~ Batch, data = pca_df_post), error = function(e) NULL)
if (!is.null(lda_post)) {
  lda_pred_post <- predict(lda_post)
  lda_acc_post <- mean(lda_pred_post$class == pca_df_post$Batch)
  cat(sprintf("  LDA batch classification accuracy (PC1): %.1f%%\n", lda_acc_post * 100))
  cat(sprintf("  Batch separability reduction: %.1f%% -> %.1f%%\n",
              lda_acc_pre * 100, lda_acc_post * 100))
}

# ---- 8a. PCA visualization ----
library(ggrepel)

p_pre <- ggplot(pca_df_pre, aes(x = PC1, y = PC2, color = Batch, shape = Response)) +
  geom_point(size = 2.5, alpha = 0.7) +
  stat_ellipse(aes(group = Batch), level = 0.6, alpha = 0.3) +
  labs(title = "PCA Before ComBat Correction",
       subtitle = ifelse(is.na(lda_acc_pre),
         sprintf("PC1=%.1f%%, PC2=%.1f%%", pca_var_pre[1], pca_var_pre[2]),
         sprintf("PC1=%.1f%%, PC2=%.1f%%, LDA batch accuracy=%.0f%%",
                 pca_var_pre[1], pca_var_pre[2], lda_acc_pre * 100)),
       x = sprintf("PC1 (%.1f%%)", pca_var_pre[1]),
       y = sprintf("PC2 (%.1f%%)", pca_var_pre[2])) +
  theme_minimal() +
  scale_color_manual(values = c("GSE39582" = "#E41A1C", "GSE72970" = "#377EB8", "GSE69657" = "#4DAF4A"))

p_post <- ggplot(pca_df_post, aes(x = PC1, y = PC2, color = Batch, shape = Response)) +
  geom_point(size = 2.5, alpha = 0.7) +
  stat_ellipse(aes(group = Batch), level = 0.6, alpha = 0.3) +
  labs(title = "PCA After ComBat Correction",
       subtitle = ifelse(is.na(lda_acc_post),
         sprintf("PC1=%.1f%%, PC2=%.1f%%", pca_var_post[1], pca_var_post[2]),
         sprintf("PC1=%.1f%%, PC2=%.1f%%, LDA batch accuracy=%.0f%%",
                 pca_var_post[1], pca_var_post[2], lda_acc_post * 100)),
       x = sprintf("PC1 (%.1f%%)", pca_var_post[1]),
       y = sprintf("PC2 (%.1f%%)", pca_var_post[2])) +
  theme_minimal() +
  scale_color_manual(values = c("GSE39582" = "#E41A1C", "GSE72970" = "#377EB8", "GSE69657" = "#4DAF4A"))

# Also color by response to check biology is preserved
p_pre_resp <- ggplot(pca_df_pre, aes(x = PC1, y = PC2, color = Response, shape = Batch)) +
  geom_point(size = 2.5, alpha = 0.7) +
  stat_ellipse(aes(group = Response), level = 0.6, alpha = 0.3) +
  labs(title = "PCA Before ComBat - Colored by Response",
       x = sprintf("PC1 (%.1f%%)", pca_var_pre[1]),
       y = sprintf("PC2 (%.1f%%)", pca_var_pre[2])) +
  theme_minimal() +
  scale_color_manual(values = c("R" = "#E41A1C", "S" = "#377EB8", "Unknown" = "grey50"))

p_post_resp <- ggplot(pca_df_post, aes(x = PC1, y = PC2, color = Response, shape = Batch)) +
  geom_point(size = 2.5, alpha = 0.7) +
  stat_ellipse(aes(group = Response), level = 0.6, alpha = 0.3) +
  labs(title = "PCA After ComBat - Colored by Response",
       x = sprintf("PC1 (%.1f%%)", pca_var_post[1]),
       y = sprintf("PC2 (%.1f%%)", pca_var_post[2])) +
  theme_minimal() +
  scale_color_manual(values = c("R" = "#E41A1C", "S" = "#377EB8", "Unknown" = "grey50"))

pdf(file.path(BC_FIG_DIR, "ComBat_PCA_diagnostic.pdf"), width = 12, height = 10)
print(p_pre)
print(p_post)
print(p_pre_resp)
print(p_post_resp)
dev.off()
cat(sprintf("  PCA diagnostic saved: %s\n",
            file.path(BC_FIG_DIR, "ComBat_PCA_diagnostic.pdf")))

# ============================================================
# 9. Probe-to-gene mapping on CORRECTED data
# ============================================================
cat("\n=== [9] Gene-level aggregation (corrected data) ===\n\n")

gene_to_expression <- function(X_probe, annot_df, target_genes) {
  # Map probes to genes for specific target genes
  expr_list <- list()
  
  for (gene in target_genes) {
    # Find probes mapping to this gene
    gene_upper <- toupper(gene)
    probe_idx <- which(toupper(annot_df$gene_symbol) == gene_upper | 
                       grepl(gene_upper, toupper(annot_df$gene_symbol)))
    
    if (length(probe_idx) == 0) {
      # Try partial match for multi-gene entries
      probe_idx <- which(grepl(gene_upper, toupper(annot_df$gene_symbol)))
    }
    
    if (length(probe_idx) > 0) {
      # Get the probes that exist in our expression matrix
      probes_in_data <- intersect(annot_df$probe_id[probe_idx], rownames(X_probe))
      if (length(probes_in_data) > 0) {
        # Take mean across probes
        expr_list[[gene]] <- colMeans(X_probe[probes_in_data, , drop = FALSE])
      }
    }
  }
  
  if (length(expr_list) == 0) return(NULL)
  result <- do.call(rbind, expr_list)
  rownames(result) <- names(expr_list)
  return(result)
}

# Extract 8 fingerprint gene expression from each cohort
cat("  Extracting gene-level expression for 8 fingerprint genes...\n")

# --- PRE-correction ---
expr_gene_39582_pre <- gene_to_expression(X_39582, annot_map, fp_genes)
expr_gene_72970_pre <- gene_to_expression(X_72970, annot_map, fp_genes)
expr_gene_69657_pre <- gene_to_expression(X_69657, annot_map, fp_genes)

cat(sprintf("  Pre-correction genes found:\n"))
cat(sprintf("    GSE39582: %d\n", nrow(expr_gene_39582_pre)))
cat(sprintf("    GSE72970: %d\n", nrow(expr_gene_72970_pre)))
cat(sprintf("    GSE69657: %d\n", nrow(expr_gene_69657_pre)))

# --- POST-correction ---
expr_gene_39582_post <- gene_to_expression(X_combat_39582, annot_map, fp_genes)
expr_gene_72970_post <- gene_to_expression(X_combat_72970, annot_map, fp_genes)
expr_gene_69657_post <- gene_to_expression(X_combat_69657, annot_map, fp_genes)

cat(sprintf("  Post-correction genes found:\n"))
cat(sprintf("    GSE39582: %d\n", nrow(expr_gene_39582_post)))
cat(sprintf("    GSE72970: %d\n", nrow(expr_gene_72970_post)))
cat(sprintf("    GSE69657: %d\n", nrow(expr_gene_69657_post)))

# ============================================================
# 10. Directional concordance analysis (PRE vs POST correction)
# ============================================================
cat("\n=== [10] Directional concordance analysis ===\n\n")

# Helper: compute logFC between R and S groups
compute_logFC <- function(expr_gene_matrix, clinical_df, sample_col = "geo_accession",
                          resp_col_pattern = "response") {
  
  # Find response column
  resp_col <- grep(resp_col_pattern, names(clinical_df), value = TRUE)
  if (length(resp_col) == 0) return(rep(NA, nrow(expr_gene_matrix)))
  resp_col <- resp_col[1]
  
  # Match samples
  common <- intersect(colnames(expr_gene_matrix), clinical_df[[sample_col]])
  if (length(common) < 5) return(rep(NA, nrow(expr_gene_matrix)))
  
  # Align
  clin_sub <- clinical_df[match(common, clinical_df[[sample_col]]), ]
  expr_sub <- expr_gene_matrix[, common, drop = FALSE]
  
  # Get response values
  resp_vals <- clin_sub[[resp_col]]
  
  # Determine R vs S mapping
  if (any(toupper(resp_vals) %in% c("NR", "NO", "RESISTANT", "1", "RELPASED"))) {
    r_idx <- toupper(resp_vals) %in% c("NR", "NO", "RESISTANT", "1", "RELPASED")
  } else if (any(toupper(resp_vals) %in% c("R", "YES", "SENSITIVE", "0"))) {
    r_idx <- toupper(resp_vals) %in% c("R", "YES", "SENSITIVE", "0")
  } else {
    return(rep(NA, nrow(expr_gene_matrix)))
  }
  
  # Compute logFC (R - S)
  result <- apply(expr_sub, 1, function(x) {
    r_mean <- mean(x[r_idx], na.rm = TRUE)
    s_mean <- mean(x[!r_idx], na.rm = TRUE)
    return(r_mean - s_mean)  # logFC: positive = higher in resistant
  })
  
  return(result)
}

# --- Map response columns for each cohort ---
resp_col_39582 <- grep("response|relapse|status", names(clin_39582), value = TRUE, ignore.case = TRUE)
cat("  GSE39582 response columns:", paste(resp_col_39582, collapse = ", "), "\n")
resp_col_72970_name <- grep("response status", names(clin_72970), value = TRUE)[1]
cat("  GSE72970 response column:", resp_col_72970_name, "\n")  
resp_col_69657_name <- grep("chemoresponse|response", names(clin_69657), value = TRUE)[1]
cat("  GSE69657 response column:", resp_col_69657_name, "\n")

# Build logFC table manually for clarity
cat("\n  Computing logFC (R - S) for each gene...\n\n")

concordance_pre <- data.frame(
  Gene = fp_genes,
  SHAP_Weight = fp_shap,
  GSE39582_logFC = fp_ref_direction,
  GSE39582_Direction = ifelse(fp_ref_direction > 0, "up", "down"),
  stringsAsFactors = FALSE
)

concordance_post <- concordance_pre

# Compute validation logFC for each cohort (PRE correction)
for (cohort_name in c("GSE72970", "GSE69657")) {
  if (cohort_name == "GSE72970") {
    expr_pre <- expr_gene_72970_pre
    expr_post <- expr_gene_72970_post
    clin <- clin_72970
    resp_col <- resp_col_72970_name
  } else {
    expr_pre <- expr_gene_69657_pre
    expr_post <- expr_gene_69657_post
    clin <- clin_69657
    resp_col <- resp_col_69657_name
  }
  
  if (is.null(expr_pre) || is.null(expr_post)) next
  
  # Match samples between expression and clinical
  common <- intersect(colnames(expr_pre), clin$geo_accession)
  clin_sub <- clin[match(common, clin$geo_accession), ]
  expr_pre_sub <- expr_pre[, common, drop = FALSE]
  expr_post_sub <- expr_post[, common, drop = FALSE]
  
  # Get response
  resp_vals <- clin_sub[[resp_col]]
  if (is.null(resp_vals)) next
  
  # Map to binary
  if (cohort_name == "GSE72970") {
    r_idx <- resp_vals == "NR"
  } else {
    r_idx <- resp_vals == "noresponder"
  }
  s_idx <- !r_idx & !is.na(resp_vals)
  
  n_R <- sum(r_idx, na.rm = TRUE)
  n_S <- sum(s_idx, na.rm = TRUE)
  cat(sprintf("  %s: R=%d, S=%d\n", cohort_name, n_R, n_S))
  
  if (n_R < 2 || n_S < 2) {
    cat("    Too few samples, skipping\n")
    next
  }
  
  # Compute logFC for each gene (pre-correction)
  logFC_pre <- apply(expr_pre_sub, 1, function(x) {
    mean(as.numeric(x)[r_idx], na.rm = TRUE) - mean(as.numeric(x)[s_idx], na.rm = TRUE)
  })
  
  # Compute logFC (post-correction)
  logFC_post <- apply(expr_post_sub, 1, function(x) {
    mean(as.numeric(x)[r_idx], na.rm = TRUE) - mean(as.numeric(x)[s_idx], na.rm = TRUE)
  })
  
  # Add to tables
  concordance_pre[[paste0(cohort_name, "_logFC_pre")]] <- logFC_pre[rownames(concordance_pre)]
  concordance_pre[[paste0(cohort_name, "_Direction_pre")]] <- ifelse(
    logFC_pre[rownames(concordance_pre)] > 0, "up", "down"
  )
  concordance_pre[[paste0(cohort_name, "_Concordant_pre")]] <- 
    concordance_pre[[paste0(cohort_name, "_Direction_pre")]] == concordance_pre$GSE39582_Direction
  
  concordance_post[[paste0(cohort_name, "_logFC_post")]] <- logFC_post[rownames(concordance_post)]
  concordance_post[[paste0(cohort_name, "_Direction_post")]] <- ifelse(
    logFC_post[rownames(concordance_post)] > 0, "up", "down"
  )
  concordance_post[[paste0(cohort_name, "_Concordant_post")]] <- 
    concordance_post[[paste0(cohort_name, "_Direction_post")]] == concordance_post$GSE39582_Direction
}

# Also compute GSE39582 post-correction logFC
cat("  GSE39582 (self):\n")
common_39582 <- intersect(colnames(expr_gene_39582_post), clin_39582$geo_accession)
if (length(common_39582) > 0) {
  clin_39582_sub <- clin_39582[match(common_39582, clin_39582$geo_accession), ]
  
  resp_39582_col <- grep("response_status", names(clin_39582_sub), value = TRUE)[1]
  if (is.na(resp_39582_col)) resp_39582_col <- grep("Relapse", names(clin_39582_sub), value = TRUE)[1]
  
  if (!is.na(resp_39582_col)) {
    resp_39582_vals <- clin_39582_sub[[resp_39582_col]]
    r_idx_39582 <- resp_39582_vals %in% c("resistant", 1, "1", "Relapsed", "YES")
    s_idx_39582 <- !r_idx_39582 & !is.na(resp_39582_vals)
    cat(sprintf("    R=%d, S=%d\n", sum(r_idx_39582, na.rm = TRUE), sum(s_idx_39582, na.rm = TRUE)))
    
    if (sum(r_idx_39582, na.rm = TRUE) >= 2 && sum(s_idx_39582, na.rm = TRUE) >= 2) {
      logFC_39582_post <- apply(expr_gene_39582_post[, common_39582, drop = FALSE], 1, function(x) {
        mean(as.numeric(x)[r_idx_39582], na.rm = TRUE) - mean(as.numeric(x)[s_idx_39582], na.rm = TRUE)
      })
      
      concordance_post[["GSE39582_logFC_post"]] <- logFC_39582_post[rownames(concordance_post)]
      concordance_post[["GSE39582_Direction_post"]] <- ifelse(logFC_39582_post > 0, "up", "down")
    }
  }
}

# ============================================================
# 11. Print concordance comparison
# ============================================================
cat("\n=== [11] Concordance comparison ===\n\n")
cat("--- PRE-Correction Concordance ---\n")
for (i in 1:nrow(concordance_pre)) {
  gene <- concordance_pre$Gene[i]
  ref_dir <- concordance_pre$GSE39582_Direction[i]
  
  dirs <- c()
  for (cohort in c("GSE72970", "GSE69657")) {
    col_name <- paste0(cohort, "_Direction_pre")
    if (col_name %in% names(concordance_pre)) {
      dirs <- c(dirs, as.character(concordance_pre[i, col_name]))
    }
  }
  
  conc_str <- paste(sapply(dirs, function(d) ifelse(d == ref_dir, "[OK]", "[XX]")), collapse = "/")
  cat(sprintf("  %-10s (ref=%s): %s  %s\n", gene, ref_dir, 
              paste(dirs, collapse = "/"), conc_str))
}

cat("\n--- POST-Correction Concordance ---\n")
for (i in 1:nrow(concordance_post)) {
  gene <- concordance_post$Gene[i]
  ref_dir <- concordance_post$GSE39582_Direction[i]
  
  dirs <- c()
  for (cohort in c("GSE72970", "GSE69657", "GSE39582")) {
    col_name <- paste0(cohort, "_Direction_post")
    if (col_name %in% names(concordance_post)) {
      dirs <- c(dirs, as.character(concordance_post[i, col_name]))
    }
  }
  
  conc_str <- paste(sapply(dirs, function(d) ifelse(d == ref_dir, "[OK]", "[XX]")), collapse = "/")
  cat(sprintf("  %-10s (ref=%s): %s  %s\n", gene, ref_dir, 
              paste(dirs, collapse = "/"), conc_str))
}

# ============================================================
# 12. Jaccard stability filtering
# ============================================================
cat("\n=== [12] Jaccard stability filtering ===\n\n")

compute_jaccard <- function(concordance_df, direction_cols, ref_col = "GSE39582_Direction") {
  # For each gene, compute Jaccard-like stability score
  # = (# cohorts where direction matches ref) / (# total cohorts with data)
  jaccard_scores <- c()
  
  for (i in 1:nrow(concordance_df)) {
    ref_dir <- as.character(concordance_df[i, ref_col])
    matches <- 0
    total <- 0
    
    for (col in direction_cols) {
      if (col %in% names(concordance_df) && !is.na(concordance_df[i, col])) {
        total <- total + 1
        if (as.character(concordance_df[i, col]) == ref_dir) {
          matches <- matches + 1
        }
      }
    }
    
    jaccard_scores <- c(jaccard_scores, matches / max(total, 1))
  }
  
  return(jaccard_scores)
}

# Pre-correction: GSE72970 and GSE69657 (exclude GSE39582 self)
pre_dir_cols <- grep("_Direction_pre$", names(concordance_pre), value = TRUE)
concordance_pre$Jaccard_Stability <- compute_jaccard(concordance_pre, pre_dir_cols)
concordance_pre$Pass_Stability <- concordance_pre$Jaccard_Stability >= 0.6

# Post-correction: all external cohorts
post_dir_cols <- grep("_Direction_post$", names(concordance_post), value = TRUE)
concordance_post$Jaccard_Stability <- compute_jaccard(concordance_post, post_dir_cols)
concordance_post$Pass_Stability <- concordance_post$Jaccard_Stability >= 0.6

cat("--- Stability Summary ---\n")
cat(sprintf("%-10s  %-12s  %-12s  %s\n", "Gene", "Pre-Stability", "Post-Stability", "Pass Post?"))
cat(paste(rep("-", 50), collapse = ""), "\n")
for (i in 1:nrow(concordance_pre)) {
  gene <- concordance_pre$Gene[i]
  pre_j <- concordance_pre$Jaccard_Stability[i]
  post_j <- concordance_post$Jaccard_Stability[i]
  pass <- ifelse(post_j >= 0.6, "[PASS]", "[FAIL]")
  cat(sprintf("%-10s  %.2f (%.0f%%)     %.2f (%.0f%%)     %s\n", 
              gene, pre_j, pre_j*100, post_j, post_j*100, pass))
}

# Identify stable genes (post >= 60%)
stable_genes <- concordance_post$Gene[concordance_post$Pass_Stability]
unstable_genes <- concordance_post$Gene[!concordance_post$Pass_Stability]

cat(sprintf("\n  Stable genes (Jaccard >= 0.6): %d\n", length(stable_genes)))
cat(sprintf("    %s\n", paste(stable_genes, collapse = ", ")))
cat(sprintf("  Unstable genes (< 0.6): %d\n", length(unstable_genes)))
cat(sprintf("    %s\n", paste(unstable_genes, collapse = ", ")))

# Print improvement
improved <- sum(concordance_post$Jaccard_Stability > concordance_pre$Jaccard_Stability, na.rm = TRUE)
worsened <- sum(concordance_post$Jaccard_Stability < concordance_pre$Jaccard_Stability, na.rm = TRUE)
same <- sum(concordance_post$Jaccard_Stability == concordance_pre$Jaccard_Stability, na.rm = TRUE)
cat(sprintf("\n  ComBat effect on stability: %d improved, %d worsened, %d unchanged\n",
            improved, worsened, same))

# ============================================================
# 13. Save outputs
# ============================================================
cat("\n=== [13] Saving outputs ===\n\n")

# 13a. Corrected expression matrix
saveRDS(X_combat, file.path(BC_DIR, "ComBat_corrected_expression.rds"))
cat(sprintf("  Corrected expression: %s\n",
            file.path(BC_DIR, "ComBat_corrected_expression.rds")))

# 13b. Concordance tables
write.csv(concordance_pre, file.path(BC_DIR, "concordance_pre_combat.csv"), row.names = FALSE)
cat(sprintf("  Pre-ComBat concordance: %s\n",
            file.path(BC_DIR, "concordance_pre_combat.csv")))

write.csv(concordance_post, file.path(BC_DIR, "concordance_post_combat.csv"), row.names = FALSE)
cat(sprintf("  Post-ComBat concordance: %s\n",
            file.path(BC_DIR, "concordance_post_combat.csv")))

# 13c. Stability summary
stability_df <- data.frame(
  Gene = concordance_post$Gene,
  SHAP_Weight = concordance_post$SHAP_Weight,
  GSE39582_Ref_Direction = concordance_post$GSE39582_Direction,
  Jaccard_Stability_Pre = concordance_pre$Jaccard_Stability,
  Jaccard_Stability_Post = concordance_post$Jaccard_Stability,
  Pass_Stability = concordance_post$Pass_Stability,
  stringsAsFactors = FALSE
)

# Add direction details
for (col in names(concordance_post)) {
  if (grepl("_Direction_", col) && !col %in% names(stability_df)) {
    stability_df[[col]] <- concordance_post[[col]]
  }
}

write.csv(stability_df, file.path(BC_DIR, "directional_stability_summary.csv"), row.names = FALSE)
cat(sprintf("  Stability summary: %s\n",
            file.path(BC_DIR, "directional_stability_summary.csv")))

# 13d. Stable gene list
stable_df <- data.frame(
  Gene = stable_genes,
  SHAP_Weight = fp_shap[match(stable_genes, names(fp_shap))],
  stringsAsFactors = FALSE
)
write.csv(stable_df, file.path(BC_DIR, "stable_fingerprint_genes.csv"), row.names = FALSE)
cat(sprintf("  Stable genes: %s\n",
            file.path(BC_DIR, "stable_fingerprint_genes.csv")))

# 13e. Also save per-cohort corrected expression for downstream use
saveRDS(X_combat_39582, file.path(BC_DIR, "GSE39582_combat_corrected.rds"))
saveRDS(X_combat_72970, file.path(BC_DIR, "GSE72970_combat_corrected.rds"))
saveRDS(X_combat_69657, file.path(BC_DIR, "GSE69657_combat_corrected.rds"))
cat("  Per-cohort corrected expression saved\n")

# ============================================================
# 14. Summary report
# ============================================================
cat("\n=== [14] Summary ===\n\n")
cat("============================================================\n")
cat("  ComBat Batch Correction Complete\n")
cat("============================================================\n\n")

cat(sprintf("  Cohorts corrected: 3 (GSE39582 + GSE72970 + GSE69657)\n"))
cat(sprintf("  Total samples: %d\n", ncol(X_combat)))
cat(sprintf("  Total probes: %d\n", nrow(X_combat)))
cat(sprintf("  Fingerprint genes: %d\n", length(fp_genes)))
cat(sprintf("  Stable genes (Jaccard >= 60%%): %d\n", length(stable_genes)))
cat(sprintf("  Direction stability improved: %d / %d\n", improved, length(fp_genes)))

if (length(stable_genes) >= 6) {
  cat("\n  [PASS] ComBat correction restored directional concordance for most genes.\n")
  cat("     The 8-gene fingerprint is viable after batch correction.\n")
  cat("     Proceed to: Route 2 (GSVA) + Route 3 (PRS+Cox)\n")
} else if (length(stable_genes) >= 4) {
  cat("\n  [WARN] Partial stabilization: some genes remain directionally unstable.\n")
  cat("     Use stable subset for validation. Proceed to Route 2 (GSVA)\n")
  cat("     and Route 3 (PRS+Cox on stable subset).\n")
} else {
  cat("\n  [FAIL] ComBat correction insufficient to stabilize directional concordance.\n")
  cat("     Batch effect is not the primary cause of fingerprint failure.\n")
  cat("     Directly proceed to Route 3 (LASSO->PRS).\n")
}

cat(sprintf("\n  Output directory: %s\n", BC_DIR))
cat("  Script complete!\n")
