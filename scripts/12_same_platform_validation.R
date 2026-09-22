# ============================================================
# 12_SAME_PLATFORM_VALIDATION.R
# XELOX Resistance - GPL570 Same-Platform Validation (GSE72970 + GSE69657)
# ============================================================
# Purpose:
#   Validate the 8-gene fingerprint on GPL570 independent cohorts
#   Strategy A: SHAP-weighted risk score (z-score normalization)
#   Strategy B: Simple mean z-score risk score
#   Strategy C: Directional concordance (logFC comparison with GSE39582)
#
# Primary target: GSE72970 FOLFOX subset (36 samples, RECIST response)
# Secondary target: GSE69657 FOLFOX4 (30 samples, pathological response)
#
# Output:
#   tables/GSE72970_FOLFOX_predictions.csv
#   tables/GSE69657_FOLFOX4_predictions.csv
#   tables/same_platform_validation_summary.csv
#   figures/GSE72970_validation_roc.pdf
#   figures/GSE69657_validation_roc.pdf
#   figures/GSE72970_direction_concordance.pdf
# ============================================================

Sys.setenv(TMPDIR = "C:/temp", TMP = "C:/temp", TEMP = "C:/temp")
.libPaths(c("C:/Rlibs", .libPaths()))

# Accept PROJECT_ROOT from environment variable (avoids Chinese path issues on cmd line)
PROJECT_ROOT <- Sys.getenv("XELOX_ROOT",
  "/path/to/xelox_project")

# Set locale AFTER defining paths (env var approach avoids Chinese char in cmd arg)
if (.Platform$OS.type == "windows") {
  tryCatch(Sys.setlocale("LC_ALL", "Chinese"), error = function(e) {
    tryCatch(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), error = function(e2) {})
  })
}
DATA_GEO_DIR     <- file.path(PROJECT_ROOT, "data", "geo")
RESULTS_TAB_DIR  <- file.path(PROJECT_ROOT, "results", "tables")
RESULTS_FIG_DIR  <- file.path(PROJECT_ROOT, "results", "figures")
ML2_TAB_DIR      <- file.path(RESULTS_TAB_DIR, "ml_phase2")

set.seed(42)

cat("============================================================\n")
cat("GPL570 Same-Platform Validation of 8-Gene Fingerprint\n")
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
# 1. Load Phase 2 outputs (SHAP weights + fingerprint genes)
# ============================================================
cat("=== [1] Loading Phase 2 outputs ===\n\n")

# 1a. SHAP summary
shap_file <- file.path(ML2_TAB_DIR, "SHAP_summary_phase2.csv")
shap_df <- read.csv(shap_file, stringsAsFactors = FALSE)
rownames(shap_df) <- shap_df$Gene
cat(sprintf("  SHAP summary: %d genes\n", nrow(shap_df)))

# 1b. Fingerprint genes (top 10)
fingerprint_file <- file.path(ML2_TAB_DIR, "XELOX_resistance_fingerprint.csv")
fingerprint_df <- read.csv(fingerprint_file, stringsAsFactors = FALSE)
fp_genes_all <- fingerprint_df$Gene
cat(sprintf("  Fingerprint genes (full): %d\n", length(fp_genes_all)))

# 1c. Define the 8 usable fingerprint genes (exclude C5ORF42, GGT1/GGT2)
# C5ORF42: poorly characterized, excluded from paper narrative
# GGT1/GGT2: multi-gene entry "GGT1 /// GGT2 /// GGTLC1 /// GGTLC2"
fp_genes <- c("CSNK1G2", "KAZN", "KLK6", "MGA", "MID2", "SOX11", "TAS2R40", "ZNF451")
cat("  8 usable fingerprint genes:\n")
for (g in fp_genes) cat(sprintf("    - %s\n", g))

# ============================================================
# 2. Load GSE39582 direction info from DEG data
# ============================================================
cat("\n=== [2] Loading GSE39582 direction information ===\n\n")

deg_file <- file.path(RESULTS_TAB_DIR, "DEG_GSE39582_mapped.csv")
deg_df <- read.csv(deg_file, stringsAsFactors = FALSE)

# Build gene -> logFC map (take most significant probe per gene)
deg_clean <- deg_df[!is.na(deg_df$gene_symbol) & deg_df$gene_symbol != "", ]
gene_order <- order(deg_clean$P.Value)
deg_clean <- deg_clean[gene_order, ]
deg_dedup <- deg_clean[!duplicated(deg_clean$gene_symbol), ]
gene_direction <- setNames(deg_dedup$logFC, deg_dedup$gene_symbol)

# Helper: resolve multi-gene entry
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

# Get directions for our 8 genes
cat("  Fingerprint gene directions (logFC in GSE39582, resistant vs sensitive):\n")
fp_directions <- sapply(fp_genes, function(g) resolve_direction(g, gene_direction))
names(fp_directions) <- fp_genes
for (i in seq_along(fp_genes)) {
  cat(sprintf("    %s: logFC = %s\n", fp_genes[i],
              ifelse(is.na(fp_directions[i]), "NA", sprintf("%+.4f", fp_directions[i]))))
}

# SHAP weights for the 8 genes (matching multi-gene entries)
get_shap_weight <- function(gene, shap_df) {
  # Try exact match
  if (gene %in% rownames(shap_df)) return(shap_df[gene, "MeanAbsSHAP"])
  # Try fuzzy match (gene might be part of a multi-gene entry)
  for (g in rownames(shap_df)) {
    parts <- trimws(strsplit(g, " /// ")[[1]])
    if (toupper(gene) %in% toupper(parts)) return(shap_df[g, "MeanAbsSHAP"])
  }
  return(NA)
}

fp_shap <- sapply(fp_genes, function(g) get_shap_weight(g, shap_df))
cat("\n  SHAP weights for 8 genes:\n")
for (i in seq_along(fp_genes)) {
  cat(sprintf("    %s: SHAP = %s\n", fp_genes[i],
              ifelse(is.na(fp_shap[i]), "NA", sprintf("%.4f", fp_shap[i]))))
}

# ============================================================
# 3. Load GPL570 probe-to-gene annotation
# ============================================================
cat("\n=== [3] Loading GPL570 probe-gene mapping ===\n\n")

annot_file <- file.path(DATA_GEO_DIR, "GPL570_probe_gene_map.csv")
annot_map <- read.csv(annot_file, stringsAsFactors = FALSE)
cat(sprintf("  GPL570 probe-gene mapping: %d probes mapped\n", nrow(annot_map)))
cat(sprintf("  Unique genes: %d\n", length(unique(annot_map$gene_symbol))))
cat(sprintf("  Sample: %s -> %s\n", annot_map$probe_id[1], annot_map$gene_symbol[1]))

# ============================================================
# 4. Helper function: validate on one dataset
# ============================================================

validate_cohort <- function(expr_matrix, clinical_df, cohort_name,
                            regimen_filter = NULL, response_col = "response status:ch1",
                            response_map = c("NR" = 1, "R" = 0),
                            pfs_col = NULL, pfs_censor_col = NULL) {

  cat(sprintf("\n%s\n", paste(rep("=", 60), collapse = "")))
  cat(sprintf("  === %s ===\n", cohort_name))
  cat(sprintf("%s\n\n", paste(rep("=", 60), collapse = "")))

  # ---- 4a. Sample matching ----
  gsm_expr <- colnames(expr_matrix)
  gsm_clin <- clinical_df$geo_accession
  common_gsm <- intersect(gsm_expr, gsm_clin)
  cat(sprintf("  Samples in expression matrix: %d\n", length(gsm_expr)))
  cat(sprintf("  Samples in clinical data: %d\n", length(gsm_clin)))
  cat(sprintf("  Matched samples: %d\n", length(common_gsm)))

  # Subset to common samples
  clin_matched <- clinical_df[match(common_gsm, gsm_clin), , drop = FALSE]
  expr_matched <- expr_matrix[, common_gsm, drop = FALSE]

  # ---- 4b. Regimen filter (if specified) ----
  if (!is.null(regimen_filter)) {
    col_regimen <- grep("regimen", colnames(clin_matched), value = TRUE)[1]
    cat(sprintf("  Regimen column: %s\n", col_regimen))

    # Check available regimens
    cat("  Available regimens:\n")
    print(table(clin_matched[[col_regimen]]))

    # Apply filter
    filter_idx <- clin_matched[[col_regimen]] %in% regimen_filter
    clin_filtered <- clin_matched[filter_idx, , drop = FALSE]
    expr_filtered <- expr_matched[, filter_idx, drop = FALSE]
    cat(sprintf("\n  After regimen filter [%s]: %d samples\n",
                paste(regimen_filter, collapse = "+"), sum(filter_idx)))
  } else {
    clin_filtered <- clin_matched
    expr_filtered <- expr_matched
    cat("  No regimen filter applied\n")
  }

  # ---- 4c. Response labels ----
  col_resp <- grep(response_col, colnames(clin_filtered), value = TRUE)[1]
  if (is.na(col_resp)) {
    # Try alternative response column names
    alt_cols <- grep("response|responder", colnames(clin_filtered), value = TRUE)
    if (length(alt_cols) > 0) col_resp <- alt_cols[1]
  }
  cat(sprintf("  Response column: %s\n", col_resp))

  # Build R/S labels
  response_vals <- clin_filtered[[col_resp]]
  cat("  Response value distribution:\n")
  print(table(response_vals))

  # Map to binary: 1 = resistant, 0 = sensitive
  y <- ifelse(response_vals %in% names(response_map[response_map == 1]),
              1,
              ifelse(response_vals %in% names(response_map[response_map == 0]), 0, NA))
  na_idx <- is.na(y)
  cat(sprintf("  Binary labels: R=%d, S=%d, NA=%d\n",
              sum(y == 1, na.rm = TRUE), sum(y == 0, na.rm = TRUE), sum(na_idx)))

  if (sum(!na_idx) < 10) {
    cat("  ERROR: Too few labeled samples\n")
    return(NULL)
  }

  clin_final <- clin_filtered[!na_idx, , drop = FALSE]
  expr_final <- expr_filtered[, !na_idx, drop = FALSE]
  y_final <- y[!na_idx]

  # PFS data (if available)
  pfs <- NULL
  if (!is.null(pfs_col)) {
    pfs_time <- clin_final[[pfs_col]]
    pfs_event <- clin_final[[pfs_censor_col]]
    # Handle format: might have "pfs:ch1" prefix in column names
    if (is.null(pfs_time) && !is.null(pfs_col)) {
      pfs_col_full <- grep(pfs_col, colnames(clin_final), value = TRUE)[1]
      if (!is.na(pfs_col_full)) pfs_time <- clin_final[[pfs_col_full]]
    }
    if (!is.null(pfs_time)) {
      pfs <- data.frame(time = as.numeric(pfs_time),
                        event = 1 - as.numeric(clin_final[[pfs_censor_col]]),
                        stringsAsFactors = FALSE)
      cat(sprintf("  PFS data: %d samples, median time = %.1f mo\n",
                  nrow(pfs), median(pfs$time, na.rm = TRUE)))
    }
  }

  # ---- 4d. Probe-to-gene mapping ----
  cat("\n  --- Probe-to-gene mapping ---\n")
  common_probes <- intersect(rownames(expr_final), annot_map$probe_id)
  cat(sprintf("  Probes in expression matrix: %d\n", nrow(expr_final)))
  cat(sprintf("  Probes matched to annotation: %d\n", length(common_probes)))

  expr_probes <- expr_final[common_probes, , drop = FALSE]
  probe_gene_map <- annot_map[match(common_probes, annot_map$probe_id), ]

  # Aggregate: for each gene, average all probes (robust approach)
  gene_expr_list <- split(seq_len(nrow(expr_probes)), probe_gene_map$gene_symbol)
  gene_names <- names(gene_expr_list)
  expr_gene_list <- lapply(gene_names, function(g) {
    idx <- gene_expr_list[[g]]
    if (length(idx) == 1) return(as.numeric(expr_probes[idx, , drop = FALSE]))
    return(colMeans(expr_probes[idx, , drop = FALSE]))
  })
  expr_gene <- do.call(rbind, expr_gene_list)
  rownames(expr_gene) <- gene_names
  colnames(expr_gene) <- colnames(expr_probes)
  cat(sprintf("  Unique genes after aggregation: %d\n", nrow(expr_gene)))

  # ---- 4e. Extract fingerprint genes ----
  cat("\n  --- Fingerprint gene extraction ---\n")
  genes_found <- intersect(fp_genes, rownames(expr_gene))
  cat(sprintf("  Fingerprint genes found: %d / %d\n", length(genes_found), length(fp_genes)))
  for (g in fp_genes) {
    cat(sprintf("    %s: %s\n", g, ifelse(g %in% genes_found, "FOUND", "MISSING")))
  }

  if (length(genes_found) < 3) {
    cat("  ERROR: Too few fingerprint genes found\n")
    return(NULL)
  }

  X_fp <- expr_gene[genes_found, , drop = FALSE]

  # ---- 4f. Strategy A: SHAP-weighted risk score ----
  cat("\n  --- Strategy A: SHAP-weighted risk score ---\n")

  weights <- sapply(genes_found, function(g) {
    w <- fp_shap[which(fp_genes == g)]
    if (is.na(w)) w <- mean(fp_shap, na.rm = TRUE)
    return(w)
  })
  dir_signs <- sapply(genes_found, function(g) {
    d <- fp_directions[which(fp_genes == g)]
    if (is.na(d)) return(0)
    return(sign(d))
  })
  names(weights) <- genes_found
  names(dir_signs) <- genes_found

  cat("  Gene weights and direction signs:\n")
  for (g in genes_found) {
    cat(sprintf("    %s: weight=%.4f, sign=%+d\n", g, weights[g], dir_signs[g]))
  }

  # Z-score normalize each gene
  z_scores <- t(apply(X_fp, 1, function(x) {
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(x)))
    return((x - mean(x, na.rm = TRUE)) / s)
  }))

  # Compute risk score
  risk_A <- colSums(z_scores * weights * dir_signs, na.rm = TRUE)
  risk_A[is.nan(risk_A) | is.infinite(risk_A)] <- 0
  cat(sprintf("  Risk Score A range: [%.4f, %.4f]\n", min(risk_A), max(risk_A)))
  cat(sprintf("  Mean (Resistant): %.4f, Mean (Sensitive): %.4f\n",
              mean(risk_A[y_final == 1]), mean(risk_A[y_final == 0])))

  # ROC for Strategy A
  auc_A <- NA; ci_A <- c(NA, NA, NA); wt_A_p <- NA
  if (length(unique(y_final)) == 2 && sd(risk_A) > 0) {
    roc_A <- roc(y_final, risk_A, quiet = TRUE)
    auc_A <- as.numeric(roc_A$auc)
    ci_A <- ci.auc(roc_A)
    wt_A <- wilcox.test(risk_A ~ y_final, alternative = "two.sided")
    wt_A_p <- wt_A$p.value
    cat(sprintf("  AUC (Strategy A): %.4f (95%% CI: %.4f - %.4f)\n", auc_A, ci_A[1], ci_A[3]))
    cat(sprintf("  Wilcoxon p-value: %.6f\n", wt_A_p))
  } else {
    cat("  SKIPPING ROC: only one class or no variation\n")
  }

  # ---- 4g. Strategy B: Simple mean z-score ----
  cat("\n  --- Strategy B: Simple mean z-score ---\n")
  risk_B <- colMeans(z_scores, na.rm = TRUE)
  cat(sprintf("  Risk Score B range: [%.4f, %.4f]\n", min(risk_B), max(risk_B)))
  cat(sprintf("  Mean (Resistant): %.4f, Mean (Sensitive): %.4f\n",
              mean(risk_B[y_final == 1]), mean(risk_B[y_final == 0])))

  # ROC for Strategy B
  auc_B <- NA; ci_B <- c(NA, NA, NA); wt_B_p <- NA
  if (length(unique(y_final)) == 2 && sd(risk_B) > 0) {
    roc_B <- roc(y_final, risk_B, quiet = TRUE)
    auc_B <- as.numeric(roc_B$auc)
    ci_B <- ci.auc(roc_B)
    wt_B <- wilcox.test(risk_B ~ y_final, alternative = "two.sided")
    wt_B_p <- wt_B$p.value
    cat(sprintf("  AUC (Strategy B): %.4f (95%% CI: %.4f - %.4f)\n", auc_B, ci_B[1], ci_B[3]))
    cat(sprintf("  Wilcoxon p-value: %.6f\n", wt_B_p))
  } else {
    cat("  SKIPPING ROC: only one class or no variation\n")
  }

  # ---- 4h. Strategy C: Directional concordance ----
  cat("\n  --- Strategy C: Directional concordance ---\n")

  # Compute logFC in validation cohort (R vs S)
  concordance <- data.frame(
    Gene = genes_found,
    GSE39582_logFC = NA_real_,
    GSE39582_Direction = "",
    Validation_Mean_R = NA_real_,
    Validation_Mean_S = NA_real_,
    Validation_logFC = NA_real_,
    Validation_Direction = "",
    Concordant = FALSE,
    SHAP_Weight = NA_real_,
    stringsAsFactors = FALSE
  )

  n_concordant <- 0
  n_total <- 0

  for (i in seq_along(genes_found)) {
    g <- genes_found[i]
    r_vals <- as.numeric(expr_gene[g, y_final == 1])
    s_vals <- as.numeric(expr_gene[g, y_final == 0])

    mean_r <- mean(r_vals, na.rm = TRUE)
    mean_s <- mean(s_vals, na.rm = TRUE)
    val_logFC <- mean_r - mean_s  # log scale (already log2 from RMA)
    val_dir <- ifelse(val_logFC > 0, "up", "down")

    gse_dir_val <- fp_directions[which(fp_genes == g)]
    gse_dir <- ifelse(is.na(gse_dir_val), "NA", ifelse(gse_dir_val > 0, "up", "down"))

    concordant <- !is.na(gse_dir_val) && gse_dir == val_dir
    if (concordant && !is.na(gse_dir_val)) {
      n_concordant <- n_concordant + 1
    }
    if (!is.na(gse_dir_val)) n_total <- n_total + 1

    concordance$GSE39582_logFC[i] <- ifelse(is.na(gse_dir_val), NA, round(gse_dir_val, 4))
    concordance$GSE39582_Direction[i] <- gse_dir
    concordance$Validation_Mean_R[i] <- round(mean_r, 4)
    concordance$Validation_Mean_S[i] <- round(mean_s, 4)
    concordance$Validation_logFC[i] <- round(val_logFC, 4)
    concordance$Validation_Direction[i] <- val_dir
    concordance$Concordant[i] <- concordant
    concordance$SHAP_Weight[i] <- round(weights[g], 4)
  }

  # Binomial test
  binom_p <- NA
  if (n_total > 0) {
    binom_p <- binom.test(n_concordant, n_total, p = 0.5, alternative = "greater")$p.value
  }
  cat(sprintf("  Concordant: %d / %d (%.1f%%)\n",
              n_concordant, n_total, 100 * n_concordant / max(n_total, 1)))
  cat(sprintf("  Binomial test (one-sided): p = %.4f\n", binom_p))
  print(concordance, row.names = FALSE)

  # ---- 4i. Compile results ----
  results <- list(
    cohort = cohort_name,
    n_total = sum(!na_idx),
    n_resistant = sum(y_final == 1, na.rm = TRUE),
    n_sensitive = sum(y_final == 0, na.rm = TRUE),
    n_genes = length(genes_found),
    auc_A = auc_A, ci_A_low = ci_A[1], ci_A_high = ci_A[3], wt_A_p = wt_A_p,
    auc_B = auc_B, ci_B_low = ci_B[1], ci_B_high = ci_B[3], wt_B_p = wt_B_p,
    n_concordant = n_concordant, n_total_concordance = n_total,
    concordance_pct = ifelse(n_total > 0, 100 * n_concordant / n_total, NA),
    binom_p = binom_p,
    risk_A = risk_A, risk_B = risk_B,
    y = y_final,
    concordance = concordance,
    pfs = pfs
  )

  return(results)
}

# ============================================================
# 5. Validate GSE72970 (FOLFOX subset)
# ============================================================
cat("\n\n")
cat("============================================================\n")
cat("  PHASE 5: GSE72970 Validation\n")
cat("============================================================\n\n")

# Load expression matrix
expr_72970 <- readRDS(file.path(DATA_GEO_DIR, "GSE72970_expression.rds"))
cat(sprintf("GSE72970 expression: %d probes x %d samples\n", nrow(expr_72970), ncol(expr_72970)))

# Load clinical data
clin_72970 <- read.csv(file.path(RESULTS_TAB_DIR, "GSE72970_clinical_data.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("GSE72970 clinical: %d samples\n", nrow(clin_72970)))

# Run validation for FOLFOX subset
res_72970 <- validate_cohort(
  expr_matrix = expr_72970,
  clinical_df = clin_72970,
  cohort_name = "GSE72970_FOLFOX",
  regimen_filter = c("FOLFOX", "FOLFOX+Bev"),
  response_col = "response status:ch1",
  response_map = c("NR" = 1, "R" = 0),
  pfs_col = "pfs:ch1",
  pfs_censor_col = "pfs censored:ch1"
)

# Also run on FOLFIRI subset for comparison (if FOLFOX results are interesting)
res_72970_folfiri <- NULL
if (!is.null(res_72970)) {
  res_72970_folfiri <- validate_cohort(
    expr_matrix = expr_72970,
    clinical_df = clin_72970,
    cohort_name = "GSE72970_FOLFIRI",
    regimen_filter = c("FOLFIRI", "FOLFIRI+Bev"),
    response_col = "response status:ch1",
    response_map = c("NR" = 1, "R" = 0),
    pfs_col = "pfs:ch1",
    pfs_censor_col = "pfs censored:ch1"
  )
}

# ============================================================
# 6. Validate GSE69657 (FOLFOX4, pathological response)
# ============================================================
cat("\n\n")
cat("============================================================\n")
cat("  PHASE 6: GSE69657 Validation\n")
cat("============================================================\n\n")

expr_69657 <- readRDS(file.path(DATA_GEO_DIR, "GSE69657_expression.rds"))
cat(sprintf("GSE69657 expression: %d probes x %d samples\n", nrow(expr_69657), ncol(expr_69657)))

clin_69657 <- read.csv(file.path(RESULTS_TAB_DIR, "GSE69657_clinical_data.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("GSE69657 clinical: %d samples\n", nrow(clin_69657)))

res_69657 <- validate_cohort(
  expr_matrix = expr_69657,
  clinical_df = clin_69657,
  cohort_name = "GSE69657_FOLFOX4",
  regimen_filter = NULL,  # All samples are FOLFOX4
  response_col = "chemoresponse:ch1",
  response_map = c("noresponder" = 1, "responder" = 0)  # noresponder=resistant, responder=sensitive
)

# ============================================================
# 7. Compile summary and save results
# ============================================================
cat("\n\n")
cat("============================================================\n")
cat("  PHASE 7: Summary & Output\n")
cat("============================================================\n\n")

# Collect all non-null results
all_results <- list()
if (!is.null(res_72970)) all_results[["GSE72970_FOLFOX"]] <- res_72970
if (!is.null(res_72970_folfiri)) all_results[["GSE72970_FOLFIRI"]] <- res_72970_folfiri
if (!is.null(res_69657)) all_results[["GSE69657_FOLFOX4"]] <- res_69657

# ---- 7a. Summary table ----
summary_df <- do.call(rbind, lapply(names(all_results), function(nm) {
  r <- all_results[[nm]]
  data.frame(
    Cohort = nm,
    N = r$n_total,
    R = r$n_resistant,
    S = r$n_sensitive,
    N_genes = r$n_genes,
    AUC_A = sprintf("%.4f", r$auc_A),
    AUC_A_CI = sprintf("%.4f-%.4f", r$ci_A_low, r$ci_A_high),
    AUC_A_p = sprintf("%.6f", r$wt_A_p),
    AUC_B = sprintf("%.4f", r$auc_B),
    AUC_B_CI = sprintf("%.4f-%.4f", r$ci_B_low, r$ci_B_high),
    AUC_B_p = sprintf("%.6f", r$wt_B_p),
    Concordance = sprintf("%d/%d (%.1f%%)", r$n_concordant, r$n_total_concordance, r$concordance_pct),
    Binom_p = sprintf("%.4f", r$binom_p),
    stringsAsFactors = FALSE
  )
}))

cat("\n--- Validation Summary ---\n")
print(summary_df, row.names = FALSE)

write.csv(summary_df,
          file.path(RESULTS_TAB_DIR, "same_platform_validation_summary.csv"),
          row.names = FALSE)
cat(sprintf("\n  Summary saved: %s\n",
            file.path(RESULTS_TAB_DIR, "same_platform_validation_summary.csv")))

# ---- 7b. Per-sample predictions ----
for (nm in names(all_results)) {
  r <- all_results[[nm]]
  pred_df <- data.frame(
    sample = names(r$risk_A),
    risk_score_A = round(r$risk_A, 6),
    risk_score_B = round(r$risk_B, 6),
    response = r$y,
    stringsAsFactors = FALSE
  )

  # Add PFS if available
  if (!is.null(r$pfs) && nrow(r$pfs) == length(r$y)) {
    pred_df$pfs_time <- r$pfs$time
    pred_df$pfs_event <- r$pfs$event
  }

  out_file <- file.path(RESULTS_TAB_DIR, paste0(nm, "_predictions.csv"))
  write.csv(pred_df, out_file, row.names = FALSE)
  cat(sprintf("  Predictions saved: %s\n", out_file))
}

# ---- 7c. Concordance tables ----
for (nm in names(all_results)) {
  r <- all_results[[nm]]
  conc_file <- file.path(RESULTS_TAB_DIR, paste0(nm, "_concordance.csv"))
  write.csv(r$concordance, conc_file, row.names = FALSE)
  cat(sprintf("  Concordance saved: %s\n", conc_file))
}

# ============================================================
# 8. ROC Curves
# ============================================================
cat("\n=== [8] Generating ROC curves ===\n\n")

# Combined ROC plot: overlay all cohorts
pdf(file.path(RESULTS_FIG_DIR, "same_platform_roc_curves.pdf"), width = 8, height = 7)
par(mar = c(5, 5, 4, 2))

# Initialize plot with first valid result
first_valid <- TRUE
colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
ltys <- c(1, 2, 3, 4)
legend_labs <- c()
legend_cols <- c()
legend_ltys <- c()

for (i in seq_along(names(all_results))) {
  nm <- names(all_results)[i]
  r <- all_results[[nm]]

  if (!is.na(r$auc_A) && length(unique(r$y)) == 2 && sd(r$risk_A) > 0) {
    roc_obj <- roc(r$y, r$risk_A, quiet = TRUE)

    if (first_valid) {
      plot(roc_obj, col = colors[i], lwd = 2, lty = ltys[i],
           main = "GPL570 Same-Platform Validation: ROC Curves",
           cex.main = 1.3, cex.lab = 1.2,
           xlab = "1 - Specificity", ylab = "Sensitivity")
      first_valid <- FALSE
    } else {
      plot(roc_obj, col = colors[i], lwd = 2, lty = ltys[i], add = TRUE)
    }

    legend_labs <- c(legend_labs,
                     sprintf("%s: AUC=%.3f", nm, r$auc_A))
    legend_cols <- c(legend_cols, colors[i])
    legend_ltys <- c(legend_ltys, ltys[i])
  }
}

if (length(legend_labs) > 0) {
  legend("bottomright", legend = legend_labs, col = legend_cols,
         lwd = 2, lty = legend_ltys, cex = 0.9, bty = "n")
}
dev.off()
cat("  Combined ROC curves saved\n")

# Individual ROC plots per cohort
for (nm in names(all_results)) {
  r <- all_results[[nm]]

  if (!is.na(r$auc_A) && length(unique(r$y)) == 2 && sd(r$risk_A) > 0) {
    pdf(file.path(RESULTS_FIG_DIR, paste0(nm, "_roc.pdf")), width = 7, height = 7)
    par(mar = c(5, 5, 4, 2))

    roc_A <- roc(r$y, r$risk_A, quiet = TRUE)
    roc_B <- roc(r$y, r$risk_B, quiet = TRUE)

    plot(roc_A, col = "#E41A1C", lwd = 2.5,
         main = sprintf("%s Validation (Strategy A)", nm),
         sub = sprintf("AUC = %.4f (95%% CI: %.4f - %.4f), p = %.6f",
                       r$auc_A, r$ci_A_low, r$ci_A_high, r$wt_A_p),
         cex.main = 1.2, cex.lab = 1.1,
         xlab = "1 - Specificity", ylab = "Sensitivity")

    # Add Strategy B
    plot(roc_B, col = "#377EB8", lwd = 2, lty = 2, add = TRUE)

    legend("bottomright",
           legend = c(sprintf("Strategy A (SHAP-weighted): AUC=%.3f", r$auc_A),
                      sprintf("Strategy B (mean z-score): AUC=%.3f", r$auc_B)),
           col = c("#E41A1C", "#377EB8"), lwd = 2, lty = c(1, 2), cex = 0.9, bty = "n")

    dev.off()
    cat(sprintf("  %s ROC saved\n", nm))
  }
}

# ============================================================
# 9. Direction Concordance Plot
# ============================================================
cat("\n=== [9] Generating direction concordance plots ===\n\n")

for (nm in names(all_results)) {
  r <- all_results[[nm]]
  conc <- r$concordance

  if (nrow(conc) > 0) {
    # Prepare data for plotting
    conc$Gene <- factor(conc$Gene, levels = conc$Gene[order(abs(conc$GSE39582_logFC), decreasing = TRUE)])
    conc$Concordant <- ifelse(conc$Concordant, "Concordant", "Discordant")

    # Melt logFC values for side-by-side comparison
    conc_plot <- rbind(
      data.frame(Gene = conc$Gene, logFC = conc$GSE39582_logFC,
                 Dataset = "GSE39582 (Training)",
                 Concordant = conc$Concordant, stringsAsFactors = FALSE),
      data.frame(Gene = conc$Gene, logFC = conc$Validation_logFC,
                 Dataset = nm,
                 Concordant = conc$Concordant, stringsAsFactors = FALSE)
    )
    conc_plot$Dataset <- factor(conc_plot$Dataset, levels = c("GSE39582 (Training)", nm))

    p <- ggplot(conc_plot, aes(x = Gene, y = logFC, fill = Dataset)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.7),
               width = 0.6, alpha = 0.85) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", size = 0.8) +
      scale_fill_manual(values = c("GSE39582 (Training)" = "#E41A1C", nm = "#377EB8")) +
      labs(title = paste("Directional Concordance:", nm),
           subtitle = sprintf("Concordant: %d/%d (%.1f%%), Binomial p = %.4f",
                              r$n_concordant, r$n_total_concordance,
                              ifelse(r$n_total_concordance > 0,
                                     100 * r$n_concordant / r$n_total_concordance, 0),
                              r$binom_p),
           x = "Fingerprint Gene", y = "logFC (Resistant vs Sensitive)") +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom",
            legend.title = element_blank())

    ggsave(file.path(RESULTS_FIG_DIR, paste0(nm, "_direction_concordance.pdf")),
           p, width = 9, height = 6)
    cat(sprintf("  %s direction concordance plot saved\n", nm))
  }
}

# ============================================================
# 10. Decision Summary
# ============================================================
cat("\n\n")
cat("============================================================\n")
cat("  DECISION SUMMARY\n")
cat("============================================================\n\n")

for (nm in names(all_results)) {
  r <- all_results[[nm]]
  cat(sprintf("--- %s ---\n", nm))
  cat(sprintf("  Strategy A (SHAP-weighted) AUC: %.4f (%.4f - %.4f)\n",
              r$auc_A, r$ci_A_low, r$ci_A_high))
  cat(sprintf("  Strategy B (mean z-score) AUC: %.4f (%.4f - %.4f)\n",
              r$auc_B, r$ci_B_low, r$ci_B_high))
  cat(sprintf("  Directional concordance: %d/%d (%.1f%%)\n",
              r$n_concordant, r$n_total_concordance, r$concordance_pct))

  # Decision tree
  if (!is.na(r$auc_A)) {
    if (r$auc_A >= 0.65) {
      cat("  >>> DECISION: ✅ Fingerprint VALID in this cohort\n")
    } else if (r$auc_A >= 0.55) {
      cat("  >>> DECISION: ⚠️ Borderline - proceed to GSVA + PRS\n")
    } else {
      cat("  >>> DECISION: ❌ Fingerprint FAILS - proceed to LASSO PRS\n")
    }
  } else {
    cat("  >>> DECISION: Unable to compute AUC\n")
  }
  cat("\n")
}

cat("=== Script complete! ===\n")
cat(sprintf("  Summary table: %s\n",
            file.path(RESULTS_TAB_DIR, "same_platform_validation_summary.csv")))
cat(sprintf("  ROC figures: %s\n",
            file.path(RESULTS_FIG_DIR, "same_platform_roc_curves.pdf")))
