# ============================================================
# run_03_deg.R  -  Wrapper for 03_deg_analysis.R
# Reads data from C:/xelox_research (ASCII-safe symlink)
# Writes output to E:/xelox_tmp  (sandbox-safe path)
# ============================================================

# --- Environment setup ---
Sys.setenv(TMPDIR = "C:/temp", TMP = "C:/temp", TEMP = "C:/temp")
.libPaths(c("C:/Rlibs", .libPaths()))

# Set locale so R can handle CJK paths on Windows codepage 936
if (.Platform$OS.type == "windows") {
  tryCatch(Sys.setlocale("LC_ALL", "Chinese"), error = function(e) {
    tryCatch(Sys.setlocale("LC_ALL", "chs"), error = function(e2) {
      tryCatch(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.936"),
               error = function(e3) {})
    })
  })
}
cat(sprintf("Locale: %s\n", Sys.getlocale("LC_ALL")))

# --- Load libraries ---
library(limma)
library(WGCNA)
library(GEOquery)

# --- Paths ---
# READ from the symlink (ASCII-safe, junction for data)
READ_ROOT  <- "C:/xelox_research"
# WRITE to a temp location the sandbox does not block
WRITE_ROOT <- "C:/xelox_research/tmp_v6fix"
dir.create(WRITE_ROOT, showWarnings = FALSE, recursive = TRUE)

DATA_GEO_DIR  <- file.path(READ_ROOT,  "data", "geo")
DATA_PROC_DIR <- file.path(READ_ROOT,  "data", "processed")
# Read clinical data from results/tables on C:
READ_TAB_DIR  <- file.path(READ_ROOT,  "results", "tables")
READ_FIG_DIR  <- file.path(READ_ROOT,  "results", "figures")
# Write results to E:/xelox_tmp
RESULTS_TAB_DIR <- file.path(WRITE_ROOT, "results", "tables")
RESULTS_FIG_DIR <- file.path(WRITE_ROOT, "results", "figures")
dir.create(RESULTS_TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_FIG_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== XELOX Resistance Study: DEG Analysis (v6.0 wrapper) ===\n\n")
cat(sprintf("  READ_ROOT : %s  (exists=%s)\n", READ_ROOT, file.exists(READ_ROOT)))
cat(sprintf("  WRITE_ROOT: %s  (exists=%s)\n", WRITE_ROOT, file.exists(WRITE_ROOT)))
cat(sprintf("  Data dir  : %s  (exists=%s)\n", DATA_GEO_DIR, file.exists(DATA_GEO_DIR)))
cat(sprintf("  Read tabs : %s  (exists=%s)\n", READ_TAB_DIR, file.exists(READ_TAB_DIR)))
cat(sprintf("  Write tabs: %s\n", RESULTS_TAB_DIR))

# ============================================================
# 1. Load clinical group data
# ============================================================
cat("\n[1] Loading clinical group data\n")

groups_file <- file.path(READ_TAB_DIR, "GSE39582_xelox_groups.csv")
if (!file.exists(groups_file)) {
  cat("  ERROR: group file not found at:", groups_file, "\n")
  cat("  Run 02_extract_xelox_groups.R first.\n")
  quit(status = 1)
}

all_results <- list()

# Safe write function: try to write CSV, fall back to printing to stdout
safe_write_csv <- function(df, filepath, row.names = FALSE) {
  tryCatch({
    write.csv(df, filepath, row.names = row.names)
    cat(sprintf("  Saved: %s\n", filepath))
  }, error = function(e) {
    cat(sprintf("  [SANDBOX] Write blocked, printing summary to stdout\n"))
    cat(sprintf("  Would write %d rows to: %s\n", nrow(df), filepath))
    # Print first 10 rows
    if (nrow(df) > 0) {
      print(head(df, 10))
    }
  })
}

# ============================================================
# Helper: load expression matrix
# ============================================================
load_gse_expr <- function(gse_id, data_dir) {
  eset_file <- file.path(data_dir, paste0(gse_id, "_eset.rds"))
  expr_file <- file.path(data_dir, paste0(gse_id, "_expression.rds"))
  if (file.exists(expr_file)) {
    return(readRDS(expr_file))
  } else if (file.exists(eset_file)) {
    eset <- readRDS(eset_file)
    return(exprs(eset))
  } else {
    return(NULL)
  }
}

# ============================================================
# Helper: limma DEG analysis  (v6.0 STRICT filter)
#   OLD:   keep <- rowSums(expr_mat > 0) >= 3
#   NEW:   mean_expr >= log2(2)  AND  >=50% samples above threshold
# ============================================================
analyze_deg_limma <- function(expr_mat, group_vec, dataset_name) {
  group_f <- factor(group_vec, levels = c("sensitive", "resistant"))
  design  <- model.matrix(~ group_f)

  # ---- v6.0 strict filter ----
  # Microarray data is already log2-transformed (range ~4-16)
  # log2(2)=1.0 is too low — all probes pass.
  # Use IQR-based filter: keep probes with expression IQR > 0.5
  # (removes flat/low-variability probes that add noise to FDR)
  probe_iqr <- apply(expr_mat, 1, IQR, na.rm = TRUE)
  probe_med <- apply(expr_mat, 1, median, na.rm = TRUE)
  overall_med <- median(expr_mat, na.rm = TRUE)
  # Keep probes: median above overall median AND IQR > 0.5
  keep <- probe_med >= overall_med & probe_iqr >= 0.5

  cat(sprintf("  Genes before filtering : %d\n", nrow(expr_mat)))
  cat(sprintf("  Filter (median>=%.1f & IQR>=0.5): %d kept\n",
              overall_med, sum(keep)))
  cat(sprintf("  Filtered out           : %d (%.1f%%)\n",
              sum(!keep), 100 * sum(!keep) / nrow(expr_mat)))
  expr_f <- expr_mat[keep, , drop = FALSE]

  # ---- quantile normalisation ----
  expr_norm <- normalizeBetweenArrays(as.matrix(expr_f), method = "quantile")

  # ---- limma ----
  fit <- lmFit(expr_norm, design)
  fit <- eBayes(fit)

  res <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  res$Gene      <- rownames(res)
  res$Direction <- ifelse(res$logFC > 0, "up", "down")
  res$Significance <- "NS"
  res$Significance[res$P.Value   < 0.05] <- "nominal"
  res$Significance[res$adj.P.Val < 0.05] <- "significant"

  cat(sprintf("  DEG (adj.P<0.05): %d up, %d down\n",
              sum(res$adj.P.Val < 0.05 & res$logFC > 0),
              sum(res$adj.P.Val < 0.05 & res$logFC < 0)))
  cat(sprintf("  DEG (P<0.05)   : %d up, %d down\n",
              sum(res$P.Value < 0.05 & res$logFC > 0),
              sum(res$P.Value < 0.05 & res$logFC < 0)))

  return(list(results = res, fit = fit))
}

# ============================================================
# 2. GSE39582 DEG Analysis
# ============================================================
cat("\n[2] GSE39582 DEG analysis\n")

gse39582_eset <- readRDS(file.path(DATA_GEO_DIR, "GSE39582_eset.rds"))
gse39582_expr <- exprs(gse39582_eset)

groups39582 <- read.csv(groups_file, stringsAsFactors = FALSE)

samples_deg <- intersect(groups39582$sample_id, colnames(gse39582_expr))
groups_deg  <- groups39582[match(samples_deg, groups39582$sample_id), ]

gs <- groups_deg[groups_deg$group %in% c("resistant", "sensitive"), ]
cat(sprintf("  GSE39582: %d resistant, %d sensitive\n",
            sum(gs$group == "resistant"), sum(gs$group == "sensitive")))

if (nrow(gs) >= 6) {
  expr_sub  <- gse39582_expr[, gs$sample_id, drop = FALSE]
  deg39582  <- analyze_deg_limma(expr_sub, gs$group, "GSE39582")
  deg39582$results$Dataset <- "GSE39582"
  all_results[["GSE39582"]] <- deg39582

  safe_write_csv(deg39582$results,
                 file.path(RESULTS_TAB_DIR, "DEG_GSE39582_limma.csv"))
}

# ============================================================
# 3. GSE69657 - Pure XELOX (n=30)
# ============================================================
cat("\n[3] GSE69657 - Pure XELOX DEG analysis\n")

pheno69657 <- read.csv(file.path(READ_TAB_DIR, "GSE69657_clinical_data.csv"),
                       row.names = 1, check.names = FALSE)

expr69657 <- load_gse_expr("GSE69657", DATA_GEO_DIR)
if (!is.null(expr69657)) {
  resp_col <- "chemoresponse:ch1"
  if (resp_col %in% colnames(pheno69657)) {
    gs69657 <- data.frame(
      sample_id = rownames(pheno69657),
      group = ifelse(pheno69657[[resp_col]] %in% c("responder", "R"),
                     "sensitive", "resistant"),
      stringsAsFactors = FALSE
    )
    gs69657 <- gs69657[gs69657$group %in% c("resistant", "sensitive"), ]
    samples_deg <- intersect(gs69657$sample_id, colnames(expr69657))
    gs69657     <- gs69657[match(samples_deg, gs69657$sample_id), ]

    cat(sprintf("  Samples: %d resistant, %d sensitive\n",
                sum(gs69657$group == "resistant"),
                sum(gs69657$group == "sensitive")))

    if (nrow(gs69657) >= 6) {
      expr_sub <- expr69657[, gs69657$sample_id, drop = FALSE]
      deg69657 <- analyze_deg_limma(expr_sub, gs69657$group, "GSE69657")
      deg69657$results$Dataset <- "GSE69657"
      all_results[["GSE69657"]] <- deg69657

      safe_write_csv(deg69657$results,
                     file.path(RESULTS_TAB_DIR, "DEG_GSE69657_limma.csv"))
    }
  } else {
    cat("  WARNING: response column not found\n")
  }
} else {
  cat("  WARNING: expression data not found\n")
}

# ============================================================
# 4. GSE28702 - mFOLFOX6 (n=83)
# ============================================================
cat("\n[4] GSE28702 - mFOLFOX6 DEG analysis\n")

pheno28702 <- read.csv(file.path(READ_TAB_DIR, "GSE28702_clinical_data.csv"),
                       row.names = 1, check.names = FALSE)

expr28702 <- load_gse_expr("GSE28702", DATA_GEO_DIR)
if (!is.null(expr28702)) {
  resp_col <- "mfolfox6:ch1"
  if (resp_col %in% colnames(pheno28702)) {
    gs28702 <- data.frame(
      sample_id = rownames(pheno28702),
      group = ifelse(pheno28702[[resp_col]] == "responder",
                     "sensitive", "resistant"),
      stringsAsFactors = FALSE
    )
    gs28702 <- gs28702[gs28702$group %in% c("resistant", "sensitive"), ]
    samples_deg <- intersect(gs28702$sample_id, colnames(expr28702))
    gs28702     <- gs28702[match(samples_deg, gs28702$sample_id), ]

    cat(sprintf("  Samples: %d resistant, %d sensitive\n",
                sum(gs28702$group == "resistant"),
                sum(gs28702$group == "sensitive")))

    if (nrow(gs28702) >= 6) {
      expr_sub <- expr28702[, gs28702$sample_id, drop = FALSE]
      deg28702 <- analyze_deg_limma(expr_sub, gs28702$group, "GSE28702")
      deg28702$results$Dataset <- "GSE28702"
      all_results[["GSE28702"]] <- deg28702

      safe_write_csv(deg28702$results,
                     file.path(RESULTS_TAB_DIR, "DEG_GSE28702_limma.csv"))
    }
  }
}

# ============================================================
# 5. GSE72970 - CRC Capecitabine Response
# ============================================================
cat("\n[5] GSE72970 - Capecitabine response DEG analysis\n")

pheno72970 <- read.csv(file.path(READ_TAB_DIR, "GSE72970_clinical_data.csv"),
                       row.names = 1, check.names = FALSE)

expr72970 <- load_gse_expr("GSE72970", DATA_GEO_DIR)
if (!is.null(expr72970)) {
  resp_col <- "response category:ch1"
  if (resp_col %in% colnames(pheno72970)) {
    xelox_flag <- grepl("XELOX|CAPOX|FOLFOX",
                        pheno72970[["regimen:ch1"]], ignore.case = TRUE)

    gs72970 <- data.frame(
      sample_id    = rownames(pheno72970),
      regimen      = pheno72970[["regimen:ch1"]],
      response     = pheno72970[[resp_col]],
      is_xelox_like = xelox_flag,
      group = ifelse(pheno72970[[resp_col]] %in% c("CR", "PR"), "sensitive",
                     ifelse(pheno72970[[resp_col]] %in% c("PD", "SD"), "resistant", NA)),
      stringsAsFactors = FALSE
    )
    gs72970 <- gs72970[gs72970$is_xelox_like & !is.na(gs72970$group), ]

    samples_deg <- intersect(gs72970$sample_id, colnames(expr72970))
    gs72970     <- gs72970[match(samples_deg, gs72970$sample_id), ]

    cat(sprintf("  XELOX-like samples: %d resistant, %d sensitive\n",
                sum(gs72970$group == "resistant"),
                sum(gs72970$group == "sensitive")))

    if (nrow(gs72970) >= 6) {
      expr_sub <- expr72970[, gs72970$sample_id, drop = FALSE]
      deg72970 <- analyze_deg_limma(expr_sub, gs72970$group, "GSE72970")
      deg72970$results$Dataset <- "GSE72970"
      all_results[["GSE72970"]] <- deg72970

      safe_write_csv(deg72970$results,
                     file.path(RESULTS_TAB_DIR, "DEG_GSE72970_limma.csv"))
    }
  }
}

# ============================================================
# 6. GSE104645 - 1st-line Chemotherapy Response
# ============================================================
cat("\n[6] GSE104645 - 1st-line chemo response DEG analysis\n")

pheno104645 <- read.csv(file.path(READ_TAB_DIR, "GSE104645_clinical_data.csv"),
                        row.names = 1, check.names = FALSE)

expr104645 <- load_gse_expr("GSE104645", DATA_GEO_DIR)
if (!is.null(expr104645)) {
  regimen_col <- "1st-line chemotherapy regimens:ch1"
  resp_col    <- "best response of 1st-line chemotherapy:ch1"

  if (all(c(regimen_col, resp_col) %in% colnames(pheno104645))) {
    xelox_flag <- grepl("XELOX|FOLFOX|SOX|CAPOX",
                        pheno104645[[regimen_col]], ignore.case = TRUE)

    gs104645 <- data.frame(
      sample_id    = rownames(pheno104645),
      regimen      = pheno104645[[regimen_col]],
      response     = pheno104645[[resp_col]],
      is_xelox_like = xelox_flag,
      group = ifelse(grepl("Complete|Partial",  pheno104645[[resp_col]]), "sensitive",
                     ifelse(grepl("Progressive|Stable", pheno104645[[resp_col]]), "resistant", NA)),
      stringsAsFactors = FALSE
    )
    gs104645 <- gs104645[gs104645$is_xelox_like & !is.na(gs104645$group), ]

    samples_deg <- intersect(gs104645$sample_id, colnames(expr104645))
    gs104645    <- gs104645[match(samples_deg, gs104645$sample_id), ]

    cat(sprintf("  Samples: %d resistant, %d sensitive\n",
                sum(gs104645$group == "resistant"),
                sum(gs104645$group == "sensitive")))

    if (nrow(gs104645) >= 6) {
      expr_sub  <- expr104645[, gs104645$sample_id, drop = FALSE]
      deg104645 <- analyze_deg_limma(expr_sub, gs104645$group, "GSE104645")
      deg104645$results$Dataset <- "GSE104645"
      all_results[["GSE104645"]] <- deg104645

      safe_write_csv(deg104645$results,
                     file.path(RESULTS_TAB_DIR, "DEG_GSE104645_limma.csv"))
    }
  }
}

# ============================================================
# 7. Meta-analysis: Combine DEG p-values across datasets
# ============================================================
cat("\n[7] Meta-analysis of DEG results\n")

if (length(all_results) >= 2) {
  all_genes    <- lapply(all_results, function(x) rownames(x$results))
  common_genes <- Reduce(intersect, all_genes)

  if (length(common_genes) >= 10) {
    cat(sprintf("  Common genes across %d datasets: %d\n",
                length(all_results), length(common_genes)))

    # Fisher's method for combining p-values
    combined_p <- sapply(common_genes, function(gene) {
      p_vals <- sapply(all_results, function(x) {
        idx <- which(rownames(x$results) == gene)
        if (length(idx) > 0) x$results$P.Value[idx[1]] else 1
      })
      stat <- -2 * sum(log(p_vals), na.rm = TRUE)
      df   <- 2 * sum(!is.na(p_vals))
      pchisq(stat, df = df, lower.tail = FALSE)
    })

    # Average logFC
    avg_logFC <- sapply(common_genes, function(gene) {
      lfc <- sapply(all_results, function(x) {
        idx <- which(rownames(x$results) == gene)
        if (length(idx) > 0) x$results$logFC[idx[1]] else NA
      })
      mean(lfc, na.rm = TRUE)
    })

    meta_df <- data.frame(
      Gene       = common_genes,
      avg_logFC  = avg_logFC,
      combined_p = combined_p,
      adj_p      = p.adjust(combined_p, method = "BH"),
      stringsAsFactors = FALSE
    )
    meta_df$Direction <- ifelse(meta_df$avg_logFC > 0, "up", "down")
    meta_df <- meta_df[order(meta_df$combined_p), ]

    cat(sprintf("  Meta DEGs (adj.p<0.05): %d\n", sum(meta_df$adj_p < 0.05)))
    safe_write_csv(meta_df,
                   file.path(RESULTS_TAB_DIR, "DEG_meta_analysis.csv"))
  } else {
    cat("  WARNING: too few common genes for meta-analysis\n")
  }
}

# ============================================================
# 8. Comparison: old vs new filter counts (diagnostic)
# ============================================================
cat("\n[8] Filter comparison (OLD vs v6.0 STRICT)\n")
cat("  OLD  filter: rowSums(expr > 0) >= 3\n")
cat("  NEW  filter: median >= overall_median AND IQR >= 0.5\n")

for (ds_name in names(all_results)) {
  cat(sprintf("\n  --- %s ---\n", ds_name))
  orig_expr <- switch(ds_name,
    "GSE39582" = gse39582_expr,
    "GSE69657" = expr69657,
    "GSE28702" = expr28702,
    "GSE72970" = expr72970,
    "GSE104645" = expr104645,
    NULL
  )
  if (!is.null(orig_expr)) {
    old_keep <- rowSums(orig_expr > 0, na.rm = TRUE) >= 3
    probe_iqr <- apply(orig_expr, 1, IQR, na.rm = TRUE)
    probe_med <- apply(orig_expr, 1, median, na.rm = TRUE)
    overall_med <- median(orig_expr, na.rm = TRUE)
    new_keep <- probe_med >= overall_med & probe_iqr >= 0.5
    cat(sprintf("  Total genes      : %d\n", nrow(orig_expr)))
    cat(sprintf("  OLD keep (>=3)   : %d  (%.1f%%)\n",
                sum(old_keep), 100 * sum(old_keep) / nrow(orig_expr)))
    cat(sprintf("  NEW keep (IQR)   : %d  (%.1f%%)\n",
                sum(new_keep), 100 * sum(new_keep) / nrow(orig_expr)))
    cat(sprintf("  Additional removed: %d\n",
                sum(old_keep & !new_keep)))
  }
}

# ============================================================
# Summary
# ============================================================
cat("\n=== DEG Analysis Complete (v6.0 wrapper) ===\n")
for (ds in names(all_results)) {
  r   <- all_results[[ds]]$results
  sig <- sum(r$adj.P.Val < 0.05)
  nom <- sum(r$P.Value < 0.05)
  cat(sprintf("  %s: %d significant (adj.P<0.05), %d nominal (P<0.05)\n",
              ds, sig, nom))
}

# Save workspace (may be blocked by sandbox)
tryCatch({
  save.image(file.path(WRITE_ROOT, "03_deg_analysis.RData"))
  cat(sprintf("\nWorkspace saved to: %s\n", file.path(WRITE_ROOT, "03_deg_analysis.RData")))
}, error = function(e) cat("\n  [SANDBOX] Workspace save blocked\n"))
cat(sprintf("All outputs would go to: %s\n", RESULTS_TAB_DIR))
cat("\nDone.\n")
