# ============================================================
# 03_DEG_ANALYSIS.R
# ------------------------------------------------------------
# STATUS (2026-09-24): superseded for the reported numbers.
#   * Line ~64 asserted "Microarray data is already log2-transformed"; that is
#     false for GSE28702 and GSE69657, whose series matrices are on a linear
#     scale. limma was therefore fitted on mixed scales in this run.
#   * The IQR probe filter below (added 2026-06-02) was never applied to the
#     tables shipped with the manuscript, and applying it leaves no gene at
#     FDR < 0.05. The manuscript states this explicitly.
# The corrected run is scripts/reanalysis/rerun_deg3.R + meta_final.R; its
# outputs are in results_tables/reanalysis_2026-09/ and they are what the
# paper reports (nine genes at Stouffer FDR < 0.05, not the 253 of this run).
# ------------------------------------------------------------
# DEG analysis - limma for microarray data across GEO datasets
# XELOX Resistance Study
# ============================================================

Sys.setenv(TMPDIR = "/tmp", TMP = "/tmp", TEMP = "/tmp")
.libPaths(c("/path/to/Rlibs", .libPaths()))

library(limma)
library(WGCNA)
library(GEOquery)

PROJECT_ROOT <- "/path/to/xelox_project"
DATA_GEO_DIR  <- file.path(PROJECT_ROOT, "data", "geo")
DATA_PROC_DIR <- file.path(PROJECT_ROOT, "data", "processed")
RESULTS_TAB_DIR <- file.path(PROJECT_ROOT, "results", "tables")
RESULTS_FIG_DIR <- file.path(PROJECT_ROOT, "results", "figures")

dir.create(RESULTS_FIG_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== XELOX Resistance Study: DEG Analysis ===\n\n")

# ============================================================
# 1. Load clinical group data
# ============================================================
cat("[1] Loading clinical group data\n")

# Load from the clinical extraction results
groups_file <- file.path(RESULTS_TAB_DIR, "GSE39582_xelox_groups.csv")
if (!file.exists(groups_file)) {
  cat("  ERROR: group file not found. Run 02_extract_xelox_groups.R first.\n")
  quit(status = 1)
}

all_results <- list()

# ============================================================
# 2. GSE39582 DEG Analysis
# ============================================================
cat("\n[2] GSE39582 DEG analysis\n")

load_gse_expr <- function(gse_id, data_dir) {
  # Try loading full eset first, then expression matrix
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

analyze_deg_limma <- function(expr_mat, group_vec, dataset_name) {
  # Ensure group is factor with two levels
  group_f <- factor(group_vec, levels = c("sensitive", "resistant"))
  design <- model.matrix(~ group_f)
  
  # Filter low-expression genes (v6.0 fix: IQR-based filter for microarray)
  # Microarray data is already log2-transformed (range ~4-16)
  # Use IQR filter: keep probes with median above overall median AND IQR > 0.5
  probe_iqr <- apply(expr_mat, 1, IQR, na.rm = TRUE)
  probe_med <- apply(expr_mat, 1, median, na.rm = TRUE)
  overall_med <- median(expr_mat, na.rm = TRUE)
  keep <- probe_med >= overall_med & probe_iqr >= 0.5
  cat(sprintf("  Genes before filtering: %d\n", nrow(expr_mat)))
  cat(sprintf("  Genes after filtering (median>=%.1f & IQR>=0.5): %d\n",
              overall_med, sum(keep)))
  cat(sprintf("  Filtered out: %d (%.1f%%)\n", sum(!keep), sum(!keep)/nrow(expr_mat)*100))
  expr_f <- expr_mat[keep, , drop = FALSE]
  
  # Normalize (quantile)
  expr_norm <- normalizeBetweenArrays(as.matrix(expr_f), method = "quantile")
  
  # Fit linear model
  fit <- lmFit(expr_norm, design)
  fit <- eBayes(fit)
  
  # Extract results
  res <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  res$Gene <- rownames(res)
  
  # Add direction
  res$Direction <- ifelse(res$logFC > 0, "up", "down")
  
  # Annotate significance
  res$Significance <- "NS"
  res$Significance[res$P.Value < 0.05] <- "nominal"
  res$Significance[res$adj.P.Val < 0.05] <- "significant"
  
  cat(sprintf("  DEG (adj.P<0.05): %d up, %d down\n",
              sum(res$adj.P.Val < 0.05 & res$logFC > 0),
              sum(res$adj.P.Val < 0.05 & res$logFC < 0)))
  cat(sprintf("  DEG (P<0.05): %d up, %d down\n",
              sum(res$P.Value < 0.05 & res$logFC > 0),
              sum(res$P.Value < 0.05 & res$logFC < 0)))
  
  return(list(results = res, fit = fit))
}

# --- GSE39582 ---
gse39582_eset <- readRDS(file.path(DATA_GEO_DIR, "GSE39582_eset.rds"))
gse39582_expr <- exprs(gse39582_eset)

groups39582 <- read.csv(groups_file, stringsAsFactors = FALSE)

# Map groups to expression samples
samples_deg <- intersect(groups39582$sample_id, colnames(gse39582_expr))
groups_deg <- groups39582[match(samples_deg, groups39582$sample_id), ]

# Only keep resistant and sensitive
gs <- groups_deg[groups_deg$group %in% c("resistant", "sensitive"), ]
cat(sprintf("  GSE39582: %d resistant, %d sensitive\n", 
            sum(gs$group == "resistant"), sum(gs$group == "sensitive")))

if (nrow(gs) >= 6) {
  expr_sub <- gse39582_expr[, gs$sample_id, drop = FALSE]
  deg39582 <- analyze_deg_limma(expr_sub, gs$group, "GSE39582")
  deg39582$results$Dataset <- "GSE39582"
  all_results[["GSE39582"]] <- deg39582
  
  write.csv(deg39582$results, 
            file.path(RESULTS_TAB_DIR, "DEG_GSE39582_limma.csv"), row.names = FALSE)
  cat("  Results saved\n")
}

# ============================================================
# 3. GSE69657 - Pure XELOX (n=30)
# ============================================================
cat("\n[3] GSE69657 - Pure XELOX DEG analysis\n")

pheno69657 <- read.csv(file.path(RESULTS_TAB_DIR, "GSE69657_clinical_data.csv"),
                        row.names = 1, check.names = FALSE)

expr69657 <- load_gse_expr("GSE69657", DATA_GEO_DIR)
if (!is.null(expr69657)) {
  # Response: "responder" vs "noresponder"
  resp_col <- "chemoresponse:ch1"
  if (resp_col %in% colnames(pheno69657)) {
    gs69657 <- data.frame(
      sample_id = rownames(pheno69657),
      group = ifelse(pheno69657[[resp_col]] %in% c("responder", "R"), "sensitive", "resistant"),
      stringsAsFactors = FALSE
    )
    gs69657 <- gs69657[gs69657$group %in% c("resistant", "sensitive"), ]
    
    samples_deg <- intersect(gs69657$sample_id, colnames(expr69657))
    gs69657 <- gs69657[match(samples_deg, gs69657$sample_id), ]
    
    cat(sprintf("  Samples: %d resistant, %d sensitive\n",
                sum(gs69657$group == "resistant"), sum(gs69657$group == "sensitive")))
    
    if (nrow(gs69657) >= 6) {
      expr_sub <- expr69657[, gs69657$sample_id, drop = FALSE]
      deg69657 <- analyze_deg_limma(expr_sub, gs69657$group, "GSE69657")
      deg69657$results$Dataset <- "GSE69657"
      all_results[["GSE69657"]] <- deg69657
      
      write.csv(deg69657$results,
                file.path(RESULTS_TAB_DIR, "DEG_GSE69657_limma.csv"), row.names = FALSE)
      cat("  Results saved\n")
    }
  } else {
    cat("  WARNING: response column not found\n")
  }
} else {
  cat("  WARNING: expression data not found\n")
}

# ============================================================
# 4. GSE28702 - mFOLFOX6 (n=83) - clean binary labels
# ============================================================
cat("\n[4] GSE28702 - mFOLFOX6 DEG analysis\n")

pheno28702 <- read.csv(file.path(RESULTS_TAB_DIR, "GSE28702_clinical_data.csv"),
                        row.names = 1, check.names = FALSE)

expr28702 <- load_gse_expr("GSE28702", DATA_GEO_DIR)
if (!is.null(expr28702)) {
  resp_col <- "mfolfox6:ch1"
  if (resp_col %in% colnames(pheno28702)) {
    gs28702 <- data.frame(
      sample_id = rownames(pheno28702),
      group = ifelse(pheno28702[[resp_col]] == "responder", "sensitive", "resistant"),
      stringsAsFactors = FALSE
    )
    gs28702 <- gs28702[gs28702$group %in% c("resistant", "sensitive"), ]
    
    samples_deg <- intersect(gs28702$sample_id, colnames(expr28702))
    gs28702 <- gs28702[match(samples_deg, gs28702$sample_id), ]
    
    cat(sprintf("  Samples: %d resistant, %d sensitive\n",
                sum(gs28702$group == "resistant"), sum(gs28702$group == "sensitive")))
    
    if (nrow(gs28702) >= 6) {
      expr_sub <- expr28702[, gs28702$sample_id, drop = FALSE]
      deg28702 <- analyze_deg_limma(expr_sub, gs28702$group, "GSE28702")
      deg28702$results$Dataset <- "GSE28702"
      all_results[["GSE28702"]] <- deg28702
      
      write.csv(deg28702$results,
                file.path(RESULTS_TAB_DIR, "DEG_GSE28702_limma.csv"), row.names = FALSE)
      cat("  Results saved\n")
    }
  }
}

# ============================================================
# 5. GSE72970 - CRC Capecitabine Response
# ============================================================
cat("\n[5] GSE72970 - Capecitabine response DEG analysis\n")

pheno72970 <- read.csv(file.path(RESULTS_TAB_DIR, "GSE72970_clinical_data.csv"),
                        row.names = 1, check.names = FALSE)

expr72970 <- load_gse_expr("GSE72970", DATA_GEO_DIR)
if (!is.null(expr72970)) {
  resp_col <- "response category:ch1"
  if (resp_col %in% colnames(pheno72970)) {
    # Only XELOX-like regimens
    xelox_flag <- grepl("XELOX|CAPOX|FOLFOX", pheno72970[["regimen:ch1"]], ignore.case = TRUE)
    
    gs72970 <- data.frame(
      sample_id = rownames(pheno72970),
      regimen = pheno72970[["regimen:ch1"]],
      response = pheno72970[[resp_col]],
      is_xelox_like = xelox_flag,
      group = ifelse(pheno72970[[resp_col]] %in% c("CR", "PR"), "sensitive",
                     ifelse(pheno72970[[resp_col]] %in% c("PD", "SD"), "resistant", NA)),
      stringsAsFactors = FALSE
    )
    
    # Filter: XELOX-like and has valid group
    gs72970 <- gs72970[gs72970$is_xelox_like & !is.na(gs72970$group), ]
    
    samples_deg <- intersect(gs72970$sample_id, colnames(expr72970))
    gs72970 <- gs72970[match(samples_deg, gs72970$sample_id), ]
    
    cat(sprintf("  XELOX-like samples: %d resistant, %d sensitive\n",
                sum(gs72970$group == "resistant"), sum(gs72970$group == "sensitive")))
    
    if (nrow(gs72970) >= 6) {
      expr_sub <- expr72970[, gs72970$sample_id, drop = FALSE]
      deg72970 <- analyze_deg_limma(expr_sub, gs72970$group, "GSE72970")
      deg72970$results$Dataset <- "GSE72970"
      all_results[["GSE72970"]] <- deg72970
      
      write.csv(deg72970$results,
                file.path(RESULTS_TAB_DIR, "DEG_GSE72970_limma.csv"), row.names = FALSE)
      cat("  Results saved\n")
    }
  }
}

# ============================================================
# 6. GSE104645 - 1st-line Chemotherapy Response
# ============================================================
cat("\n[6] GSE104645 - 1st-line chemo response DEG analysis\n")

pheno104645 <- read.csv(file.path(RESULTS_TAB_DIR, "GSE104645_clinical_data.csv"),
                         row.names = 1, check.names = FALSE)

expr104645 <- load_gse_expr("GSE104645", DATA_GEO_DIR)
if (!is.null(expr104645)) {
  regimen_col <- "1st-line chemotherapy regimens:ch1"
  resp_col <- "best response of 1st-line chemotherapy:ch1"
  
  if (all(c(regimen_col, resp_col) %in% colnames(pheno104645))) {
    xelox_flag <- grepl("XELOX|FOLFOX|SOX|CAPOX", pheno104645[[regimen_col]], ignore.case = TRUE)
    
    gs104645 <- data.frame(
      sample_id = rownames(pheno104645),
      regimen = pheno104645[[regimen_col]],
      response = pheno104645[[resp_col]],
      is_xelox_like = xelox_flag,
      group = ifelse(grepl("Complete|Partial", pheno104645[[resp_col]]), "sensitive",
                     ifelse(grepl("Progressive|Stable", pheno104645[[resp_col]]), "resistant", NA)),
      stringsAsFactors = FALSE
    )
    
    gs104645 <- gs104645[gs104645$is_xelox_like & !is.na(gs104645$group), ]
    
    samples_deg <- intersect(gs104645$sample_id, colnames(expr104645))
    gs104645 <- gs104645[match(samples_deg, gs104645$sample_id), ]
    
    cat(sprintf("  Samples: %d resistant, %d sensitive\n",
                sum(gs104645$group == "resistant"), sum(gs104645$group == "sensitive")))
    
    if (nrow(gs104645) >= 6) {
      expr_sub <- expr104645[, gs104645$sample_id, drop = FALSE]
      deg104645 <- analyze_deg_limma(expr_sub, gs104645$group, "GSE104645")
      deg104645$results$Dataset <- "GSE104645"
      all_results[["GSE104645"]] <- deg104645
      
      write.csv(deg104645$results,
                file.path(RESULTS_TAB_DIR, "DEG_GSE104645_limma.csv"), row.names = FALSE)
      cat("  Results saved\n")
    }
  }
}

# ============================================================
# 7. Meta-analysis: Combine DEG p-values across datasets
# ============================================================
cat("\n[7] Meta-analysis of DEG results\n")

if (length(all_results) >= 2) {
  # Extract common genes across datasets
  all_genes <- lapply(all_results, function(x) rownames(x$results))
  common_genes <- Reduce(intersect, all_genes)
  
  if (length(common_genes) >= 10) {
    cat(sprintf("  Common genes across %d datasets: %d\n", length(all_results), length(common_genes)))
    
    # Fisher's method for combining p-values
    combined_p <- sapply(common_genes, function(gene) {
      p_vals <- sapply(all_results, function(x) {
        idx <- which(rownames(x$results) == gene)
        if (length(idx) > 0) x$results$P.Value[idx[1]] else 1
      })
      # Fisher's method: -2*sum(log(p)) ~ chi-squared(2k)
      stat <- -2 * sum(log(p_vals), na.rm = TRUE)
      df <- 2 * sum(!is.na(p_vals))
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
      Gene = common_genes,
      avg_logFC = avg_logFC,
      combined_p = combined_p,
      adj_p = p.adjust(combined_p, method = "BH"),
      stringsAsFactors = FALSE
    )
    meta_df$Direction <- ifelse(meta_df$avg_logFC > 0, "up", "down")
    meta_df <- meta_df[order(meta_df$combined_p), ]
    
    cat(sprintf("  Meta DEGs (adj.p<0.05): %d\n", sum(meta_df$adj_p < 0.05)))
    write.csv(meta_df, file.path(RESULTS_TAB_DIR, "DEG_meta_analysis.csv"), row.names = FALSE)
    cat("  Meta-analysis results saved\n")
  } else {
    cat("  WARNING: too few common genes for meta-analysis\n")
  }
}

# ============================================================
# Summary
# ============================================================
cat("\n=== DEG Analysis Complete ===\n")
for (ds in names(all_results)) {
  r <- all_results[[ds]]$results
  sig <- sum(r$adj.P.Val < 0.05)
  nom <- sum(r$P.Value < 0.05)
  cat(sprintf("  %s: %d significant, %d nominal DEGs\n", ds, sig, nom))
}

save.image(file.path(DATA_PROC_DIR, "03_deg_analysis.RData"))
cat("Workspace saved\n")
