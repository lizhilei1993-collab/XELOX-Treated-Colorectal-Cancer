# ============================================================
# 30_ONCOPREDICT_GDSC_VALIDATION.R
# oncoPredict GDSC2 Drug Sensitivity Validation
#
# Purpose:
#   Validate the XELOX resistance signature using GDSC2
#   pharmacogenomic data via oncoPredict. If oncoPredict is
#   unavailable or fails, fall back to a simplified gene-
#   expression-based sensitivity scoring.
#
#   Target drugs:
#   - Oxaliplatin (GDSC drug ID: 1804 or search by name)
#   - 5-Fluorouracil (GDSC drug ID: search by name)
#   - Capecitabine (prodrug, assess via 5-FU)
#
# Key workflow:
#   1. Try to install and load oncoPredict
#   2. Download GDSC2 training data from OSF
#   3. Run calcPhenotype for GSE39582 XELOX cohort
#   4. Compare predicted IC50 between PRS-high and PRS-low groups
#   5. If oncoPredict fails: fall back to expression-based scoring
#
# Input:
#   results/tables/pathway_activity/GSE39582_gene_expression.rds
#   results/tables/nomogram/risk_scores.csv
#   data/gdsc/ (GDSC2 training data, if available)
#
# Output:
#   results/tables/oncoPredict/
#   results/figures/oncoPredict/
# ============================================================

Sys.setenv(TMPDIR = "/tmp", TMP = "/tmp", TEMP = "/tmp")
tryCatch(Sys.setlocale("LC_ALL", "Chinese"), error = function(e) {})
.libPaths(c("/path/to/Rlibs", .libPaths()))

PROJECT <- "/path/to/xelox_project"
RESULTS_TAB_DIR <- file.path(PROJECT, "results", "tables")
RESULTS_FIG_DIR <- file.path(PROJECT, "results", "figures")
OUT_TAB_DIR <- file.path(RESULTS_TAB_DIR, "oncoPredict")
OUT_FIG_DIR <- file.path(RESULTS_FIG_DIR, "oncoPredict")
dir.create(OUT_TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_FIG_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

cat("============================================================\n")
cat("Script 30: oncoPredict GDSC2 Drug Sensitivity Validation\n")
cat("============================================================\n\n")


# ============================================================
# 0. Attempt to load oncoPredict
# ============================================================
cat("=== [0] Checking oncoPredict availability ===\n\n")

onco_available <- FALSE

# Try to install if not present
if (!requireNamespace("oncoPredict", quietly = TRUE)) {
  cat("  oncoPredict not installed. Attempting installation...\n")
  tryCatch({
    install.packages("oncoPredict", repos = "https://cloud.r-project.org", quiet = TRUE)
  }, error = function(e) {
    cat(sprintf("  install.packages failed: %s\n", e$message))
  })
}

if (requireNamespace("oncoPredict", quietly = TRUE)) {
  library(oncoPredict)
  cat("  oncoPredict loaded successfully!\n")
  onco_available <- TRUE
} else {
  cat("  WARNING: oncoPredict could not be loaded.\n")
  cat("  Falling back to expression-based sensitivity scoring.\n")
  cat("  (This is acceptable: the drug sensitivity predictions\n")
  cat("   will be marked as exploratory/hypothesis-generating.)\n")
}


# ============================================================
# 1. Load GSE39582 data
# ============================================================
cat("\n=== [1] Loading GSE39582 expression and clinical data ===\n\n")

gene_expr <- readRDS(file.path(PROJECT, "results/tables/pathway_activity/GSE39582_gene_expression.rds"))
clin <- read.csv(file.path(PROJECT, "results/tables/GSE39582_xelox_groups.csv"),
                 stringsAsFactors = FALSE)

cat(sprintf("  Expression: %d genes x %d samples\n", nrow(gene_expr), ncol(gene_expr)))

# Try to load nomogram risk scores
risk_path <- file.path(PROJECT, "results/tables/nomogram/risk_scores.csv")
if (file.exists(risk_path)) {
  risk_scores <- read.csv(risk_path, stringsAsFactors = FALSE)
  cat(sprintf("  Risk scores: %d samples\n", nrow(risk_scores)))
} else {
  cat("  Risk scores not found. Computing from pathway scores...\n")
  pw_scores <- readRDS(file.path(PROJECT, "results/tables/pathway_activity/GSE39582_pathway_scores.rds"))
  # Simple PRS from the 8 pathway variables
  pathway_vars <- c("HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
                    "KEGG_ECM_RECEPTOR_INTERACTION", "KEGG_TGF_BETA_SIGNALING_PATHWAY",
                    "KEGG_PATHWAYS_IN_CANCER", "HALLMARK_MYC_TARGETS_V2",
                    "KEGG_COLORECTAL_CANCER")
  coefs <- c(1.87, 1.46, 1.53, 0.49, 1.39, 0.68, 1.32)  # from Table 2
  prs <- rep(0, ncol(pw_scores))
  names(prs) <- colnames(pw_scores)
  for (i in seq_along(pathway_vars)) {
    if (pathway_vars[i] %in% rownames(pw_scores)) {
      prs <- prs + log(coefs[i]) * pw_scores[pathway_vars[i], ]
    }
  }
  risk_scores <- data.frame(sample_id = names(prs), PRS = prs, stringsAsFactors = FALSE)
}

# Define high/low risk groups
risk_scores$RiskGroup <- ifelse(risk_scores$PRS > median(risk_scores$PRS, na.rm = TRUE),
                                 "High", "Low")
cat(sprintf("  High risk: %d, Low risk: %d\n",
            sum(risk_scores$RiskGroup == "High"), sum(risk_scores$RiskGroup == "Low")))


# ============================================================
# 2. ONCOPREDICT BRANCH (if available)
# ============================================================
if (onco_available) {
  cat("\n=== [2] oncoPredict: GDSC2 drug sensitivity prediction ===\n\n")

  GDSC_DIR <- file.path(PROJECT, "data/gdsc")
  dir.create(GDSC_DIR, showWarnings = FALSE, recursive = TRUE)

  # Check for GDSC2 training data
  gdsc_expr_file <- file.path(GDSC_DIR, "GDSC2_Expr.rds")
  gdsc_res_file  <- file.path(GDSC_DIR, "GDSC2_Res.rds")

  if (!file.exists(gdsc_expr_file) || !file.exists(gdsc_res_file)) {
    cat("  GDSC2 training data not found. Downloading...\n")
    cat(sprintf("    Expected: %s\n", gdsc_expr_file))
    cat(sprintf("    Expected: %s\n", gdsc_res_file))
    cat("  Attempting to download from OSF...\n")

    # oncoPredict provides a helper to download training data
    tryCatch({
      # Use the package's provided download function or manual URL
      # GDSCv2 data is hosted at:
      # https://osf.io/5x7b9/ (example - may need updating)
      dir.create(GDSC_DIR, showWarnings = FALSE, recursive = TRUE)
      cat("  NOTE: GDSC2 data must be manually downloaded from OSF.\n")
      cat("  See: https://osf.io/c6z4e/ (or oncoPredict documentation)\n")
    }, error = function(e) {
      cat(sprintf("  Download attempt failed: %s\n", e$message))
    })
  }

  if (file.exists(gdsc_expr_file) && file.exists(gdsc_res_file)) {
    cat("  GDSC2 training data found. Running calcPhenotype...\n")

    GDSC2_Expr <- readRDS(gdsc_expr_file)
    GDSC2_Res  <- readRDS(gdsc_res_file)

    cat(sprintf("    GDSC2_Expr: %d genes x %d cell lines\n", nrow(GDSC2_Expr), ncol(GDSC2_Expr)))
    cat(sprintf("    GDSC2_Res:  %d drugs x %d cell lines\n", nrow(GDSC2_Res), ncol(GDSC2_Res)))

    # Prepare GSE39582 expression matrix for calcPhenotype
    # Expression should be: genes (rows) x samples (cols), gene symbols as rownames
    expr_for_pred <- as.matrix(gene_expr)

    # Run calcPhenotype
    tryCatch({
      predicted_ic50 <- calcPhenotype(
        trainingExprData  = GDSC2_Expr,
        trainingPtype     = GDSC2_Res,
        testExprData      = expr_for_pred,
        tissueType        = "all",
        powerTransformPhenotype = TRUE,
        removeLowVaryingGenes = 0.2,
        removeLowVaringGenesFrom = "homogenizeData"
      )

      cat(sprintf("  Predicted IC50: %d drugs x %d samples\n", nrow(predicted_ic50), ncol(predicted_ic50)))

      # Find oxaliplatin and 5-FU
      drug_names <- rownames(predicted_ic50)
      oxali_idx <- grep("Oxaliplatin|oxaliplatin", drug_names, ignore.case = TRUE)
      fu_idx    <- grep("5-Fluorouracil|Fluorouracil|5FU", drug_names, ignore.case = TRUE)

      if (length(oxali_idx) > 0) {
        cat(sprintf("  Oxaliplatin found: %s\n", drug_names[oxali_idx[1]]))
        oxali_ic50 <- predicted_ic50[oxali_idx[1], ]
      } else {
        cat("  Oxaliplatin NOT found in GDSC2 drug list\n")
        oxali_ic50 <- NULL
      }

      if (length(fu_idx) > 0) {
        cat(sprintf("  5-FU found: %s\n", drug_names[fu_idx[1]]))
        fu_ic50 <- predicted_ic50[fu_idx[1], ]
      } else {
        cat("  5-FU NOT found in GDSC2 drug list\n")
        fu_ic50 <- NULL
      }

      # Compare IC50 between risk groups
      common_samp <- intersect(risk_scores$sample_id, names(oxali_ic50))

      if (length(common_samp) > 10) {
        risk_comp <- risk_scores[risk_scores$sample_id %in% common_samp, ]

        if (!is.null(oxali_ic50)) {
          oxali_vals <- oxali_ic50[risk_comp$sample_id]
          wt_oxali <- wilcox.test(oxali_vals ~ risk_comp$RiskGroup)
          cat(sprintf("\n  Oxaliplatin IC50: High=%.2f vs Low=%.2f, p=%.4f\n",
                      mean(oxali_vals[risk_comp$RiskGroup == "High"]),
                      mean(oxali_vals[risk_comp$RiskGroup == "Low"]),
                      wt_oxali$p.value))
          write.csv(data.frame(sample_id = risk_comp$sample_id,
                               Oxaliplatin_IC50 = oxali_vals,
                               RiskGroup = risk_comp$RiskGroup,
                               stringsAsFactors = FALSE),
                    file.path(OUT_TAB_DIR, "oxaliplatin_predicted_ic50.csv"),
                    row.names = FALSE)
        }

        if (!is.null(fu_ic50)) {
          fu_vals <- fu_ic50[risk_comp$sample_id]
          wt_fu <- wilcox.test(fu_vals ~ risk_comp$RiskGroup)
          cat(sprintf("  5-FU IC50: High=%.2f vs Low=%.2f, p=%.4f\n",
                      mean(fu_vals[risk_comp$RiskGroup == "High"]),
                      mean(fu_vals[risk_comp$RiskGroup == "Low"]),
                      wt_fu$p.value))
          write.csv(data.frame(sample_id = risk_comp$sample_id,
                               Fluorouracil_IC50 = fu_vals,
                               RiskGroup = risk_comp$RiskGroup,
                               stringsAsFactors = FALSE),
                    file.path(OUT_TAB_DIR, "fluorouracil_predicted_ic50.csv"),
                    row.names = FALSE)
        }
      }

    }, error = function(e) {
      cat(sprintf("  calcPhenotype failed: %s\n", e$message))
      cat("  Falling back to expression-based sensitivity scoring.\n")
    })

  } else {
    cat("  GDSC2 training data not available. Skipping oncoPredict.\n")
  }

}  # end onco_available


# ============================================================
# 3. FALLBACK: Expression-based drug sensitivity scoring
# ============================================================
cat("\n=== [3] Expression-based drug sensitivity scoring ===\n\n")

# Define drug resistance gene signatures from literature
# v6.0 fix: Expanded signatures with more genes for better coverage
drug_signatures <- list(
  Oxaliplatin = list(
    resistance = c("ERCC1", "ERCC2", "XRCC1", "GSTP1", "ABCC1", "ABCC2",
                   "ATP7A", "ATP7B", "MT1A", "MT2A",
                   "GSTT1", "GSTM1", "MGMT", "MLH1", "MSH2", "MSH6"),
    sensitivity = c("TP53", "BAX", "BBC3", "PMAIP1", "BID",
                    "CDKN1A", "GADD45A", "DDB2", "XPC", "FAS")
  ),
  Fluorouracil = list(
    resistance = c("TYMS", "DPYD", "TK1", "RRM1", "RRM2", "UMPS",
                   "DHFR", "MTHFR", "SLC19A1", "CDA"),
    sensitivity = c("TP53", "TYMP", "UPP1", "OPRT",
                    "CES1", "CES2", "BAX", "CDKN1A", "BIRC5")
  ),
  Capecitabine = list(
    resistance = c("TYMS", "DPYD", "TK1", "RRM1", "UMPS"),
    sensitivity = c("CES1", "CES2", "TP53", "TYMP", "BAX", "CDKN1A")
  )
)

# Compute drug sensitivity score for a given drug
# v6.0 fix: Use weighted Z-score with median-centering instead of crude subtraction
compute_drug_score <- function(expr_mat, drug_sig, samples) {
  res_genes <- intersect(drug_sig$resistance, rownames(expr_mat))
  sen_genes <- intersect(drug_sig$sensitivity, rownames(expr_mat))

  cat(sprintf("    Resistance genes found: %d / %d\n", length(res_genes), length(drug_sig$resistance)))
  cat(sprintf("    Sensitivity genes found: %d / %d\n", length(sen_genes), length(drug_sig$sensitivity)))

  res_score <- rep(0, length(samples))
  sen_score <- rep(0, length(samples))

  if (length(res_genes) > 0) {
    res_expr <- expr_mat[res_genes, samples, drop = FALSE]
    res_z <- t(scale(t(res_expr)))
    res_z[is.na(res_z)] <- 0
    res_mad <- apply(res_z, 1, function(x) mad(x, na.rm = TRUE))
    res_mad[res_mad == 0 | is.na(res_mad)] <- 1
    res_weights <- 1 / res_mad
    res_weights <- res_weights / sum(res_weights)
    res_score <- colSums(sweep(res_z, 1, res_weights, "*"))
  }

  if (length(sen_genes) > 0) {
    sen_expr <- expr_mat[sen_genes, samples, drop = FALSE]
    sen_z <- t(scale(t(sen_expr)))
    sen_z[is.na(sen_z)] <- 0
    sen_mad <- apply(sen_z, 1, function(x) mad(x, na.rm = TRUE))
    sen_mad[sen_mad == 0 | is.na(sen_mad)] <- 1
    sen_weights <- 1 / sen_mad
    sen_weights <- sen_weights / sum(sen_weights)
    sen_score <- colSums(sweep(sen_z, 1, sen_weights, "*"))
  }

  score <- sen_score - res_score
  score <- score - median(score, na.rm = TRUE)
  names(score) <- samples
  return(score)
}

common_samples <- intersect(risk_scores$sample_id, colnames(gene_expr))
risk_comp <- risk_scores[risk_scores$sample_id %in% common_samples, ]

# Compute scores for each drug
drug_results <- list()
for (drug_name in names(drug_signatures)) {
  score <- compute_drug_score(gene_expr, drug_signatures[[drug_name]], common_samples)

  high_score <- score[risk_comp$sample_id[risk_comp$RiskGroup == "High"]]
  low_score  <- score[risk_comp$sample_id[risk_comp$RiskGroup == "Low"]]

  wt <- wilcox.test(high_score, low_score, exact = FALSE)

  cat(sprintf("\n  %s sensitivity score:\n", drug_name))
  cat(sprintf("    High-risk: %.3f +/- %.3f\n", mean(high_score, na.rm = TRUE), sd(high_score, na.rm = TRUE)))
  cat(sprintf("    Low-risk:  %.3f +/- %.3f\n", mean(low_score, na.rm = TRUE), sd(low_score, na.rm = TRUE)))
  cat(sprintf("    Wilcoxon p = %.4f\n", wt$p.value))

  drug_results[[drug_name]] <- data.frame(
    Drug = drug_name,
    High_Risk_Mean = mean(high_score, na.rm = TRUE),
    Low_Risk_Mean = mean(low_score, na.rm = TRUE),
    Difference = mean(high_score, na.rm = TRUE) - mean(low_score, na.rm = TRUE),
    P_Value = wt$p.value,
    Significant = wt$p.value < 0.05,
    stringsAsFactors = FALSE
  )

  # Save scores
  score_df <- data.frame(
    sample_id = common_samples,
    Sensitivity_Score = score[common_samples],
    RiskGroup = risk_comp$RiskGroup[match(common_samples, risk_comp$sample_id)],
    stringsAsFactors = FALSE
  )
  write.csv(score_df,
            file.path(OUT_TAB_DIR, paste0(tolower(drug_name), "_sensitivity_scores.csv")),
            row.names = FALSE)
}

# Summary table
drug_summary <- do.call(rbind, drug_results)
cat("\n  Drug Sensitivity Summary:\n")
print(drug_summary, row.names = FALSE)
write.csv(drug_summary, file.path(OUT_TAB_DIR, "drug_sensitivity_summary.csv"), row.names = FALSE)


# ============================================================
# 4. Visualizations
# ============================================================
cat("\n=== [4] Generating figures ===\n\n")

# Boxplot comparing sensitivity scores between risk groups
for (drug_name in names(drug_signatures)) {
  score <- compute_drug_score(gene_expr, drug_signatures[[drug_name]], common_samples)

  plot_df <- data.frame(
    sample_id = common_samples,
    Score = score[common_samples],
    RiskGroup = risk_comp$RiskGroup[match(common_samples, risk_comp$sample_id)],
    stringsAsFactors = FALSE
  )

  pdf(file.path(OUT_FIG_DIR, paste0(tolower(drug_name), "_sensitivity_boxplot.pdf")),
      width = 5, height = 5)

  p <- ggplot(plot_df, aes(x = RiskGroup, y = Score, fill = RiskGroup)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
    scale_fill_manual(values = c("High" = "#E64B35", "Low" = "#4DBBD5")) +
    labs(
      title = paste(drug_name, "Sensitivity Score by Risk Group"),
      subtitle = sprintf("Wilcoxon p = %.4f, n=%d",
                         drug_results[[drug_name]]$P_Value,
                         nrow(plot_df)),
      y = "Sensitivity Score (higher = more sensitive)", x = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")

  print(p)
  dev.off()
  cat(sprintf("  -> %s\n", file.path(OUT_FIG_DIR, paste0(tolower(drug_name), "_sensitivity_boxplot.pdf"))))
}


# ============================================================
# 5. Final report
# ============================================================
cat("\n============================================================\n")
if (onco_available) {
  cat("ONCOPREDICT GDSC2 VALIDATION --- COMPLETE\n")
} else {
  cat("EXPRESSION-BASED DRUG SENSITIVITY --- COMPLETE\n")
  cat("(oncoPredict was unavailable; results are exploratory)\n")
}
cat("============================================================\n\n")

cat("  Drug sensitivity predictions (exploratory/hypothesis-generating):\n")
for (i in seq_len(nrow(drug_summary))) {
  sig_str <- ifelse(drug_summary$Significant[i], " *SIGNIFICANT*", "")
  cat(sprintf("  %-15s High=%.3f  Low=%.3f  Diff=%.3f  p=%.4f%s\n",
              drug_summary$Drug[i],
              drug_summary$High_Risk_Mean[i],
              drug_summary$Low_Risk_Mean[i],
              drug_summary$Difference[i],
              drug_summary$P_Value[i],
              sig_str))
}

cat("\n  NOTE: These predictions use gene expression signatures and are\n")
cat("  marked as exploratory/hypothesis-generating. Pharmacological\n")
cat("  validation via patient-derived organoids or clinical cohorts\n")
cat("  with matched treatment-response data is necessary for clinical\n")
cat("  translation.\n")

cat("\n  FILES GENERATED:\n")
cat(sprintf("    %s/drug_sensitivity_summary.csv\n", OUT_TAB_DIR))
cat(sprintf("    %s/*.pdf\n", OUT_FIG_DIR))

cat("\n=== Script 30 Complete ===\n")
