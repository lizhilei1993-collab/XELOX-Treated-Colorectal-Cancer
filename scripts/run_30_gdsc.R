# ============================================================
# run_30_gdsc.R - GDSC Drug Sensitivity Validation (v6.0 fix)
# Core fix: Weighted Z-score + MAD weighting + median centering
# Expanded drug signatures
# ============================================================

Sys.setenv(TMPDIR = "C:/temp", TMP = "C:/temp", TEMP = "C:/temp")
.libPaths(c("C:/Rlibs", .libPaths()))
if (.Platform$OS.type == "windows") {
  tryCatch(Sys.setlocale("LC_ALL", "Chinese"), error = function(e) {
    tryCatch(Sys.setlocale("LC_ALL", "chs"), error = function(e2) {})
  })
}

set.seed(42)
READ_ROOT <- "C:/xelox_research"
READ_TAB_DIR <- file.path(READ_ROOT, "results", "tables")

cat("=== GDSC Drug Sensitivity Validation (v6.0 wrapper) ===\n\n")

# ============================================================
# 1. Load data
# ============================================================
cat("[1] Loading data\n")

gene_expr <- readRDS(file.path(READ_TAB_DIR, "pathway_activity", "GSE39582_gene_expression.rds"))
cat(sprintf("  Expression: %d genes x %d samples\n", nrow(gene_expr), ncol(gene_expr)))

risk_path <- file.path(READ_TAB_DIR, "nomogram", "risk_scores.csv")
if (file.exists(risk_path)) {
  risk_scores <- read.csv(risk_path, stringsAsFactors = FALSE)
  cat(sprintf("  Risk scores: %d samples\n", nrow(risk_scores)))
  # Normalize column names
  if ("risk_score" %in% names(risk_scores)) names(risk_scores)[names(risk_scores) == "risk_score"] <- "PRS"
  if ("risk_group" %in% names(risk_scores)) {
    risk_scores$RiskGroup <- ifelse(grepl("High", risk_scores$risk_group), "High", "Low")
  }
} else {
  cat("  Risk scores not found, computing from pathway scores...\n")
  pw_scores <- readRDS(file.path(READ_TAB_DIR, "pathway_activity", "GSE39582_pathway_scores.rds"))
  pathway_vars <- c("HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
                    "KEGG_ECM_RECEPTOR_INTERACTION", "KEGG_TGF_BETA_SIGNALING_PATHWAY",
                    "KEGG_PATHWAYS_IN_CANCER", "HALLMARK_MYC_TARGETS_V2",
                    "KEGG_COLORECTAL_CANCER")
  coefs <- c(1.87, 1.46, 1.53, 0.49, 1.39, 0.68, 1.32)
  prs <- rep(0, ncol(pw_scores))
  names(prs) <- colnames(pw_scores)
  for (i in seq_along(pathway_vars)) {
    if (pathway_vars[i] %in% rownames(pw_scores)) {
      prs <- prs + log(coefs[i]) * pw_scores[pathway_vars[i], ]
    }
  }
  risk_scores <- data.frame(sample_id = names(prs), PRS = prs, stringsAsFactors = FALSE)
  risk_scores$RiskGroup <- ifelse(risk_scores$PRS > median(risk_scores$PRS, na.rm = TRUE), "High", "Low")
}
cat(sprintf("  High risk: %d, Low risk: %d\n",
            sum(risk_scores$RiskGroup == "High"), sum(risk_scores$RiskGroup == "Low")))

# ============================================================
# 2. Drug signatures (v6.0: expanded)
# ============================================================
cat("\n[2] Drug signatures\n")

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

for (nm in names(drug_signatures)) {
  sig <- drug_signatures[[nm]]
  cat(sprintf("  %s: %d resistance + %d sensitivity genes\n",
              nm, length(sig$resistance), length(sig$sensitivity)))
}

# ============================================================
# 3. Compute drug scores - OLD vs NEW
# ============================================================
cat("\n[3] Computing drug scores\n")

common_samples <- intersect(risk_scores$sample_id, colnames(gene_expr))
risk_comp <- risk_scores[risk_scores$sample_id %in% common_samples, ]

# OLD: simple Z-score subtraction
compute_drug_score_OLD <- function(expr_mat, drug_sig, samples) {
  res_genes <- intersect(drug_sig$resistance, rownames(expr_mat))
  sen_genes <- intersect(drug_sig$sensitivity, rownames(expr_mat))
  res_score <- rep(0, length(samples))
  sen_score <- rep(0, length(samples))
  if (length(res_genes) > 0) {
    res_expr <- expr_mat[res_genes, samples, drop = FALSE]
    res_z <- t(scale(t(res_expr)))
    res_score <- colMeans(res_z, na.rm = TRUE)
  }
  if (length(sen_genes) > 0) {
    sen_expr <- expr_mat[sen_genes, samples, drop = FALSE]
    sen_z <- t(scale(t(sen_expr)))
    sen_score <- colMeans(sen_z, na.rm = TRUE)
  }
  score <- sen_score - res_score
  names(score) <- samples
  return(score)
}

# NEW: weighted Z-score with MAD + median centering
compute_drug_score_NEW <- function(expr_mat, drug_sig, samples) {
  res_genes <- intersect(drug_sig$resistance, rownames(expr_mat))
  sen_genes <- intersect(drug_sig$sensitivity, rownames(expr_mat))
  res_score <- rep(0, length(samples))
  sen_score <- rep(0, length(samples))
  if (length(res_genes) > 0) {
    res_expr <- expr_mat[res_genes, samples, drop = FALSE]
    res_z <- t(scale(t(res_expr)))
    res_z[is.na(res_z)] <- 0
    res_mad <- apply(res_z, 1, function(x) mad(x, na.rm = TRUE))
    res_mad[res_mad == 0 | is.na(res_mad)] <- 1
    res_w <- 1 / res_mad; res_w <- res_w / sum(res_w)
    res_score <- colSums(sweep(res_z, 1, res_w, "*"))
  }
  if (length(sen_genes) > 0) {
    sen_expr <- expr_mat[sen_genes, samples, drop = FALSE]
    sen_z <- t(scale(t(sen_expr)))
    sen_z[is.na(sen_z)] <- 0
    sen_mad <- apply(sen_z, 1, function(x) mad(x, na.rm = TRUE))
    sen_mad[sen_mad == 0 | is.na(sen_mad)] <- 1
    sen_w <- 1 / sen_mad; sen_w <- sen_w / sum(sen_w)
    sen_score <- colSums(sweep(sen_z, 1, sen_w, "*"))
  }
  score <- sen_score - res_score
  score <- score - median(score, na.rm = TRUE)
  names(score) <- samples
  return(score)
}

cat("\n  === OLD vs NEW comparison ===\n")
results <- list()

for (drug_name in names(drug_signatures)) {
  cat(sprintf("\n  --- %s ---\n", drug_name))
  
  score_old <- compute_drug_score_OLD(gene_expr, drug_signatures[[drug_name]], common_samples)
  score_new <- compute_drug_score_NEW(gene_expr, drug_signatures[[drug_name]], common_samples)
  
  high_old <- score_old[risk_comp$sample_id[risk_comp$RiskGroup == "High"]]
  low_old  <- score_old[risk_comp$sample_id[risk_comp$RiskGroup == "Low"]]
  high_new <- score_new[risk_comp$sample_id[risk_comp$RiskGroup == "High"]]
  low_new  <- score_new[risk_comp$sample_id[risk_comp$RiskGroup == "Low"]]
  
  wt_old <- wilcox.test(high_old, low_old, exact = FALSE)
  wt_new <- wilcox.test(high_new, low_new, exact = FALSE)
  
  cat(sprintf("  OLD: High=%.3f, Low=%.3f, Diff=%.3f, p=%.4f%s\n",
              mean(high_old), mean(low_old), mean(high_old)-mean(low_old),
              wt_old$p.value, ifelse(wt_old$p.value < 0.05, " *", "")))
  cat(sprintf("  NEW: High=%.3f, Low=%.3f, Diff=%.3f, p=%.4f%s\n",
              mean(high_new), mean(low_new), mean(high_new)-mean(low_new),
              wt_new$p.value, ifelse(wt_new$p.value < 0.05, " *", "")))
  cat(sprintf("  Median OLD: %.4f, NEW: %.4f\n", median(score_old), median(score_new)))
  
  results[[drug_name]] <- data.frame(
    Drug = drug_name,
    Old_High = mean(high_old), Old_Low = mean(low_old), Old_P = wt_old$p.value,
    New_High = mean(high_new), New_Low = mean(low_new), New_P = wt_new$p.value,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# 4. Summary
# ============================================================
cat("\n=== Drug Sensitivity Summary (v6.0) ===\n")
cat(sprintf("%-15s  %8s  %8s  %8s  |  %8s  %8s  %8s\n",
            "Drug", "Old_H", "Old_L", "Old_P", "New_H", "New_L", "New_P"))
cat(paste(rep("-", 75), collapse = ""), "\n")
for (nm in names(results)) {
  r <- results[[nm]]
  cat(sprintf("%-15s  %8.3f  %8.3f  %8.4f  |  %8.3f  %8.3f  %8.4f%s\n",
              r$Drug, r$Old_High, r$Old_Low, r$Old_P,
              r$New_High, r$New_Low, r$New_P,
              ifelse(r$New_P < 0.05, " *", "")))
}

cat("\n=== Script 30 Complete ===\n")
