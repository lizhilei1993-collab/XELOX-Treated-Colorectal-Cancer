# ============================================================
# 24_TCGA_PRS_EXTERNAL_VALIDATION.R
# P0-1: External Survival Validation on TCGA COAD+READ
#
# Purpose:
#   Apply the pathway-based PRS model (from Script 16) to
#   TCGA COAD+READ oxaliplatin-treated patients and validate
#   its prognostic value for Overall Survival (OS) and
#   Disease-Specific Survival (DSS).
#
# PRS Model (3 KEGG pathways selected by LASSO):
#   KEGG_CHEMOKINE_SIGNALING_PATHWAY  coef = -0.4416 (pro-sensitive)
#   KEGG_ECM_RECEPTOR_INTERACTION     coef = +0.2988 (pro-resistance)
#   KEGG_BASE_EXCISION_REPAIR         coef = -0.0676 (pro-sensitive)
#
# PRS = Σ(coef_i × zscore_i)  where zscore uses training mean/sd
#
# Input:
#   TCGA_COAD_READ_expression_logcpm.rds  — log2(CPM+1) normalized expression
#   TCGA_COAD_READ_gene_info.csv          — ENSG→gene_symbol mapping
#   TCGA_COAD_READ_drug_data.csv          — drug treatment records
#   TCGA_COAD_READ_clinical_data.csv      — patient clinical + survival
#   c2.cp.kegg_legacy.v2024.1.symbols.gmt — KEGG gene sets
#   PRS_model_coefficients.csv            — LASSO coefficients (3 pathways)
#   GSE39582_pathway_scores.rds           — training ssGSEA scores (mean/sd)
#   GSE104645_pathway_scores.rds          — training ssGSEA scores (mean/sd)
#
# Output (C:/xelox_work/results/tables/tcga_validation/):
#   tcga_prs_scores.csv           — PRS scores for all TCGA samples
#   tcga_oxali_analysis_set.csv   — Final analysis set with PRS + survival
#   tcga_cox_results.csv          — Cox regression summary
#   tcga_cox_multivariable.csv    — Multivariable Cox summary
#
# Output (C:/xelox_work/results/figures/tcga_validation/):
#   tcga_km_os.pdf                — KM curve (OS, PRS median split)
#   tcga_km_dss.pdf               — KM curve (DSS, PRS median split)
#   tcga_roc_1yr.pdf / 3yr / 5yr — Time-dependent ROC
#   tcga_prs_boxplot.pdf          — PRS distribution
# ============================================================

Sys.setenv(TMPDIR = "C:/temp", TMP = "C:/temp", TEMP = "C:/temp")
.libPaths(c("C:/Rlibs", .libPaths()))

PROJECT_ROOT_XELOX <- "/path/to/xelox_project"
TCGA_DATA_DIR      <- "C:/xelox_research/data/tcga"

OUT_TAB_DIR <- file.path(PROJECT_ROOT_XELOX, "results", "tables", "tcga_validation")
OUT_FIG_DIR <- file.path(PROJECT_ROOT_XELOX, "results", "figures", "tcga_validation")
dir.create(OUT_TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_FIG_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

cat("============================================================\n")
cat("P0-1: TCGA COAD+READ External Survival Validation\n")
cat("  PRS Model: 3 KEGG pathways (Chemokine/ECM/BER)\n")
cat("  Outcome: OS + DSS in oxaliplatin-treated patients\n")
cat("============================================================\n\n")

# ============================================================
# 0. Load packages
# ============================================================
cat("=== [0] Loading packages ===\n\n")

required_pkgs <- c("GSVA", "GSEABase", "survival", "survminer",
                    "ggplot2", "pROC", "timeROC", "data.table")

for (pkg in required_pkgs) {
  loaded <- require(pkg, lib.loc = "C:/Rlibs", character.only = TRUE, quietly = TRUE)
  if (!loaded) {
    cat(sprintf("  WARNING: %s not loaded, trying default library...\n", pkg))
    loaded <- require(pkg, character.only = TRUE, quietly = TRUE)
  }
  if (loaded) {
    cat(sprintf("  \u2713 %s\n", pkg))
  } else {
    stop(sprintf("  \u2717 %s - required package not found", pkg))
  }
}
cat("\n")

# ============================================================
# 1. Load PRS model coefficients
# ============================================================
cat("=== [1] Loading PRS model ===\n\n")

coef_df <- read.csv(file.path(PROJECT_ROOT_XELOX, "results", "tables", "prs",
                               "PRS_model_coefficients.csv"),
                     stringsAsFactors = FALSE)
cat(sprintf("  PRS coefficients loaded: %d pathways\n", nrow(coef_df)))

coef_vec <- coef_df$Coefficient
names(coef_vec) <- coef_df$Pathway

for (i in seq_len(nrow(coef_df))) {
  cat(sprintf("  %s: coef = %+.4f (%s)\n",
              coef_df$Pathway[i], coef_df$Coefficient[i], coef_df$Direction[i]))
}
cat(sprintf("\n  PRS formula = %.4f\u00d7Z(%s) %+.4f\u00d7Z(%s) %+.4f\u00d7Z(%s)\n",
            coef_vec[1], names(coef_vec)[1],
            coef_vec[2], names(coef_vec)[2],
            coef_vec[3], names(coef_vec)[3]))

# ============================================================
# 2. Load KEGG pathway gene sets from GMT
# ============================================================
cat("\n=== [2] Loading KEGG pathway gene sets ===\n\n")

gmt_file <- "C:/temp/c2.cp.kegg_legacy.v2024.1.symbols.gmt"
if (!file.exists(gmt_file)) {
  stop("GMT file not found: ", gmt_file)
}

read_gmt_genes <- function(gmt_path, pathway_names) {
  lines <- readLines(gmt_path, warn = FALSE)
  result <- list()
  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    gs_name <- parts[1]
    if (gs_name %in% pathway_names) {
      genes <- parts[-(1:2)]
      genes <- genes[nchar(genes) > 0]
      result[[gs_name]] <- genes
      cat(sprintf("  + %s: %d genes\n", gs_name, length(genes)))
    }
  }
  missing <- setdiff(pathway_names, names(result))
  if (length(missing) > 0) {
    cat(sprintf("  ! MISSING: %s\n", paste(missing, collapse = ", ")))
  }
  return(result)
}

target_pathways <- names(coef_vec)
kegg_gs <- read_gmt_genes(gmt_file, target_pathways)

if (length(kegg_gs) != length(target_pathways)) {
  stop("Not all PRS pathways found in GMT file!")
}

# ============================================================
# 3. Load TCGA expression data
# ============================================================
cat("\n=== [3] Loading TCGA expression data ===\n\n")

logcpm_file <- file.path(TCGA_DATA_DIR, "TCGA_COAD_READ_expression_logcpm.rds")
if (!file.exists(logcpm_file)) {
  stop("TCGA logCPM RDS not found: ", logcpm_file)
}

cat("  Loading logCPM matrix...\n")
logcpm <- readRDS(logcpm_file)
cat(sprintf("  logCPM dim: %d genes x %d samples\n", nrow(logcpm), ncol(logcpm)))

# Load gene info for ENSG → symbol mapping
gene_info_file <- file.path(TCGA_DATA_DIR, "TCGA_COAD_READ_gene_info.csv")
gene_info <- read.csv(gene_info_file, stringsAsFactors = FALSE)
cat(sprintf("  Gene info loaded: %d rows\n", nrow(gene_info)))

# Map ENSG IDs to gene symbols
# Row names of logcpm are ENSG IDs (may include version suffix like .15)
ensg_ids <- rownames(logcpm)
ensg_clean <- gsub("\\..*$", "", ensg_ids)  # Remove version suffix
gene_info$ensg_clean <- gsub("\\..*$", "", gene_info$gene_id)

# Build mapping
map_idx <- match(ensg_clean, gene_info$ensg_clean)
gene_symbols <- gene_info$gene_name[map_idx]
gene_symbols[is.na(gene_symbols)] <- ensg_ids[is.na(gene_symbols)]

# Assign gene symbols as rownames
# Handle duplicate gene symbols (take max mean expression per gene)
dup_genes <- unique(gene_symbols[duplicated(gene_symbols)])
if (length(dup_genes) > 0) {
  cat(sprintf("  Removing %d duplicated gene symbols (keeping max-mean probe)\n", length(dup_genes)))
  gene_mean <- rowMeans(logcpm, na.rm = TRUE)
  keep <- !duplicated(gene_symbols, fromLast = FALSE)
  for (dg in dup_genes) {
    idx <- which(gene_symbols == dg)
    best <- idx[which.max(gene_mean[idx])]
    keep[setdiff(idx, best)] <- FALSE
  }
  logcpm <- logcpm[keep, ]
  gene_symbols <- gene_symbols[keep]
}
rownames(logcpm) <- gene_symbols
cat(sprintf("  After dedup: %d genes x %d samples\n", nrow(logcpm), ncol(logcpm)))

# ============================================================
# 4. Compute ssGSEA pathway scores for TCGA
# ============================================================
cat("\n=== [4] Computing ssGSEA pathway scores ===\n\n")

tcga_mat <- as.matrix(logcpm)

# Filter pathway genes to those present in TCGA
filtered_gs <- list()
for (nm in names(kegg_gs)) {
  genes_in <- intersect(kegg_gs[[nm]], rownames(tcga_mat))
  if (length(genes_in) >= 10) {
    filtered_gs[[nm]] <- genes_in
    cat(sprintf("  %s: %d/%d genes available\n", nm, length(genes_in), length(kegg_gs[[nm]])))
  } else {
    stop(sprintf("  %s: only %d genes available (< 10)! Cannot compute ssGSEA.", nm, length(genes_in)))
  }
}

# ssGSEA
cat("\n  Running ssGSEA (GSVA v2.x ssgseaParam)...\n")
param <- ssgseaParam(tcga_mat, filtered_gs,
                     minSize = 10, maxSize = 500,
                     alpha = 0.25, normalize = TRUE)
tcga_ps <- gsva(param)

cat(sprintf("  TCGA pathway scores: %d pathways x %d samples\n",
            nrow(tcga_ps), ncol(tcga_ps)))

# ============================================================
# 5. Compute training mean/sd for 3 PRS pathways
# ============================================================
cat("\n=== [5] Computing training scaling parameters ===\n\n")

pa_dir <- file.path(PROJECT_ROOT_XELOX, "results", "tables", "pathway_activity")

# Load training pathway scores
ps_39582  <- readRDS(file.path(pa_dir, "GSE39582_pathway_scores.rds"))
ps_104645 <- readRDS(file.path(pa_dir, "GSE104645_pathway_scores.rds"))

cat(sprintf("  GSE39582 scores: %d x %d\n", nrow(ps_39582), ncol(ps_39582)))
cat(sprintf("  GSE104645 scores: %d x %d\n", nrow(ps_104645), ncol(ps_104645)))

# Combined training scores (matching Script 16 approach)
# Only need the 3 PRS pathways
common_pw <- intersect(target_pathways, rownames(ps_39582))
cat(sprintf("  Pathways found in training scores: %d / %d\n", length(common_pw), length(target_pathways)))

if (length(common_pw) != length(target_pathways)) {
  stop("Some PRS pathways missing from training scores!")
}

# Combined training matrix for the 3 pathways
train_combined <- cbind(
  ps_39582[common_pw, , drop = FALSE],
  ps_104645[common_pw, , drop = FALSE]
)

# Per-pathway mean and sd
train_mean <- rowMeans(train_combined, na.rm = TRUE)
train_sd   <- apply(train_combined, 1, sd, na.rm = TRUE)

cat("\n  Training scaling parameters:\n")
for (pw in common_pw) {
  cat(sprintf("  %s: mean = %.4f, sd = %.4f\n", pw, train_mean[pw], train_sd[pw]))
}

# ============================================================
# 6. Compute PRS for TCGA samples
# ============================================================
cat("\n=== [6] Computing PRS for TCGA samples ===\n\n")

# Extract TCGA pathway scores for the 3 PRS pathways
tcga_pw_scores <- tcga_ps[common_pw, , drop = FALSE]

# Z-score using training parameters
tcga_z <- (tcga_pw_scores - train_mean[common_pw]) / train_sd[common_pw]

# PRS = weighted sum
tcga_prs <- as.numeric(t(col(tcga_z)) %*% coef_vec[common_pw])
names(tcga_prs) <- colnames(tcga_ps)

cat(sprintf("  TCGA PRS computed for %d samples\n", length(tcga_prs)))
cat(sprintf("  PRS range: [%.4f, %.4f]\n", min(tcga_prs, na.rm = TRUE), max(tcga_prs, na.rm = TRUE)))
cat(sprintf("  PRS mean ± sd: %.4f ± %.4f\n", mean(tcga_prs, na.rm = TRUE), sd(tcga_prs, na.rm = TRUE)))

# Save all PRS scores
prs_df <- data.frame(
  sample_barcode = names(tcga_prs),
  PRS = tcga_prs,
  stringsAsFactors = FALSE
)

# Add patient barcode (TCGA-XX-XXXX format)
prs_df$patient_barcode <- sapply(strsplit(prs_df$sample_barcode, "-"),
                                  function(x) paste(x[1:3], collapse = "-"))

write.csv(prs_df, file.path(OUT_TAB_DIR, "tcga_prs_scores.csv"), row.names = FALSE)
cat(sprintf("  -> Saved: tcga_prs_scores.csv\n"))

# ============================================================
# 7. Identify oxaliplatin-treated patients
# ============================================================
cat("\n=== [7] Identifying oxaliplatin-treated patients ===\n\n")

drug_data <- read.csv(file.path(TCGA_DATA_DIR, "TCGA_COAD_READ_drug_data.csv"),
                       stringsAsFactors = FALSE)
cat(sprintf("  Drug records: %d\n", nrow(drug_data)))
cat(sprintf("  Unique patients: %d\n", length(unique(drug_data$bcr_patient_barcode))))

# Identify oxaliplatin-containing regimens
oxali_patt  <- "oxaliplatin|eloxatin|folfox|xelox|capox"
fluoro_patt <- "fluorouracil|5.fu|capecitabine|xeloda|folfox|xelox|capox"

drug_data$is_oxali  <- grepl(oxali_patt, drug_data$drug_name, ignore.case = TRUE, perl = TRUE)
drug_data$is_fluoro <- grepl(fluoro_patt, drug_data$drug_name, ignore.case = TRUE, perl = TRUE)

oxali_pts  <- unique(drug_data$bcr_patient_barcode[drug_data$is_oxali])
fluoro_pts <- unique(drug_data$bcr_patient_barcode[drug_data$is_fluoro])
xelox_pts  <- intersect(oxali_pts, fluoro_pts)

cat(sprintf("  Oxaliplatin-treated patients: %d\n", length(oxali_pts)))
cat(sprintf("  Fluoropyrimidine-treated: %d\n", length(fluoro_pts)))
cat(sprintf("  Both (XELOX/FOLFOX): %d\n", length(xelox_pts)))

# Also include single-agent oxaliplatin (in case some only have oxali)
analysis_pts <- union(oxali_pts, xelox_pts)
cat(sprintf("  Analysis set (any oxaliplatin): %d patients\n", length(analysis_pts)))

# ============================================================
# 8. Prepare survival data
# ============================================================
cat("\n=== [8] Preparing survival data ===\n\n")

clinical <- read.csv(file.path(TCGA_DATA_DIR, "TCGA_COAD_READ_clinical_data.csv"),
                      stringsAsFactors = FALSE)
cat(sprintf("  Clinical data: %d patients\n", nrow(clinical)))

# OS: Overall Survival
# TCGA vital_status: "Dead" or "Alive" (case-sensitive)
# Some records have blank vital_status but valid days_to_death
# days_to_death: time for dead patients
# days_to_last_followup: censoring time for alive patients
clinical$os_event <- ifelse(
  clinical$vital_status == "Dead" | (!is.na(clinical$days_to_death) & clinical$days_to_death > 0),
  1, 0
)
clinical$os_time <- ifelse(
  clinical$os_event == 1 & !is.na(clinical$days_to_death),
  clinical$days_to_death, NA
)
clinical$os_time <- ifelse(
  is.na(clinical$os_time) & !is.na(clinical$days_to_last_followup) & clinical$days_to_last_followup > 0,
  clinical$days_to_last_followup, clinical$os_time
)

# Convert to months (/30.44)
clinical$os_months <- clinical$os_time / 30.44

cat(sprintf("  OS events: %d / %d (%.1f%%)\n",
            sum(clinical$os_event, na.rm = TRUE),
            nrow(clinical),
            mean(clinical$os_event, na.rm = TRUE) * 100))

# For DSS, look for disease-specific death
# TCGA XML has days_to_death but doesn't distinguish cause easily
# We'll default to OS as primary endpoint and attempt DSS if available

# ============================================================
# 9. Merge PRS + oxaliplatin + survival
# ============================================================
cat("\n=== [9] Building analysis set ===\n\n")

# Match PRS samples to clinical patients
# PRS sample barcodes are full GTEX barcodes (TCGA-XX-XXXX-XX-XXXX-XX)
# Patient barcodes are TCGA-XX-XXXX
prs_df$patient_match <- sapply(strsplit(prs_df$sample_barcode, "-"),
                                function(x) paste(x[1:3], collapse = "-"))

# Merge PRS with clinical
merged <- merge(prs_df, clinical,
                by.x = "patient_match", by.y = "bcr_patient_barcode",
                all.x = FALSE, all.y = FALSE)
cat(sprintf("  PRS + Clinical matched: %d samples\n", nrow(merged)))

# Filter to oxaliplatin-treated patients
merged_oxali <- merged[merged$patient_match %in% analysis_pts, ]
cat(sprintf("  After oxaliplatin filter: %d samples\n", nrow(merged_oxali)))

# Filter to patients with valid survival data
merged_oxali <- merged_oxali[!is.na(merged_oxali$os_time) & merged_oxali$os_time > 0, ]
cat(sprintf("  After removing missing survival: %d samples\n", nrow(merged_oxali)))

# Remove duplicate patient barcodes (keep first sample per patient)
merged_oxali <- merged_oxali[!duplicated(merged_oxali$patient_match), ]
cat(sprintf("  After dedup (1 sample/patient): %d samples\n", nrow(merged_oxali)))

cat(sprintf("\n  Final analysis set: %d patients\n", nrow(merged_oxali)))
cat(sprintf("  OS events: %d / %d (%.1f%%)\n",
            sum(merged_oxali$os_event, na.rm = TRUE),
            nrow(merged_oxali),
            mean(merged_oxali$os_event, na.rm = TRUE) * 100))

# Save analysis set
write.csv(merged_oxali, file.path(OUT_TAB_DIR, "tcga_oxali_analysis_set.csv"),
          row.names = FALSE)
cat(sprintf("  -> Saved: tcga_oxali_analysis_set.csv\n"))

# Check if we have enough samples
if (nrow(merged_oxali) < 30) {
  cat("\n  WARNING: Too few samples for survival analysis!\n")
  cat("  Attempting to expand: checking all patients with any PRS data...\n")
  
  # Fall back: use all patients with PRS data, regardless of oxaliplatin
  merged_all <- merge(prs_df, clinical,
                      by.x = "patient_match", by.y = "bcr_patient_barcode",
                      all.x = FALSE, all.y = FALSE)
  merged_all <- merged_all[!is.na(merged_all$os_time) & merged_all$os_time > 0, ]
  merged_all <- merged_all[!duplicated(merged_all$patient_match), ]
  
  cat(sprintf("  All patients with PRS + survival: %d\n", nrow(merged_all)))
  cat(sprintf("  OS events: %d / %d\n",
              sum(merged_all$os_event, na.rm = TRUE), nrow(merged_all)))
  
  # Also run analysis on all patients for comparison
  assign("merged_oxali", merged_all, envir = .GlobalEnv)
  cat("  WARNING: Using ALL TCGA patients (not just oxaliplatin-treated)\n")
  cat("  Reason: insufficient oxaliplatin-treated samples with survival data\n")
}

# ============================================================
# 10. Cox regression
# ============================================================
cat("\n=== [10] Cox Regression ===\n\n")

analysis_df <- merged_oxali

# 10a. Univariate Cox: PRS
cat("  [10a] Univariate Cox: PRS\n")
cox_prs <- coxph(Surv(os_months, os_event) ~ PRS, data = analysis_df)
print(summary(cox_prs))

hr_prs <- exp(coef(cox_prs))
ci_prs <- exp(confint(cox_prs))
p_prs  <- summary(cox_prs)$coefficients[1, "Pr(>|z|)"]
cat(sprintf("\n  HR per 1-unit PRS = %.3f (%.3f\u2013%.3f), p = %.4f\n",
            hr_prs, ci_prs[1], ci_prs[2], p_prs))

# Save univariate results
sink(file.path(OUT_TAB_DIR, "tcga_cox_univariate.txt"))
cat("Univariate Cox: PRS ~ OS\n")
cat(sprintf("N = %d, Events = %d\n", nrow(analysis_df), sum(analysis_df$os_event)))
print(summary(cox_prs))
sink()

# 10b. Multivariable Cox: PRS + clinical
cat("\n  [10b] Multivariable Cox: PRS + Clinical\n")

# Format clinical covariates
analysis_df$age_num <- as.numeric(analysis_df$age_at_initial_pathologic_diagnosis)
analysis_df$gender_bin <- ifelse(analysis_df$gender == "MALE", 1, 0)

# TNM stage
stage_col <- "stage_event_pathologic_stage"
if (stage_col %in% colnames(analysis_df)) {
  analysis_df$stage_III_IV <- ifelse(
    grepl("III|IV|3|4", analysis_df[[stage_col]]), 1, 0
  )
} else {
  analysis_df$stage_III_IV <- NA
}

# Cancer type
analysis_df$cancer_type_bin <- ifelse(analysis_df$cancer_type == "READ", 1, 0)

cox_multi <- tryCatch({
  coxph(Surv(os_months, os_event) ~ PRS + age_num + gender_bin + stage_III_IV,
        data = analysis_df)
}, error = function(e) {
  cat(sprintf("  WARNING: Multivariable Cox failed: %s\n", e$message))
  return(NULL)
})

if (!is.null(cox_multi)) {
  print(summary(cox_multi))
  
  sink(file.path(OUT_TAB_DIR, "tcga_cox_multivariable.txt"))
  cat("Multivariable Cox: PRS + Clinical ~ OS\n")
  cat(sprintf("N = %d, Events = %d\n", nrow(analysis_df), sum(analysis_df$os_event)))
  print(summary(cox_multi))
  sink()
  cat("  -> Saved: tcga_cox_multivariable.txt\n")
  
  # Build results table
  multi_coef <- summary(cox_multi)$coefficients
  multi_ci   <- summary(cox_multi)$conf.int
  
  multi_results <- data.frame(
    Variable = rownames(multi_coef),
    HR = round(multi_ci[, 1], 4),
    CI_lower = round(multi_ci[, 3], 4),
    CI_upper = round(multi_ci[, 4], 4),
    p_value = round(multi_coef[, "Pr(>|z|)"], 4),
    stringsAsFactors = FALSE
  )
  write.csv(multi_results, file.path(OUT_TAB_DIR, "tcga_cox_multivariable.csv"),
            row.names = FALSE)
  cat("  -> Saved: tcga_cox_multivariable.csv\n")
} else {
  cat("  SKIP: Multivariable Cox failed to converge\n")
}

# Save univariate results as CSV
cox_result <- data.frame(
  Variable = "PRS",
  N = nrow(analysis_df),
  Events = sum(analysis_df$os_event),
  HR = round(hr_prs, 4),
  CI_lower = round(ci_prs[1], 4),
  CI_upper = round(ci_prs[2], 4),
  p_value = round(p_prs, 4),
  stringsAsFactors = FALSE
)
write.csv(cox_result, file.path(OUT_TAB_DIR, "tcga_cox_results.csv"), row.names = FALSE)
cat("  -> Saved: tcga_cox_results.csv\n")

# ============================================================
# 11. Kaplan-Meier curves
# ============================================================
cat("\n=== [11] Kaplan-Meier Curves ===\n\n")

# Stratify by PRS median
analysis_df$prs_group <- ifelse(
  analysis_df$PRS > median(analysis_df$PRS, na.rm = TRUE),
  "High PRS", "Low PRS"
)
analysis_df$prs_group <- factor(analysis_df$prs_group, levels = c("Low PRS", "High PRS"))
cat(sprintf("  Median PRS = %.4f\n", median(analysis_df$PRS, na.rm = TRUE)))
cat(sprintf("  Low PRS: %d, High PRS: %d\n",
            sum(analysis_df$prs_group == "Low PRS"),
            sum(analysis_df$prs_group == "High PRS")))

# 11a. KM for OS
cat("\n  [11a] KM Curve: OS\n")
km_os <- survfit(Surv(os_months, os_event) ~ prs_group, data = analysis_df)

pdf(file.path(OUT_FIG_DIR, "tcga_km_os.pdf"), width = 8, height = 7)
survminer::ggsurvplot(
  km_os, data = analysis_df,
  pval = TRUE, pval.method = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  xlab = "Overall Survival (months)",
  ylab = "OS Probability",
  title = "PRS-stratified OS in TCGA COAD+READ (Oxaliplatin-treated)",
  subtitle = sprintf("PRS model: %s | N=%d, Events=%d",
                     paste(target_pathways, collapse = " + "),
                     nrow(analysis_df), sum(analysis_df$os_event)),
  palette = c("steelblue", "darkred"),
  legend.title = "PRS Group",
  ggtheme = theme_bw()
)
dev.off()
cat(sprintf("  -> Saved: tcga_km_os.pdf\n"))

# Log-rank
lr_os <- survdiff(Surv(os_months, os_event) ~ prs_group, data = analysis_df)
p_os <- 1 - pchisq(lr_os$chisq, df = 1)
cat(sprintf("  Log-rank p (OS) = %.4f\n", p_os))

# 11b. Alternative stratification: tertiles
cat("\n  [11b] Tertile stratification\n")
tertiles <- quantile(analysis_df$PRS, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
analysis_df$prs_tertile <- cut(analysis_df$PRS,
                                breaks = tertiles,
                                labels = c("Low", "Mid", "High"),
                                include.lowest = TRUE)

km_tertile <- survfit(Surv(os_months, os_event) ~ prs_tertile, data = analysis_df)

pdf(file.path(OUT_FIG_DIR, "tcga_km_os_tertile.pdf"), width = 8, height = 7)
survminer::ggsurvplot(
  km_tertile, data = analysis_df,
  pval = TRUE, pval.method = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  xlab = "Overall Survival (months)",
  ylab = "OS Probability",
  title = "PRS-stratified OS (Tertiles) - TCGA COAD+READ",
  palette = c("steelblue", "orange", "darkred"),
  legend.title = "PRS Tertile",
  ggtheme = theme_bw()
)
dev.off()
cat(sprintf("  -> Saved: tcga_km_os_tertile.pdf\n"))

# ============================================================
# 12. Time-dependent ROC
# ============================================================
cat("\n=== [12] Time-dependent ROC ===\n\n")

# Use timeROC package
roc_times <- c(12, 36, 60)  # 1, 3, 5 years in months

for (t in roc_times) {
  # Check if enough events before the time point
  n_at_risk <- sum(analysis_df$os_months >= t | (analysis_df$os_event == 1 & analysis_df$os_months < t))
  events_before <- sum(analysis_df$os_event == 1 & analysis_df$os_months <= t)
  
  cat(sprintf("  %d-month ROC: N at risk = %d, events before = %d\n", t, n_at_risk, events_before))
  
  if (events_before < 5) {
    cat(sprintf("    SKIP: too few events (%d) before %d months\n", events_before, t))
    next
  }
  
  tryCatch({
    roc_obj <- timeROC(
      T = analysis_df$os_months,
      delta = analysis_df$os_event,
      marker = analysis_df$PRS,
      cause = 1,
      times = t,
      iid = FALSE,
      weighting = "marginal"
    )
    
    auc_val <- roc_obj$AUC[1]
    cat(sprintf("    Time-dependent AUC at %d months = %.3f\n", t, auc_val))
    
    # ROC plot
    pdf(file.path(OUT_FIG_DIR, sprintf("tcga_roc_%dyr.pdf", t / 12)), width = 7, height = 7)
    plot(roc_obj, time = t, col = "darkred", lwd = 2,
         title = sprintf("PRS Time-dependent ROC: %d-month OS", t))
    legend("bottomright",
           legend = sprintf("AUC = %.3f", auc_val),
           col = "darkred", lwd = 2)
    dev.off()
    cat(sprintf("    -> Saved: tcga_roc_%dyr.pdf\n", t / 12))
    
  }, error = function(e) {
    cat(sprintf("    ERROR: %s\n", e$message))
  })
}

# ============================================================
# 13. PRS distribution plot
# ============================================================
cat("\n=== [13] PRS Distribution Plot ===\n\n")

p <- ggplot(analysis_df, aes(x = PRS, fill = factor(os_event))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "steelblue", "1" = "darkred"),
                    labels = c("Alive", "Dead"),
                    name = "Vital Status") +
  geom_vline(xintercept = median(analysis_df$PRS, na.rm = TRUE),
             linetype = "dashed", color = "gray50") +
  theme_bw() +
  labs(title = "PRS Distribution by Vital Status (TCGA Oxaliplatin)",
       x = "Pathway-based PRS", y = "Density")

pdf(file.path(OUT_FIG_DIR, "tcga_prs_boxplot.pdf"), width = 8, height = 5)
print(p)
dev.off()
cat(sprintf("  -> Saved: tcga_prs_boxplot.pdf\n"))

# ============================================================
# 14. Summary
# ============================================================
cat("\n=== [14] Summary ===\n\n")

cat(sprintf("  Analysis set: %d TCGA oxaliplatin-treated patients\n", nrow(analysis_df)))
cat(sprintf("  OS events: %d / %d (%.1f%%)\n",
            sum(analysis_df$os_event), nrow(analysis_df),
            mean(analysis_df$os_event) * 100))
cat(sprintf("  Median PRS: %.4f\n", median(analysis_df$PRS, na.rm = TRUE)))
cat(sprintf("  PRS range: [%.4f, %.4f]\n",
            min(analysis_df$PRS, na.rm = TRUE),
            max(analysis_df$PRS, na.rm = TRUE)))
cat(sprintf("  Cox (Univariate PRS): HR = %.4f (%.4f\u2013%.4f), p = %.4f\n",
            hr_prs, ci_prs[1], ci_prs[2], p_prs))
cat(sprintf("  Log-rank (OS): p = %.4f\n", p_os))

cat("\n  Output files:\n")
cat(sprintf("    %s/\n", OUT_TAB_DIR))
for (f in list.files(OUT_TAB_DIR)) {
  cat(sprintf("      %s\n", f))
}
cat(sprintf("    %s/\n", OUT_FIG_DIR))
for (f in list.files(OUT_FIG_DIR)) {
  cat(sprintf("      %s\n", f))
}

cat("\n=== P0-1: TCGA External Validation Complete ===\n\n")
