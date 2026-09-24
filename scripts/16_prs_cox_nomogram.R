# ============================================================
# 16_PRS_COX_NOMOGRAM.R
# Phase III: PRS Construction, Cox Regression & Nomogram
#
# Purpose:
#   Build a pathway-based Polygenic Risk Score (PRS) for XELOX
#   resistance prediction using a merged training cohort
#   (GSE39582 XELOX-like + GSE104645 oxaliplatin).
#   Two-step feature reduction:
#     Step 1: Univariate AUC pre-filter (AUC > 0.55 or < 0.45)
#     Step 2: LASSO-Logistic regression for final selection
#   Then:
#     - PRS = coefficient-weighted linear combination
#     - AUC validation in training + 3 external cohorts
#     - Cox regression (PRS + clinical) on GSE39582 RFS data
#     - Nomogram visualization
#
# Training cohort (Path C):
#   GSE39582 XELOX-like: 119S + 45R = 164 (RFS proxy)
#   GSE104645 oxaliplatin: 66S + 47R = 113 (true RECIST response)
#   Total: 185S + 92R = 277
#
# Validation cohorts:
#   GSE28702: 83 pure FOLFOX (42R + 41NR)
#   GSE69657: 30 pure XELOX (13R + 17NR)
#
# Input:
#   results/tables/pathway_activity/*_pathway_scores.rds  (44 pathways)
#   results/tables/pathway_activity/*_gene_expression.rds (for genes)
#   results/tables/GSE39582_xelox_groups.csv
#   results/tables/GSE*_clinical_data.csv
#
# Output:
#   tables/prs/PRS_model_coefficients.csv     — LASSO-selected pathways + weights
#   tables/prs/PRS_scores_training.csv        — PRS for training samples
#   tables/prs/PRS_validation_summary.csv     — AUC for all cohorts
#   tables/prs/Cox_PRS_results.csv            — Cox regression HR + CI
#   tables/prs/Cox_multivariable_results.csv  — Multivariable Cox HR + CI
#   figures/cox_prs*_forest.pdf               — Forest plot
#   figures/prs*_roc.pdf                      — ROC curves
#   figures/prs*_km.pdf                       — Kaplan-Meier curves
#   figures/nomogram_prs.pdf                  — Nomogram
#   figures/prs*_calibration.pdf              — Calibration plot
# ============================================================

Sys.setenv(TMPDIR = "/tmp", TMP = "/tmp", TEMP = "/tmp")
.libPaths(c("/path/to/Rlibs", .libPaths()))

# PROJECT_ROOT from environment variable
PROJECT_ROOT <- Sys.getenv("XELOX_ROOT",
  "/path/to/xelox_project基于可解释性机器学习的XELOX耐药分子指纹研究")

RESULTS_TAB_DIR  <- file.path(PROJECT_ROOT, "results", "tables")
RESULTS_FIG_DIR  <- file.path(PROJECT_ROOT, "results", "figures")
PA_DIR           <- file.path(RESULTS_TAB_DIR, "pathway_activity")
PRS_TAB_DIR      <- file.path(RESULTS_TAB_DIR, "prs")
PRS_FIG_DIR      <- file.path(RESULTS_FIG_DIR, "prs")

suppressWarnings(dir.create(PRS_TAB_DIR, showWarnings = FALSE, recursive = TRUE))
suppressWarnings(dir.create(PRS_FIG_DIR, showWarnings = FALSE, recursive = TRUE))

set.seed(42)

cat("============================================================\n")
cat("Phase III: PRS Construction, Cox Regression & Nomogram\n")
cat("Path C — Merged Training (GSE39582 + GSE104645)\n")
cat("============================================================\n\n")

# ============================================================
# 0. Load required packages
# ============================================================
cat("=== [0] Loading packages ===\n\n")
required_pkgs <- c("glmnet", "pROC", "survival", "survminer",
                    "rms", "ggplot2", "data.table")
for (pkg in required_pkgs) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  cat(sprintf("  %s v%s loaded\n", pkg, packageVersion(pkg)))
}
cat("\n")

# ============================================================
# 1. Load pathway scores for all datasets
# ============================================================
cat("=== [1] Loading pathway scores ===\n\n")

load_pathway_data <- function(prefix) {
  score_file <- file.path(PA_DIR, paste0(prefix, "_pathway_scores.rds"))
  ps <- readRDS(score_file)
  cat(sprintf("  %s pathway scores: %d pathways x %d samples\n",
              prefix, nrow(ps), ncol(ps)))
  return(ps)
}

ps_39582   <- load_pathway_data("GSE39582")
ps_104645  <- load_pathway_data("GSE104645")
ps_28702   <- load_pathway_data("GSE28702")
ps_69657   <- load_pathway_data("GSE69657")
ps_72970   <- load_pathway_data("GSE72970")

# Verify pathway alignment
pathways <- rownames(ps_39582)
stopifnot(all(rownames(ps_104645) == pathways))
stopifnot(all(rownames(ps_28702) == pathways))
stopifnot(all(rownames(ps_69657) == pathways))
stopifnot(all(rownames(ps_72970) == pathways))
cat(sprintf("\n  All datasets aligned: %d common pathways\n\n", length(pathways)))

# ============================================================
# 2. Annotate samples with response labels
# ============================================================
cat("=== [2] Annotating sample groups ===\n\n")

# --- 2a. GSE39582 (RFS proxy) ---
groups_39582 <- read.csv(file.path(RESULTS_TAB_DIR, "GSE39582_xelox_groups.csv"),
                          stringsAsFactors = FALSE)
rownames(groups_39582) <- groups_39582$sample_id

# Map to binary response (exclude intermediate)
g39582_binary <- groups_39582$group
names(g39582_binary) <- groups_39582$sample_id
g39582_keep <- g39582_binary %in% c("sensitive", "resistant")
g39582_y <- ifelse(g39582_binary[g39582_keep] == "resistant", 1, 0)

# Match to pathway scores
g39582_common <- intersect(names(g39582_y), colnames(ps_39582))
g39582_y <- g39582_y[g39582_common]
cat(sprintf("  GSE39582 XELOX: Sensitive=%d, Resistant=%d (RFS proxy)\n",
            sum(g39582_y == 0), sum(g39582_y == 1)))

# --- 2b. GSE104645 (RECIST response) ---
clin_104645 <- read.csv(file.path(RESULTS_TAB_DIR, "GSE104645_clinical_data.csv"),
                         row.names = 1, stringsAsFactors = FALSE)

regimen_col_104645 <- "X1st.line.chemotherapy.regimens.ch1"
response_col_104645 <- "best.response.of.1st.line.chemotherapy.ch1"

is_oxali <- grepl("XELOX|FOLFOX|SOX|CAPOX",
                   clin_104645[[regimen_col_104645]], ignore.case = TRUE)
oxali_samples <- rownames(clin_104645)[is_oxali]

resp <- clin_104645[[response_col_104645]]
names(resp) <- rownames(clin_104645)

g104645_y <- rep(NA, length(oxali_samples))
names(g104645_y) <- oxali_samples
is_sens <- resp[oxali_samples] %in% c("Complete response", "Partial response")
is_res  <- resp[oxali_samples] %in% c("Progressive disease", "Stable disease")
g104645_y[is_sens] <- 0  # Sensitive
g104645_y[is_res]  <- 1  # Resistant

# Remove NA (Not evaluated)
g104645_y <- g104645_y[!is.na(g104645_y)]

# Match to pathway scores
g104645_common <- intersect(names(g104645_y), colnames(ps_104645))
g104645_y <- g104645_y[g104645_common]
cat(sprintf("  GSE104645 oxali: Sensitive=%d, Resistant=%d (RECIST)\n",
            sum(g104645_y == 0), sum(g104645_y == 1)))

# --- 2c. Validation: GSE28702 (mFOLFOX6) ---
clin_28702 <- read.csv(file.path(RESULTS_TAB_DIR, "GSE28702_clinical_data.csv"),
                        row.names = 1, stringsAsFactors = FALSE)
resp_28702 <- clin_28702[["mfolfox6.ch1"]]
names(resp_28702) <- rownames(clin_28702)
g28702_y <- ifelse(resp_28702 == "responder", 0,
                    ifelse(resp_28702 == "non-responder", 1, NA))
g28702_common <- intersect(names(g28702_y[!is.na(g28702_y)]), colnames(ps_28702))
g28702_y <- g28702_y[g28702_common]
cat(sprintf("  GSE28702 FOLFOX: Sensitive=%d, Resistant=%d\n",
            sum(g28702_y == 0), sum(g28702_y == 1)))

# --- 2d. Validation: GSE69657 (FOLFOX4) ---
clin_69657 <- read.csv(file.path(RESULTS_TAB_DIR, "GSE69657_clinical_data.csv"),
                        row.names = 1, stringsAsFactors = FALSE)
resp_69657 <- clin_69657[["chemoresponse.ch1"]]
names(resp_69657) <- rownames(clin_69657)
g69657_y <- ifelse(resp_69657 == "responder", 0,
                    ifelse(resp_69657 == "noresponder", 1, NA))
g69657_common <- intersect(names(g69657_y[!is.na(g69657_y)]), colnames(ps_69657))
g69657_y <- g69657_y[g69657_common]
cat(sprintf("  GSE69657 XELOX: Sensitive=%d, Resistant=%d\n",
            sum(g69657_y == 0), sum(g69657_y == 1)))

cat("\n")

# ============================================================
# 3. Build merged training matrix
# ============================================================
cat("=== [3] Building merged training matrix ===\n\n")

# Training: GSE39582 + GSE104645
train_ps <- cbind(
  ps_39582[, names(g39582_y), drop = FALSE],
  ps_104645[, names(g104645_y), drop = FALSE]
)
train_y <- c(g39582_y, g104645_y)
train_dataset <- c(rep("GSE39582", length(g39582_y)),
                    rep("GSE104645", length(g104645_y)))

cat(sprintf("  Training matrix: %d pathways x %d samples\n",
            nrow(train_ps), ncol(train_ps)))
cat(sprintf("  Outcome: Resistant=%d, Sensitive=%d\n",
            sum(train_y == 1), sum(train_y == 0)))
cat(sprintf("  Resistance rate: %.1f%%\n\n", mean(train_y) * 100))

# Transpose: samples x pathways
train_x <- t(train_ps)
cat(sprintf("  Training X dim: %d samples x %d features\n",
            nrow(train_x), ncol(train_x)))

# Store pathway names for reference
pathway_names <- colnames(train_x)

# ============================================================
# 4. Step 1: Univariate AUC pre-filter
# ============================================================
cat("=== [4] Step 1: Univariate AUC pre-filter ===\n\n")

compute_auc <- function(x, y) {
  tryCatch({
    roc_obj <- roc(y, x, direction = "<", quiet = TRUE)
    as.numeric(auc(roc_obj))
  }, error = function(e) NA)
}

# Compute AUC for each pathway
auc_vec <- apply(train_ps, 1, function(pw) {
  compute_auc(pw[train_dataset == "GSE39582"],
              train_y[train_dataset == "GSE39582"])
})

names(auc_vec) <- pathway_names

# Use GSE39582 AUC for filtering (pure RFS proxy, no data leakage)
# Filter: AUC > 0.55 (pro-resistance) or AUC < 0.45 (pro-sensitive)
auc_filter <- abs(auc_vec - 0.5) > 0.05
auc_filter[is.na(auc_filter)] <- FALSE

cat(sprintf("  Passed univariate filter (|AUC-0.5| > 0.05): %d / %d pathways\n",
            sum(auc_filter), length(auc_filter)))

# Print top pathways
auc_sorted <- sort(auc_vec[auc_filter], decreasing = TRUE)
cat("\n  Top 10 pro-resistance pathways (AUC > 0.55):\n")
top_res <- head(auc_sorted[auc_sorted > 0.5], 10)
for (i in seq_along(top_res)) {
  cat(sprintf("    %2d. %s  (AUC = %.3f)\n", i, names(top_res)[i], top_res[i]))
}
cat("\n  Top 10 pro-sensitive pathways (AUC < 0.45):\n")
top_sens <- head(auc_sorted[auc_sorted < 0.5], 10)
for (i in seq_along(top_sens)) {
  cat(sprintf("    %2d. %s  (AUC = %.3f)\n", i, names(top_sens)[i], top_sens[i]))
}

# Filter training data
train_x_filtered <- train_x[, auc_filter, drop = FALSE]
cat(sprintf("\n  Filtered training X: %d samples x %d features\n",
            nrow(train_x_filtered), ncol(train_x_filtered)))

# ============================================================
# 5. Step 2: LASSO-Logistic feature selection
# ============================================================
cat("\n=== [5] Step 2: LASSO-Logistic feature selection ===\n\n")

# Scale features for LASSO
train_x_scaled <- scale(train_x_filtered)

# Cross-validated LASSO (binomial)
set.seed(42)
cv_lasso <- cv.glmnet(
  x = train_x_scaled,
  y = train_y,
  family = "binomial",
  alpha = 1,
  nfolds = 10,
  type.measure = "deviance"
)

# Extract coefficients at lambda.1se
coef_lasso <- as.matrix(coef(cv_lasso, s = "lambda.1se"))
selected_idx <- which(coef_lasso[-1, 1] != 0)  # exclude intercept
selected_pathways <- pathway_names[auc_filter][selected_idx]

cat(sprintf("  Number of selected pathways (lambda.1se): %d\n", length(selected_pathways)))

if (length(selected_pathways) == 0) {
  cat("\n  WARNING: No pathways selected at lambda.1se!\n")
  cat("  Falling back to lambda.min...\n")
  coef_lasso <- as.matrix(coef(cv_lasso, s = "lambda.min"))
  selected_idx <- which(coef_lasso[-1, 1] != 0)
  selected_pathways <- pathway_names[auc_filter][selected_idx]
  cat(sprintf("  Number of selected pathways (lambda.min): %d\n", length(selected_pathways)))
}

# Build coefficient table
coef_table <- data.frame(
  Pathway = selected_pathways,
  Coefficient = coef_lasso[selected_idx + 1, 1],
  stringsAsFactors = FALSE
)
# Add direction
coef_table$Direction <- ifelse(coef_table$Coefficient > 0, "Pro-Resistance", "Pro-Sensitive")
# Sort by absolute coefficient
coef_table <- coef_table[order(abs(coef_table$Coefficient), decreasing = TRUE), ]

cat("\n  Selected pathways:\n")
for (i in seq_len(nrow(coef_table))) {
  cat(sprintf("    %2d. %s  (coef = %+.4f, %s)\n",
              i, coef_table$Pathway[i],
              coef_table$Coefficient[i],
              coef_table$Direction[i]))
}

# Save coefficient table
write.csv(coef_table, file.path(PRS_TAB_DIR, "PRS_model_coefficients.csv"),
          row.names = FALSE)
cat(sprintf("\n  -> Saved: %s\n",
            file.path(PRS_TAB_DIR, "PRS_model_coefficients.csv")))

# ============================================================
# 6. PRS construction
# ============================================================
cat("\n=== [6] PRS Construction ===\n\n")

build_prs <- function(expr_matrix, coef_vec) {
  # expr_matrix: samples x pathways (scaled to training mean/sd)
  # coef_vec: named vector of coefficients
  common_pw <- intersect(names(coef_vec), colnames(expr_matrix))
  if (length(common_pw) == 0) return(rep(NA, nrow(expr_matrix)))
  score <- as.numeric(expr_matrix[, common_pw, drop = FALSE] %*% coef_vec[common_pw])
  return(score)
}

# Extract coefficient vector (named)
coef_vec <- coef_table$Coefficient
names(coef_vec) <- coef_table$Pathway

# --- 6a. Scale all datasets using training mean/sd ---
train_mean <- attr(train_x_scaled, "scaled:center")
train_sd   <- attr(train_x_scaled, "scaled:scale")

scale_with_train <- function(x_mat) {
  scale(x_mat, center = train_mean[colnames(x_mat)],
        scale = train_sd[colnames(x_mat)])
}

# Build PRS for each dataset
prs_train <- build_prs(train_x_scaled, coef_vec)

# GSE39582 (all samples, including intermediate for survival)
prs_39582 <- build_prs(
  scale_with_train(t(ps_39582)),
  coef_vec
)
names(prs_39582) <- colnames(ps_39582)

# Validation datasets
prs_104645_val <- build_prs(
  scale_with_train(t(ps_104645)),
  coef_vec
)
names(prs_104645_val) <- colnames(ps_104645)

prs_28702 <- build_prs(
  scale_with_train(t(ps_28702)),
  coef_vec
)
names(prs_28702) <- colnames(ps_28702)

prs_69657 <- build_prs(
  scale_with_train(t(ps_69657)),
  coef_vec
)
names(prs_69657) <- colnames(ps_69657)

prs_72970 <- build_prs(
  scale_with_train(t(ps_72970)),
  coef_vec
)
names(prs_72970) <- colnames(ps_72970)

cat("  PRS computed for all datasets\n")
cat(sprintf("  GSE39582: %d samples\n", length(prs_39582)))
cat(sprintf("  GSE104645: %d samples\n", length(prs_104645_val)))
cat(sprintf("  GSE28702: %d samples\n", length(prs_28702)))
cat(sprintf("  GSE69657: %d samples\n", length(prs_69657)))
cat(sprintf("  GSE72970: %d samples\n", length(prs_72970)))

# Save training PRS
prs_train_df <- data.frame(
  sample_id = names(train_y),
  PRS = prs_train,
  response = ifelse(train_y == 1, "Resistant", "Sensitive"),
  dataset = train_dataset,
  stringsAsFactors = FALSE
)
write.csv(prs_train_df, file.path(PRS_TAB_DIR, "PRS_scores_training.csv"),
          row.names = FALSE)
cat(sprintf("  -> Saved: %s\n",
            file.path(PRS_TAB_DIR, "PRS_scores_training.csv")))

# ============================================================
# 7. AUC evaluation in training + validation cohorts
# ============================================================
cat("\n=== [7] AUC Evaluation ===\n\n")

evaluate_auc <- function(prs_vec, y_vec, label) {
  if (length(unique(y_vec)) < 2 || sd(prs_vec, na.rm = TRUE) == 0) {
    cat(sprintf("  %-15s — Cannot compute (insufficient variation)\n", label))
    return(data.frame(Cohort = label, N = length(y_vec),
                       AUC = NA, CI_lower = NA, CI_upper = NA))
  }
  roc_obj <- roc(y_vec, prs_vec, direction = "<", quiet = TRUE)
  auc_val <- auc(roc_obj)
  ci_val <- ci.auc(roc_obj)
  cat(sprintf("  %-15s AUC = %.3f (%.3f–%.3f), N=%d\n",
              label, auc_val, ci_val[1], ci_val[3], length(y_vec)))
  data.frame(Cohort = label, N = length(y_vec),
             AUC = round(auc_val, 3),
             CI_lower = round(ci_val[1], 3),
             CI_upper = round(ci_val[3], 3))
}

auc_results <- rbind(
  evaluate_auc(prs_train, train_y, "Training (merged)"),
  evaluate_auc(prs_39582[names(g39582_y)], g39582_y, "GSE39582"),
  evaluate_auc(prs_104645_val[names(g104645_y)], g104645_y, "GSE104645"),
  evaluate_auc(prs_28702[names(g28702_y)], g28702_y, "GSE28702"),
  evaluate_auc(prs_69657[names(g69657_y)], g69657_y, "GSE69657")
)

write.csv(auc_results, file.path(PRS_TAB_DIR, "PRS_validation_summary.csv"),
          row.names = FALSE)
cat(sprintf("\n  -> Saved: %s\n",
            file.path(PRS_TAB_DIR, "PRS_validation_summary.csv")))

# --- 7b. ROC plots ---
cat("\n  Generating ROC plots...\n")

plot_roc <- function(prs_vec, y_vec, label, color) {
  roc_obj <- roc(y_vec, prs_vec, direction = "<", quiet = TRUE)
  lines(roc_obj, col = color, lwd = 2)
  text(0.8, seq(0.1, 0.5, length.out = nrow(auc_results))[i], 
       paste0(label, " AUC = ", round(auc(roc_obj), 3)), col = color, cex = 0.9)
}

# Combined ROC for training cohorts
pdf(file.path(PRS_FIG_DIR, "prs_roc_combined.pdf"), width = 8, height = 7)
plot(roc(train_y, prs_train, direction = "<", quiet = TRUE),
     col = "darkred", lwd = 2.5, main = "PRS Performance: Training Cohorts")
lines(roc(g39582_y, prs_39582[names(g39582_y)], direction = "<", quiet = TRUE),
      col = "orange", lwd = 2)
lines(roc(g104645_y, prs_104645_val[names(g104645_y)], direction = "<", quiet = TRUE),
      col = "steelblue", lwd = 2)
legend("bottomright",
       legend = c(sprintf("Merged (AUC=%.3f)", auc(roc(train_y, prs_train, direction="<", quiet=TRUE))),
                  sprintf("GSE39582 (AUC=%.3f)", auc(roc(g39582_y, prs_39582[names(g39582_y)], direction="<", quiet=TRUE))),
                  sprintf("GSE104645 (AUC=%.3f)", auc(roc(g104645_y, prs_104645_val[names(g104645_y)], direction="<", quiet=TRUE)))),
       col = c("darkred", "orange", "steelblue"), lwd = 2, cex = 0.9)
dev.off()
cat(sprintf("  -> Saved: %s\n", file.path(PRS_FIG_DIR, "prs_roc_combined.pdf")))

# Combined ROC for validation cohorts
if (length(unique(g28702_y)) == 2 && sd(prs_28702[names(g28702_y)]) > 0 &&
    length(unique(g69657_y)) == 2 && sd(prs_69657[names(g69657_y)]) > 0) {
  pdf(file.path(PRS_FIG_DIR, "prs_roc_validation.pdf"), width = 8, height = 7)
  plot(roc(g28702_y, prs_28702[names(g28702_y)], direction = "<", quiet = TRUE),
       col = "forestgreen", lwd = 2.5, main = "PRS Performance: Validation Cohorts")
  lines(roc(g69657_y, prs_69657[names(g69657_y)], direction = "<", quiet = TRUE),
        col = "purple", lwd = 2)
  legend("bottomright",
         legend = c(sprintf("GSE28702 FOLFOX (AUC=%.3f)", auc(roc(g28702_y, prs_28702[names(g28702_y)], direction="<", quiet=TRUE))),
                    sprintf("GSE69657 XELOX (AUC=%.3f)", auc(roc(g69657_y, prs_69657[names(g69657_y)], direction="<", quiet=TRUE)))),
         col = c("forestgreen", "purple"), lwd = 2, cex = 0.9)
  dev.off()
  cat(sprintf("  -> Saved: %s\n", file.path(PRS_FIG_DIR, "prs_roc_validation.pdf")))
}

# ============================================================
# 8. Cox regression on GSE39582 (RFS)
# ============================================================
cat("\n=== [8] Cox Regression (GSE39582 RFS) ===\n\n")

# Build GSE39582 survival data
cox_df <- groups_39582[, c("sample_id", "rfs_event", "rfs_delay",
                             "sex", "age", "tnm_stage", "tumor_location",
                             "mmr_status", "kras", "braf", "cit_subtype", "group")]
cox_df$PRS <- prs_39582[cox_df$sample_id]

# Remove missing survival data
cox_df <- cox_df[!is.na(cox_df$rfs_event) & !is.na(cox_df$rfs_delay) &
                   cox_df$rfs_delay > 0, ]
cox_df$rfs_event <- as.numeric(cox_df$rfs_event)

cat(sprintf("  GSE39582 samples with valid RFS: %d\n", nrow(cox_df)))
cat(sprintf("  RFS events: %d / %d (%.1f%%)\n",
            sum(cox_df$rfs_event), nrow(cox_df),
            mean(cox_df$rfs_event) * 100))

# --- 8a. Univariate Cox: PRS only ---
cat("\n  [8a] Univariate Cox: PRS\n")
cox_prs <- coxph(Surv(rfs_delay, rfs_event) ~ PRS, data = cox_df)
sink(file.path(PRS_TAB_DIR, "Cox_univariate_PRS.txt"))
print(summary(cox_prs))
sink()
cat("  -> Saved Cox summary\n")

# Extract HR
hr_prs <- exp(coef(cox_prs))
ci_prs <- exp(confint(cox_prs))
cat(sprintf("  HR per 1-unit PRS = %.3f (%.3f–%.3f), p = %.4f\n",
            hr_prs, ci_prs[1], ci_prs[2], summary(cox_prs)$coefficients[1, "Pr(>|z|)"]))

# --- 8b. Multivariable Cox: PRS + clinical ---
cat("\n  [8b] Multivariable Cox: PRS + Clinical\n")

# Format clinical variables
cox_df$sex_bin <- ifelse(cox_df$sex == "Female", 0, 1)
cox_df$age_num <- as.numeric(cox_df$age)
cox_df$stage_III_IV <- ifelse(grepl("3|4|III|IV", cox_df$tnm_stage), 1, 0)
cox_df$mmr_dmmr <- ifelse(cox_df$mmr_status == "dMMR", 1, 0)
cox_df$location_proximal <- ifelse(cox_df$tumor_location == "proximal", 1, 0)

cox_multi <- coxph(
  Surv(rfs_delay, rfs_event) ~ PRS + age_num + sex_bin + stage_III_IV + 
    location_proximal + mmr_dmmr,
  data = cox_df
)

sink(file.path(PRS_TAB_DIR, "Cox_multivariable.txt"))
print(summary(cox_multi))
sink()
cat("  -> Saved multivariable Cox summary\n")

# Build results table
multi_coef <- summary(cox_multi)$coefficients
multi_ci <- summary(cox_multi)$conf.int
multi_results <- data.frame(
  Variable = rownames(multi_coef),
  HR = round(multi_ci[, 1], 3),
  CI_lower = round(multi_ci[, 3], 3),
  CI_upper = round(multi_ci[, 4], 3),
  p_value = round(multi_coef[, "Pr(>|z|)"], 4),
  stringsAsFactors = FALSE
)
write.csv(multi_results, file.path(PRS_TAB_DIR, "Cox_multivariable_results.csv"),
          row.names = FALSE)
cat(sprintf("  -> Saved: %s\n",
            file.path(PRS_TAB_DIR, "Cox_multivariable_results.csv")))

# --- 8c. Forest plot ---
pdf(file.path(PRS_FIG_DIR, "cox_prs_forest.pdf"), width = 8, height = 5)
par(mar = c(4, 8, 3, 4))
plot(NA, xlim = c(0, max(multi_ci[, 3:4], na.rm = TRUE) * 1.2),
     ylim = c(0.5, nrow(multi_results) + 1),
     xlab = "Hazard Ratio (95% CI)", ylab = "", yaxt = "n",
     main = "Multivariable Cox Regression: PRS + Clinical")
axis(2, at = nrow(multi_results):1, labels = rev(multi_results$Variable), las = 2)
abline(v = 1, lty = 2, col = "gray")
for (i in seq_len(nrow(multi_results))) {
  idx <- nrow(multi_results) - i + 1
  points(multi_ci[idx, 1], i, pch = 15, cex = 1.2)
  segments(multi_ci[idx, 3], i, multi_ci[idx, 4], i, lwd = 2)
  text(max(multi_ci[, 3:4], na.rm = TRUE) * 1.05, i,
       sprintf("%.2f (%.2f–%.2f)", multi_ci[idx, 1], multi_ci[idx, 3], multi_ci[idx, 4]),
       cex = 0.8, adj = 0)
}
dev.off()
cat(sprintf("  -> Saved: %s\n", file.path(PRS_FIG_DIR, "cox_prs_forest.pdf")))

# ============================================================
# 9. Kaplan-Meier curves
# ============================================================
cat("\n=== [9] Kaplan-Meier Curves ===\n\n")

# Stratify PRS by median
cox_df$prs_group <- ifelse(cox_df$PRS > median(cox_df$PRS, na.rm = TRUE),
                            "High PRS", "Low PRS")
cox_df$prs_group <- factor(cox_df$prs_group, levels = c("Low PRS", "High PRS"))

# KM fit
km_fit <- survfit(Surv(rfs_delay, rfs_event) ~ prs_group, data = cox_df)

pdf(file.path(PRS_FIG_DIR, "prs_km_curve.pdf"), width = 8, height = 7)
survminer::ggsurvplot(
  km_fit,
  data = cox_df,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  xlab = "Relapse-Free Survival (months)",
  ylab = "RFS Probability",
  title = "PRS-stratified RFS in GSE39582 XELOX-like Cohort",
  palette = c("steelblue", "darkred"),
  legend.title = "PRS Group",
  ggtheme = theme_bw()
)
dev.off()
cat(sprintf("  -> Saved: %s\n", file.path(PRS_FIG_DIR, "prs_km_curve.pdf")))

# Log-rank test
lr_test <- survdiff(Surv(rfs_delay, rfs_event) ~ prs_group, data = cox_df)
lr_p <- 1 - pchisq(lr_test$chisq, df = 1)
cat(sprintf("  Log-rank test p = %.4f\n", lr_p))

# ============================================================
# 10. Nomogram (PRS-based)
# ============================================================
cat("\n=== [10] Nomogram ===\n\n")

# Use rms package functions
dd <- datadist(cox_df[, c("PRS", "age_num", "sex_bin", "stage_III_IV",
                           "location_proximal", "mmr_dmmr", "rfs_delay", "rfs_event")])
options(datadist = "dd")

# Fit Cox model using rms::cph
cox_cph <- cph(Surv(rfs_delay, rfs_event) ~ PRS + age_num + sex_bin + 
                 stage_III_IV + location_proximal + mmr_dmmr,
               data = cox_df, x = TRUE, y = TRUE, surv = TRUE)

# Define survival time points (12-month, 36-month, 60-month)
surv_times <- c(12, 36, 60)

# Generate nomogram using rms::Survival function
surv_fn_12 <- function(lp) Survival(cox_cph)(times = surv_times[1], lp = lp)
surv_fn_36 <- function(lp) Survival(cox_cph)(times = surv_times[2], lp = lp)
surv_fn_60 <- function(lp) Survival(cox_cph)(times = surv_times[3], lp = lp)

pdf(file.path(PRS_FIG_DIR, "nomogram_prs.pdf"), width = 12, height = 8)
nom <- nomogram(cox_cph, fun = list(surv_fn_12, surv_fn_36, surv_fn_60),
  funlabel = c(
    paste0(surv_times[1], "-month RFS"),
    paste0(surv_times[2], "-month RFS"),
    paste0(surv_times[3], "-month RFS")),
  maxscale = 100,
  lp = TRUE,
  lp.at = seq(-2, 2, by = 0.5))
plot(nom, xfrac = 0.35, cex.var = 0.8, cex.axis = 0.7)
dev.off()
cat(sprintf("  -> Saved: %s\n", file.path(PRS_FIG_DIR, "nomogram_prs.pdf")))

# --- 10b. Calibration plot ---
cat("\n  [10b] Calibration Plot\n")

# Bootstrap calibration for 12-month and 36-month
for (t_idx in seq_along(surv_times)) {
  st <- surv_times[t_idx]
  cal_file <- file.path(PRS_FIG_DIR, sprintf("calibration_%dmo_prs.pdf", st))
  
  tryCatch({
    pdf(cal_file, width = 7, height = 7)
    cal <- calibrate(cox_cph, u = st, B = 200)
    plot(cal, subtitles = FALSE,
         xlab = paste0("Predicted ", st, "-month RFS"),
         ylab = paste0("Observed ", st, "-month RFS"),
         main = paste0("Calibration (", st, " months)"))
    dev.off()
    cat(sprintf("  -> Saved: %s\n", cal_file))
  }, error = function(e) {
    cat(sprintf("  WARNING: Calibration for %dmo failed: %s\n", st, e$message))
    dev.off()
  })
}

# ============================================================
# 11. PRS distribution comparison
# ============================================================
cat("\n=== [11] PRS Distribution Comparison ===\n\n")

# Compare PRS between Sensitive and Resistant in each dataset
compare_prs <- function(prs_vec, y_vec, label) {
  if (length(unique(y_vec)) < 2 || sum(!is.na(prs_vec)) < 3) {
    cat(sprintf("  %-15s — Cannot compute\n", label))
    return(NULL)
  }
  t <- t.test(prs_vec[y_vec == 0], prs_vec[y_vec == 1])
  cat(sprintf("  %-15s Sens=%.3f±%.3f, Res=%.3f±%.3f, p=%.4f\n",
              label,
              mean(prs_vec[y_vec == 0], na.rm = TRUE),
              sd(prs_vec[y_vec == 0], na.rm = TRUE),
              mean(prs_vec[y_vec == 1], na.rm = TRUE),
              sd(prs_vec[y_vec == 1], na.rm = TRUE),
              t$p.value))
}

compare_prs(prs_train, train_y, "Training")
compare_prs(prs_39582[names(g39582_y)], g39582_y, "GSE39582")
compare_prs(prs_104645_val[names(g104645_y)], g104645_y, "GSE104645")
compare_prs(prs_28702[names(g28702_y)], g28702_y, "GSE28702")
compare_prs(prs_69657[names(g69657_y)], g69657_y, "GSE69657")

# --- 11b. Boxplot ---
pdf(file.path(PRS_FIG_DIR, "prs_boxplot_all.pdf"), width = 10, height = 6)

# Build combined data frame for boxplot
all_prs_list <- list(
  Training = prs_train,
  GSE39582 = prs_39582[names(g39582_y)],
  GSE104645 = prs_104645_val[names(g104645_y)],
  GSE28702 = prs_28702[names(g28702_y)],
  GSE69657 = prs_69657[names(g69657_y)]
)
all_y_list <- list(
  Training = train_y,
  GSE39582 = g39582_y,
  GSE104645 = g104645_y,
  GSE28702 = g28702_y,
  GSE69657 = g69657_y
)

boxplot_df <- do.call(rbind, lapply(names(all_prs_list), function(nm) {
  data.frame(
    Dataset = nm,
    PRS = all_prs_list[[nm]],
    Response = ifelse(all_y_list[[nm]] == 1, "Resistant", "Sensitive"),
    stringsAsFactors = FALSE
  )
}))

boxplot_df$Dataset <- factor(boxplot_df$Dataset,
                              levels = names(all_prs_list))
boxplot_df$Response <- factor(boxplot_df$Response,
                               levels = c("Sensitive", "Resistant"))

p <- ggplot(boxplot_df, aes(x = Dataset, y = PRS, fill = Response)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15),
              size = 0.8, alpha = 0.5) +
  scale_fill_manual(values = c("Sensitive" = "steelblue", "Resistant" = "darkred")) +
  theme_bw() +
  labs(title = "PRS Distribution by Response Status Across Cohorts",
       y = "Pathway-based PRS", x = "") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(p)
dev.off()
cat(sprintf("  -> Saved: %s\n", file.path(PRS_FIG_DIR, "prs_boxplot_all.pdf")))

# ============================================================
# 12. Summary output
# ============================================================
cat("\n=== [12] Analysis Summary ===\n\n")

cat(sprintf("Training merged cohort: %d samples (%d Sensitive, %d Resistant)\n",
            length(train_y), sum(train_y == 0), sum(train_y == 1)))
cat(sprintf("Pathways passed univariate filter: %d / %d\n",
            sum(auc_filter), length(auc_filter)))
cat(sprintf("Pathways selected by LASSO: %d\n", nrow(coef_table)))
cat("\nAUC Summary:\n")
for (i in seq_len(nrow(auc_results))) {
  cat(sprintf("  %-15s AUC = %.3f (%.3f–%.3f)\n",
              auc_results$Cohort[i],
              auc_results$AUC[i],
              auc_results$CI_lower[i],
              auc_results$CI_upper[i]))
}
cat(sprintf("\nCox (Univariate PRS): HR = %.3f (%.3f–%.3f), p = %.4f\n",
            hr_prs, ci_prs[1], ci_prs[2], summary(cox_prs)$coefficients[1, "Pr(>|z|)"]))
cat(sprintf("Log-rank test (PRS high vs low): p = %.4f\n", lr_p))

cat("\n=== Script 16 complete ===\n\n")
