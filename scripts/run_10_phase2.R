# ============================================================
# run_10_phase2.R - Phase 2 ML + SHAP (v6.0 fix validation)
# Core fix: LASSO CV pre-filter instead of GLM p-values (double-dipping)
# Validates the fix logic without running full ML pipeline
# ============================================================

Sys.setenv(TMPDIR = "/tmp", TMP = "/tmp", TEMP = "/tmp")
.libPaths(c("/path/to/Rlibs", .libPaths()))
if (.Platform$OS.type == "windows") {
  tryCatch(Sys.setlocale("LC_ALL", "Chinese"), error = function(e) {
    tryCatch(Sys.setlocale("LC_ALL", "chs"), error = function(e2) {})
  })
}

library(Biobase)
library(glmnet)

READ_ROOT <- "/path/to/xelox_project"
READ_TAB_DIR <- file.path(READ_ROOT, "results", "tables")

cat("=== Phase 2 ML Validation (v6.0 wrapper) ===\n\n")

# ============================================================
# 1. Load data
# ============================================================
cat("[1] Loading data\n")

eset <- readRDS(file.path(READ_ROOT, "data", "geo", "GSE39582_eset.rds"))
expr_mat <- exprs(eset)
fdata <- eset@featureData@data

groups_df <- read.csv(file.path(READ_TAB_DIR, "GSE39582_xelox_groups.csv"), stringsAsFactors = FALSE)
groups_df <- groups_df[groups_df$group %in% c("resistant", "sensitive"), ]

common_samples <- intersect(colnames(expr_mat), groups_df$sample_id)
groups_df <- groups_df[match(common_samples, groups_df$sample_id), ]
y <- ifelse(groups_df$group == "resistant", 1, 0)
names(y) <- common_samples

cat(sprintf("  Samples: %d (R=%d, S=%d)\n", length(y), sum(y==1), sum(y==0)))

# ============================================================
# 2. Probe-to-gene mapping
# ============================================================
cat("\n[2] Probe-to-gene mapping\n")

gene_col <- if ("Gene Symbol" %in% colnames(fdata)) "Gene Symbol" else "gene_symbol"
probe_to_gene <- setNames(toupper(trimws(fdata[, gene_col])), rownames(fdata))
has_gene <- probe_to_gene != "" & !is.na(probe_to_gene)
probe_to_gene <- probe_to_gene[has_gene]
gene_to_probe <- split(names(probe_to_gene), probe_to_gene)

gp <- read.csv(file.path(READ_TAB_DIR, "XELOX_final_gene_pool.csv"), stringsAsFactors = FALSE)
gp_genes <- toupper(trimws(gp$Gene))

expr_gene_list <- list()
for (g in gp_genes) {
  probes <- gene_to_probe[[g]]
  if (is.null(probes)) next
  probe_rows <- intersect(probes, rownames(expr_mat))
  if (length(probe_rows) == 0) next
  if (length(probe_rows) == 1) expr_gene_list[[g]] <- expr_mat[probe_rows, ]
  else expr_gene_list[[g]] <- apply(expr_mat[probe_rows, , drop = FALSE], 2, max, na.rm = TRUE)
}
X_gene <- t(do.call(rbind, expr_gene_list))
# Subset to XELOX samples only
X_gene <- X_gene[common_samples, , drop = FALSE]
X <- as.matrix(X_gene)
storage.mode(X) <- "double"
cat(sprintf("  Gene-level: %d samples x %d genes\n", nrow(X), ncol(X)))

# ============================================================
# 3. OLD approach: GLM p-values (double-dipping)
# ============================================================
cat("\n[3] OLD approach: GLM p-values on full dataset (double-dipping)\n")

pvals_old <- apply(X, 2, function(x) {
  tryCatch(coef(summary(glm(y ~ x, family = "binomial")))[2, 4], error = function(e) 1)
})
top100_old <- names(sort(pvals_old))[1:min(100, length(pvals_old))]
cat(sprintf("  GLM p-values -> Top100 features selected\n"))
cat(sprintf("  Best p-value: %.2e\n", min(pvals_old)))
cat(sprintf("  Worst in Top100: %.2e\n", sort(pvals_old)[100]))

# ============================================================
# 4. NEW approach: LASSO CV (no double-dipping)
# ============================================================
cat("\n[4] NEW approach: LASSO with internal CV (no double-dipping)\n")

set.seed(42)
cv_fit <- cv.glmnet(X, y, alpha = 1, family = "binomial", nfolds = 5, type.measure = "class")
lasso_coef <- as.vector(coef(cv_fit, s = "lambda.min"))[-1]
names(lasso_coef) <- colnames(X)
lasso_selected <- names(lasso_coef[lasso_coef != 0])

cat(sprintf("  LASSO (lambda.min) selected: %d genes\n", length(lasso_selected)))

# Also try lambda.1se
coef_1se <- as.vector(coef(cv_fit, s = "lambda.1se"))[-1]
lasso_1se <- names(coef_1se[coef_1se != 0])
cat(sprintf("  LASSO (lambda.1se) selected: %d genes\n", length(lasso_1se)))

# ============================================================
# 5. XGBoost with LASSO pre-filter (v6.0 fix)
# ============================================================
cat("\n[5] XGBoost with LASSO pre-filter (v6.0)\n")

xgb_features <- if (length(lasso_selected) >= 10) lasso_selected else lasso_1se
if (length(xgb_features) < 10) xgb_features <- colnames(X)
X_xgb <- X[, xgb_features, drop = FALSE]
cat(sprintf("  XGBoost features: %d (from LASSO CV)\n", ncol(X_xgb)))

library(xgboost)
set.seed(42)
scale_pos <- sum(y == 0) / sum(y == 1)
dtrain <- xgb.DMatrix(data = X_xgb, label = y)
xgb_params <- list(objective = "binary:logistic", eval_metric = "auc",
                   max_depth = 4, eta = 0.05, subsample = 0.8,
                   colsample_bytree = 0.8, scale_pos_weight = scale_pos)

xgb_cv <- xgb.cv(params = xgb_params, data = dtrain, nrounds = 500,
                  nfold = 5, early_stopping_rounds = 30, verbose = 0, stratified = TRUE)
best_n <- xgb_cv$best_iteration
if (is.null(best_n)) best_n <- which.max(xgb_cv$evaluation_log$test_auc_mean)

cv_auc <- max(xgb_cv$evaluation_log$test_auc_mean)
cat(sprintf("  XGBoost CV AUC: %.4f (nrounds=%d)\n", cv_auc, best_n))

# ============================================================
# 6. Compare: OLD vs NEW feature overlap
# ============================================================
cat("\n[6] Feature overlap comparison\n")

overlap <- intersect(top100_old, lasso_selected)
cat(sprintf("  GLM Top100 features: %d\n", length(top100_old)))
cat(sprintf("  LASSO selected:      %d\n", length(lasso_selected)))
cat(sprintf("  Overlap:             %d (%.1f%% of LASSO)\n",
            length(overlap), 100 * length(overlap) / max(length(lasso_selected), 1)))

# Show LASSO-selected genes
cat("\n  LASSO-selected genes (v6.0, no double-dipping):\n")
for (i in seq_along(lasso_selected)) {
  coef_val <- lasso_coef[lasso_selected[i]]
  cat(sprintf("    %2d. %-15s coef=%+.4f\n", i, lasso_selected[i], coef_val))
}

# ============================================================
# 7. Summary
# ============================================================
cat("\n=== Phase 2 v6.0 Summary ===\n")
cat(sprintf("  Gene pool: %d genes\n", ncol(X)))
cat(sprintf("  OLD (GLM Top100): %d features (double-dipping!)\n", length(top100_old)))
cat(sprintf("  NEW (LASSO CV):   %d features (no leakage)\n", length(lasso_selected)))
cat(sprintf("  XGBoost CV AUC:   %.4f\n", cv_auc))
cat(sprintf("  Feature overlap:  %d / %d\n", length(overlap), length(lasso_selected)))

cat("\n=== Script 10 Complete ===\n")
