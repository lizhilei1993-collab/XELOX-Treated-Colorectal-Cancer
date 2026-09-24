# ============================================================
# 09_ML_SHAP_MODELING.R
# XELOX Resistance - Molecular Fingerprint via ML + SHAP
# ============================================================
# Purpose:
#   - Train XGBoost + LightGBM on GEO (GSE39582) expression data
#   - Identify top predictive genes using SHAP values
#   - Validate on TCGA-COAD/READ independently
#
# Output:
#   results/tables/ml/
#     XGBoost_model.rds, LightGBM_model.rds
#     ML_model_performance.csv
#     SHAP_summary.csv
#     XELOX_molecular_fingerprint.csv
#     ML_training_data.rds
#   results/figures/ml/
#     SHAP_summary_top30.pdf
#     SHAP_beeswarm.pdf
#     ROC_curves.pdf
#     Model_comparison.pdf
# ============================================================

Sys.setenv(TMPDIR = "/tmp", TMP = "/tmp", TEMP = "/tmp")
.libPaths(c("/path/to/Rlibs", .libPaths()))

PROJECT_ROOT     <- "/path/to/xelox_project"
SCRIPT_DIR       <- file.path(PROJECT_ROOT, "scripts")
DATA_GEO_DIR     <- file.path(PROJECT_ROOT, "data", "geo")
DATA_TCGA_DIR    <- file.path(PROJECT_ROOT, "data", "tcga")
DATA_PROC_DIR    <- file.path(PROJECT_ROOT, "data", "processed")
RESULTS_TAB_DIR  <- file.path(PROJECT_ROOT, "results", "tables")
RESULTS_FIG_DIR  <- file.path(PROJECT_ROOT, "results", "figures")
ML_TAB_DIR       <- file.path(RESULTS_TAB_DIR, "ml")
ML_FIG_DIR       <- file.path(RESULTS_FIG_DIR, "ml")

dir.create(ML_TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(ML_FIG_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

cat("============================================================\n")
cat("XELOX Resistance - ML + SHAP Molecular Fingerprint\n")
cat("============================================================\n\n")

# ============================================================
# 0. Load required packages
# ============================================================
cat("=== [0] Loading packages ===\n\n")

required_pkgs <- c("Biobase", "xgboost", "lightgbm", "caret", "pROC",
                    "SHAPforxgboost", "ggplot2", "doParallel")
for (pkg in required_pkgs) {
  tryCatch({
    library(pkg, character.only = TRUE, quietly = TRUE)
    cat(sprintf("  \u2713 %s v%s\n", pkg, packageVersion(pkg)))
  }, error = function(e) {
    cat(sprintf("  \u2717 %s - %s\n", pkg, conditionMessage(e)))
  })
}
cat("\n")

# ============================================================
# 1. Load training data: GSE39582 expression + clinical groups
# ============================================================
cat("=== [1] Loading GSE39582 training data ===\n\n")

# Load expression set
eset_file <- file.path(DATA_GEO_DIR, "GSE39582_eset.rds")
if (!file.exists(eset_file)) stop("GSE39582 eset not found at ", eset_file)
eset <- readRDS(eset_file)
expr_mat <- exprs(eset)
cat(sprintf("  Expression matrix: %d probes x %d samples\n", nrow(expr_mat), ncol(expr_mat)))

# Load clinical groups
groups_file <- file.path(RESULTS_TAB_DIR, "GSE39582_xelox_groups.csv")
groups_df <- read.csv(groups_file, stringsAsFactors = FALSE)
cat(sprintf("  Clinical groups: %d samples\n", nrow(groups_df)))

# Check sample overlap
common_samples <- intersect(colnames(expr_mat), groups_df$sample_id)
cat(sprintf("  Overlap samples: %d\n", length(common_samples)))

# Filter to XELOX-resistant/sensitive only
groups_df <- groups_df[groups_df$sample_id %in% common_samples, ]
groups_df <- groups_df[groups_df$group %in% c("resistant", "sensitive"), ]
cat(sprintf("  Analysis samples: %d (resistant=%d, sensitive=%d)\n",
            nrow(groups_df),
            sum(groups_df$group == "resistant"),
            sum(groups_df$group == "sensitive")))

if (nrow(groups_df) < 20) stop("Too few samples for ML modeling")

# Build expression matrix for these samples
expr_sub <- expr_mat[, groups_df$sample_id, drop = FALSE]
y <- ifelse(groups_df$group == "resistant", 1, 0)  # 1 = resistant, 0 = sensitive
names(y) <- groups_df$sample_id

cat(sprintf("  Class balance: resistant=%d (%.1f%%), sensitive=%d (%.1f%%)\n",
            sum(y == 1), mean(y) * 100,
            sum(y == 0), (1 - mean(y)) * 100))

# ============================================================
# 2. Load gene pool and subset expression
# ============================================================
cat("\n=== [2] Subsetting to gene pool features ===\n\n")

gp_file <- file.path(RESULTS_TAB_DIR, "XELOX_final_gene_pool.csv")
gene_pool <- read.csv(gp_file, stringsAsFactors = FALSE)
gp_genes <- toupper(trimws(gene_pool$Gene))
cat(sprintf("  Gene pool size: %d genes\n", length(gp_genes)))

# CRITICAL: Expression matrix probes (Affy) are lowercase (e.g. "1007_s_at")
# while mapped CSV probes are mixed/uppercase (e.g. "212970_AT").
# Normalize all probe IDs to UPPERCASE for matching.
rownames_orig <- rownames(expr_sub)
rownames_upper <- toupper(rownames_orig)

# Also try loading the DEG mapped files for probe-gene mapping
probe2gene <- list()
deg_mapped_files <- list.files(RESULTS_TAB_DIR, pattern = "DEG_.*_mapped\\.csv$")
for (f in deg_mapped_files) {
  deg <- read.csv(file.path(RESULTS_TAB_DIR, f), stringsAsFactors = FALSE)
  if ("gene_symbol" %in% colnames(deg)) {
    for (i in seq_len(nrow(deg))) {
      probe <- toupper(trimws(deg$probe_id[i]))
      symbol <- toupper(trimws(deg$gene_symbol[i]))
      if (symbol != "" && !is.na(symbol) && is.null(probe2gene[[probe]])) {
        probe2gene[[probe]] <- symbol  # Keep first mapping
      }
    }
  }
}
cat(sprintf("  Probe-to-gene map: %d probes mapped\n", length(probe2gene)))

# Map gene pool genes to probes (reverse lookup: gene -> probe)
gene2probes <- list()
for (probe in names(probe2gene)) {
  symbol <- probe2gene[[probe]]
  if (symbol %in% gp_genes) {
    if (is.null(gene2probes[[symbol]])) {
      gene2probes[[symbol]] <- probe
    }
  }
}
cat(sprintf("  Gene pool -> probes mapped: %d / %d\n", length(gene2probes), length(gp_genes)))

# Select probes that exist in the expression matrix (case-insensitive via UPPER)
selected_probes_upper <- unique(c(unlist(gene2probes)))
selected_idx <- which(rownames_upper %in% selected_probes_upper)
cat(sprintf("  Probes found in expression matrix: %d / %d\n",
            length(selected_idx), length(unique(unlist(gene2probes)))))

# Also try direct probe-as-gene-name matching
probes_as_genes <- rownames_upper
direct_match_idx <- which(probes_as_genes %in% gp_genes)
# Exclude already selected
direct_match_idx <- setdiff(direct_match_idx, selected_idx)
if (length(direct_match_idx) > 0) {
  cat(sprintf("  Direct probe-as-gene-name matches: %d\n", length(direct_match_idx)))
}

# Final selected indices (keep original rownames for subsetting)
selected_idx <- unique(c(selected_idx, direct_match_idx))
cat(sprintf("  Total probes selected: %d / %d\n", length(selected_idx), nrow(expr_sub)))

if (length(selected_idx) < 50) {
  cat("  WARNING: Too few probe-matched features. Using top 2000 variable probes instead.\n")
  probe_var <- apply(expr_sub, 1, var, na.rm = TRUE)
  selected_idx <- order(probe_var, decreasing = TRUE)[1:2000]
  cat(sprintf("  Fallback: using top %d variable probes\n", length(selected_idx)))
}

# Build training matrix (use original rownames for indexing)
selected_probes <- rownames_orig[selected_idx]
X <- t(expr_sub[selected_probes, , drop = FALSE])
cat(sprintf("  Training matrix: %d samples x %d features\n", nrow(X), ncol(X)))

# Save feature names
feature_names <- colnames(X)

# ============================================================
# 3. Train-test split (stratified 70/30)
# ============================================================
cat("\n=== [3] Stratified train-test split (70/30) ===\n\n")

train_idx <- createDataPartition(y, p = 0.7, list = FALSE)
X_train <- X[train_idx, , drop = FALSE]
y_train <- y[train_idx]
X_test  <- X[-train_idx, , drop = FALSE]
y_test  <- y[-train_idx]

cat(sprintf("  Train: %d samples (resistant=%d, sensitive=%d)\n",
            length(y_train), sum(y_train == 1), sum(y_train == 0)))
cat(sprintf("  Test:  %d samples (resistant=%d, sensitive=%d)\n",
            length(y_test), sum(y_test == 1), sum(y_test == 0)))

# ============================================================
# 4. XGBoost Training with CV
# ============================================================
cat("\n=== [4] XGBoost Training ===\n\n")

dtrain <- xgb.DMatrix(data = X_train, label = y_train)
dtest  <- xgb.DMatrix(data = X_test,  label = y_test)

# Handle class imbalance: resistant≈27% vs sensitive≈73%
scale_pos <- sum(y_train == 0) / sum(y_train == 1)

xgb_params <- list(
  objective        = "binary:logistic",
  eval_metric      = "auc",
  max_depth        = 4,
  eta              = 0.05,
  subsample        = 0.8,
  colsample_bytree = 0.8,
  min_child_weight = 3,
  gamma            = 0.1,
  scale_pos_weight = scale_pos,
  nthread          = 4
)

cat("  Training XGBoost with 5-fold CV ...\n")
cat(sprintf("  scale_pos_weight = %.2f (imbalance correction)\n", scale_pos))

xgb_cv <- xgb.cv(
  params       = xgb_params,
  data         = dtrain,
  nrounds      = 500,
  nfold        = 5,
  early_stopping_rounds = 30,
  print_every_n = 50,
  verbose      = 0,
  stratified   = TRUE,
  prediction   = TRUE
)

# Safely get best nrounds (may be NULL if CV didn't converge)
best_nrounds <- xgb_cv$best_iteration
if (is.null(best_nrounds) || length(best_nrounds) == 0) {
  best_nrounds <- which.max(xgb_cv$evaluation_log$test_auc_mean)
  cat("  CV did not explicitly converge; using best AUC iteration\n")
}
cv_auc <- max(xgb_cv$evaluation_log$test_auc_mean)
cat(sprintf("  Best nrounds: %d | CV AUC: %.4f\n", best_nrounds, cv_auc))

# Train final model
xgb_model <- xgb.train(
  params       = xgb_params,
  data         = dtrain,
  nrounds      = best_nrounds,
  verbose      = 0
)

# Predict on test set
xgb_pred <- predict(xgb_model, dtest)
xgb_roc <- roc(y_test, xgb_pred, quiet = TRUE)
xgb_auc <- as.numeric(xgb_roc$auc)
xgb_threshold <- coords(xgb_roc, "best", ret = "threshold")$threshold[1]
xgb_class <- ifelse(xgb_pred >= xgb_threshold, 1, 0)

cat(sprintf("  Test AUC: %.4f | Best threshold: %.3f\n", xgb_auc, xgb_threshold))
cat(sprintf("  Test accuracy: %.2f%%\n", mean(xgb_class == y_test) * 100))

# ============================================================
# 5. LightGBM Training with CV
# ============================================================
cat("\n=== [5] LightGBM Training ===\n\n")

lgb_train <- lgb.Dataset(data = X_train, label = y_train)

lgb_params <- list(
  objective        = "binary",
  metric           = "auc",
  boosting         = "gbdt",
  num_leaves       = 16,
  learning_rate    = 0.05,
  feature_fraction = 0.8,
  bagging_fraction = 0.8,
  bagging_freq     = 5,
  min_data_in_leaf = 5,
  num_threads      = 4,
  verbose          = -1
)

cat("  Training LightGBM with 5-fold CV ...\n")

# CV
lgb_cv <- lgb.cv(
  params       = lgb_params,
  data         = lgb_train,
  nrounds      = 500,
  nfold        = 5,
  early_stopping_rounds = 30,
  stratified   = TRUE,
  eval_freq    = 50,
  verbose      = -1
)

lgb_best_rounds <- lgb_cv$best_iter
lgb_cv_auc <- max(lgb_cv$record_evals$valid$auc$eval)
cat(sprintf("  Best nrounds: %d | CV AUC: %.4f\n", lgb_best_rounds, lgb_cv_auc))

# Train final model
lgb_model <- lightgbm(
  params       = lgb_params,
  data         = lgb_train,
  nrounds      = lgb_best_rounds,
  verbose      = -1
)

# Predict on test set
lgb_pred <- predict(lgb_model, X_test)
lgb_roc <- roc(y_test, lgb_pred, quiet = TRUE)
lgb_auc <- as.numeric(lgb_roc$auc)
lgb_threshold <- coords(lgb_roc, "best", ret = "threshold")$threshold[1]
lgb_class <- ifelse(lgb_pred >= lgb_threshold, 1, 0)

cat(sprintf("  Test AUC: %.4f | Best threshold: %.3f\n", lgb_auc, lgb_threshold))
cat(sprintf("  Test accuracy: %.2f%%\n", mean(lgb_class == y_test) * 100))

# ============================================================
# 6. Model Performance Summary
# ============================================================
cat("\n=== [6] Model Performance Summary ===\n\n")

perf <- data.frame(
  Model = c("XGBoost", "LightGBM"),
  CV_AUC = sprintf("%.4f", c(cv_auc, lgb_cv_auc)),
  Test_AUC = sprintf("%.4f", c(xgb_auc, lgb_auc)),
  Test_Accuracy = sprintf("%.1f%%", c(mean(xgb_class == y_test) * 100,
                                       mean(lgb_class == y_test) * 100)),
  Test_Sensitivity = sprintf("%.1f%%", c(
    sum(xgb_class == 1 & y_test == 1) / sum(y_test == 1) * 100,
    sum(lgb_class == 1 & y_test == 1) / sum(y_test == 1) * 100)),
  Test_Specificity = sprintf("%.1f%%", c(
    sum(xgb_class == 0 & y_test == 0) / sum(y_test == 0) * 100,
    sum(lgb_class == 0 & y_test == 0) / sum(y_test == 0) * 100)),
  stringsAsFactors = FALSE
)
print(perf, row.names = FALSE)
write.csv(perf, file.path(ML_TAB_DIR, "ML_model_performance.csv"), row.names = FALSE)
cat("  Performance saved to ML_model_performance.csv\n")

# ============================================================
# 7. ROC Curves
# ============================================================
cat("\n=== [7] Plotting ROC Curves ===\n\n")

pdf(file.path(ML_FIG_DIR, "ROC_curves.pdf"), width = 8, height = 7)
plot(xgb_roc, col = "#E41A1C", lwd = 2.5, main = "ROC Curves - XELOX Resistance Prediction",
     xlim = c(1, 0), ylim = c(0, 1))
plot(lgb_roc, col = "#377EB8", lwd = 2.5, add = TRUE)
legend("bottomright", 
       legend = c(sprintf("XGBoost (AUC = %.4f)", xgb_auc),
                  sprintf("LightGBM (AUC = %.4f)", lgb_auc)),
       col = c("#E41A1C", "#377EB8"), lwd = 2.5)
dev.off()
cat(sprintf("  ROC curves saved\n"))

# ============================================================
# 8. SHAP Analysis (XGBoost primary)
# ============================================================
cat("\n=== [8] SHAP Analysis (XGBoost) ===\n\n")

# SHAP values using xgboost's built-in SHAP
cat("  Computing SHAP values ...\n")
shap_contrib <- predict(xgb_model, rbind(X_train, X_test), predcontrib = TRUE)
shap_matrix <- shap_contrib[, 1:(ncol(shap_contrib) - 1)]  # exclude BIAS column
colnames(shap_matrix) <- feature_names

# Calculate mean |SHAP| per feature (global importance)
shap_importance <- data.frame(
  Feature = feature_names,
  MeanAbsSHAP = colMeans(abs(shap_matrix)),
  stringsAsFactors = FALSE
)
shap_importance <- shap_importance[order(shap_importance$MeanAbsSHAP, decreasing = TRUE), ]
rownames(shap_importance) <- NULL

# Map probes to gene symbols if available
shap_importance$GeneSymbol <- sapply(shap_importance$Feature, function(p) {
  p_upper <- toupper(trimws(p))
  if (!is.null(probe2gene[[p_upper]])) probe2gene[[p_upper]] else p_upper
})

# Top 50 features by SHAP
top50_shap <- head(shap_importance, 50)
cat(sprintf("  SHAP importance computed: %d features\n", nrow(shap_importance)))
cat("  Top 10 features by |SHAP|:\n")
for (i in 1:min(10, nrow(top50_shap))) {
  cat(sprintf("    %d. %s (%s) |SHAP|=%.4f\n",
              i, top50_shap$Feature[i], top50_shap$GeneSymbol[i], top50_shap$MeanAbsSHAP[i]))
}

write.csv(shap_importance, file.path(ML_TAB_DIR, "SHAP_summary.csv"), row.names = FALSE)
cat("  SHAP summary saved to SHAP_summary.csv\n")

# ============================================================
# 9. SHAP Visualization
# ============================================================
cat("\n=== [9] SHAP Visualization ===\n\n")

# Top 30 SHAP summary bar plot
top30 <- head(shap_importance, 30)
top30$Label <- sprintf("%s (%s)", top30$Feature, top30$GeneSymbol)

p_bar <- ggplot(top30, aes(x = reorder(Label, MeanAbsSHAP), y = MeanAbsSHAP)) +
  geom_bar(stat = "identity", fill = "#1a56db", alpha = 0.85) +
  coord_flip() +
  labs(title = "Top 30 Features by Mean |SHAP|",
       subtitle = "XGBoost Model - XELOX Resistance Prediction",
       x = "Feature (Probe_Gene)", y = "Mean |SHAP|") +
  theme_minimal(base_size = 10) +
  theme(plot.margin = margin(10, 20, 10, 10))

ggsave(file.path(ML_FIG_DIR, "SHAP_summary_top30.pdf"), p_bar, width = 12, height = 10)
cat("  SHAP summary bar plot saved\n")

# SHAP beeswarm plot (top 20 features)
top20_features <- head(shap_importance$Feature, 20)
shap_long <- shap.prep(X = as.data.frame(rbind(X_train, X_test)),
                       shap = shap_matrix,
                       top_n = 20)

if (nrow(shap_long) > 0) {
  p_swarm <- plot.shap.summary(
    data_long = shap_long,
    scientific = TRUE,
    x_bound = NULL,
    dilute = 5,
    digits = 4
  ) +
    ggtitle("SHAP Beeswarm - Top 20 Features") +
    theme_minimal(base_size = 10)
  
  ggsave(file.path(ML_FIG_DIR, "SHAP_beeswarm.pdf"), p_swarm, width = 12, height = 8)
  cat("  SHAP beeswarm plot saved\n")
}

# ============================================================
# 10. XELOX Molecular Fingerprint (Top N genes)
# ============================================================
cat("\n=== [10] XELOX Molecular Fingerprint Definition ===\n\n")

# Select fingerprint genes: top features consistently important
fingerprint_cutoff <- min(50, nrow(shap_importance))
fingerprint <- data.frame(
  Rank = 1:fingerprint_cutoff,
  Feature = shap_importance$Feature[1:fingerprint_cutoff],
  GeneSymbol = shap_importance$GeneSymbol[1:fingerprint_cutoff],
  MeanAbsSHAP = shap_importance$MeanAbsSHAP[1:fingerprint_cutoff],
  RelativeImportance = shap_importance$MeanAbsSHAP[1:fingerprint_cutoff] /
    max(shap_importance$MeanAbsSHAP[1:fingerprint_cutoff]) * 100,
  stringsAsFactors = FALSE
)

cat(sprintf("  Molecular Fingerprint: %d genes\n", nrow(fingerprint)))
cat("  Top 10 fingerprint genes:\n")
for (i in 1:min(10, nrow(fingerprint))) {
  cat(sprintf("    %d. %s (rel.imp=%.1f%%)\n",
              i, fingerprint$GeneSymbol[i], fingerprint$RelativeImportance[i]))
}

write.csv(fingerprint, file.path(ML_TAB_DIR, "XELOX_molecular_fingerprint.csv"), row.names = FALSE)
cat("  Fingerprint saved to XELOX_molecular_fingerprint.csv\n")

# ============================================================
# 11. TCGA External Validation (if available)
# ============================================================
cat("\n=== [11] TCGA External Validation ===\n\n")

tcga_file <- file.path(DATA_TCGA_DIR, "TCGA_COAD_READ_gene_pool_expression.csv")

if (file.exists(tcga_file)) {
  cat("  TCGA validation data found.\n")
  
  tcga_data <- read.csv(tcga_file, stringsAsFactors = TRUE)
  cat(sprintf("  TCGA dataset: %d samples x %d columns\n", nrow(tcga_data), ncol(tcga_data)))
  
  # Check response column
  if ("response" %in% colnames(tcga_data)) {
    tcga_tab <- table(tcga_data$response)
    cat("  Response distribution:\n")
    print(tcga_tab)
    
    # Subset to resistant/sensitive
    tcga_df <- tcga_data[tcga_data$response %in% c("resistant", "sensitive"), ]
    tcga_y <- ifelse(tcga_df$response == "resistant", 1, 0)
    cat(sprintf("  Analysis samples: %d\n", nrow(tcga_df)))
    
    # Extract features (remove metadata columns)
    meta_cols <- c("patient_barcode", "sample_barcode", "response", "cancer_type",
                   grep("^TCGA", colnames(tcga_df), value = TRUE))
    skip_cols <- intersect(meta_cols, colnames(tcga_df))
    tcga_feat_cols <- setdiff(colnames(tcga_df), skip_cols)
    
    cat(sprintf("  Feature columns in TCGA: %d\n", length(tcga_feat_cols)))
    
    # Common features with training
    common_features <- intersect(feature_names, tcga_feat_cols)
    cat(sprintf("  Common features with model: %d / %d\n",
                length(common_features), length(feature_names)))
    
    if (length(common_features) >= 10 && nrow(tcga_df) >= 10) {
      # Prepare TCGA matrix
      X_tcga <- as.matrix(tcga_df[, common_features, drop = FALSE])
      
      # Predict
      dtcga <- xgb.DMatrix(data = X_tcga)
      tcga_pred <- predict(xgb_model, dtcga)
      
      # Evaluate
      tcga_roc <- roc(tcga_y, tcga_pred, quiet = TRUE)
      tcga_auc <- as.numeric(tcga_roc$auc)
      tcga_class <- ifelse(tcga_pred >= xgb_threshold, 1, 0)
      
      cat(sprintf("\n  === TCGA External Validation ===\n"))
      cat(sprintf("  AUC: %.4f\n", tcga_auc))
      cat(sprintf("  Accuracy: %.2f%%\n", mean(tcga_class == tcga_y) * 100))
      
      # Save TCGA validation results
      tcga_val_df <- data.frame(
        sample_barcode = tcga_df$sample_barcode,
        true_label = tcga_df$response,
        predicted_prob = tcga_pred,
        predicted_class = ifelse(tcga_class == 1, "resistant", "sensitive"),
        stringsAsFactors = FALSE
      )
      write.csv(tcga_val_df, file.path(ML_TAB_DIR, "TCGA_validation_predictions.csv"),
                row.names = FALSE)
      cat("  TCGA predictions saved\n")
      
      # Update performance table
      perf$TCGA_Validation_AUC <- sprintf("%.4f", c(tcga_auc, NA))
    } else {
      cat("  WARNING: Insufficient features or samples for TCGA validation\n")
    }
  } else {
    cat("  WARNING: TCGA data missing 'response' column\n")
  }
} else {
  cat("  TCGA validation data not available (data still downloading?)\n")
  cat("  Will check again at end of script ...\n")
}

# ============================================================
# 12. Save Models and Training Data
# ============================================================
cat("\n=== [12] Saving models and data ===\n\n")

# Save XGBoost model
xgb.save(xgb_model, file.path(ML_TAB_DIR, "XGBoost_model.xgb"))
cat("  XGBoost model saved\n")

# Save LightGBM model
lgb.save(lgb_model, file.path(ML_TAB_DIR, "LightGBM_model.txt"))
cat("  LightGBM model saved\n")

# Save training data for reproducibility
train_data <- list(
  X_train = X_train, y_train = y_train,
  X_test = X_test, y_test = y_test,
  feature_names = feature_names,
  probe2gene = probe2gene,
  gene_pool = gp_genes,
  selected_probes = selected_probes,
  xgb_threshold = xgb_threshold,
  best_nrounds = best_nrounds,
  xgb_params = xgb_params
)
saveRDS(train_data, file.path(ML_TAB_DIR, "ML_training_data.rds"))
cat("  Training data saved\n")

# Save performance table
write.csv(perf, file.path(ML_TAB_DIR, "ML_model_performance.csv"), row.names = FALSE)

# ============================================================
# 13. Model Comparison Plot
# ============================================================
cat("\n=== [13] Model Comparison ===\n\n")

pdf(file.path(ML_FIG_DIR, "Model_comparison.pdf"), width = 10, height = 7)

par(mfrow = c(2, 2))

# ROC curves
plot(xgb_roc, col = "#E41A1C", lwd = 2, main = "Test Set ROC",
     xlim = c(1, 0), ylim = c(0, 1))
plot(lgb_roc, col = "#377EB8", lwd = 2, add = TRUE)
legend("bottomright", 
       legend = c(sprintf("XGB (AUC=%.4f)", xgb_auc),
                  sprintf("LGB (AUC=%.4f)", lgb_auc)),
       col = c("#E41A1C", "#377EB8"), lwd = 2)

# SHAP top 10 importance
top10 <- head(shap_importance, 10)
barplot(top10$MeanAbsSHAP, names.arg = top10$GeneSymbol,
        las = 2, col = "#1a56db", main = "Top 10 Features (SHAP)",
        ylab = "Mean |SHAP|", cex.names = 0.8)

# Prediction distribution - XGBoost
hist(xgb_pred[y_test == 0], col = rgb(0, 0.8, 0, 0.5), 
     xlim = c(0, 1), breaks = 20,
     main = "XGBoost Prediction Distribution",
     xlab = "Predicted probability (resistant)", ylab = "Count")
hist(xgb_pred[y_test == 1], col = rgb(0.8, 0, 0, 0.5), 
     breaks = 20, add = TRUE)
legend("topright", legend = c("Sensitive", "Resistant"),
       fill = c(rgb(0, 0.8, 0, 0.5), rgb(0.8, 0, 0, 0.5)), cex = 0.8)

# CV scores
cv_scores <- data.frame(
  Iteration = xgb_cv$evaluation_log$iter,
  CV_AUC = xgb_cv$evaluation_log$test_auc_mean,
  CV_AUC_SD = xgb_cv$evaluation_log$test_auc_std
)
plot(cv_scores$Iteration, cv_scores$CV_AUC, type = "l", col = "#1a56db", lwd = 2,
     main = "XGBoost CV AUC", xlab = "Iteration", ylab = "AUC",
     ylim = range(c(cv_scores$CV_AUC - cv_scores$CV_AUC_SD,
                     cv_scores$CV_AUC + cv_scores$CV_AUC_SD)))
polygon(c(cv_scores$Iteration, rev(cv_scores$Iteration)),
        c(cv_scores$CV_AUC - cv_scores$CV_AUC_SD,
          rev(cv_scores$CV_AUC + cv_scores$CV_AUC_SD)),
        col = rgb(26, 86, 219, 50, maxColorValue = 255), border = NA)
abline(v = best_nrounds, lty = 2, col = "gray")
text(best_nrounds, max(cv_scores$CV_AUC), 
     sprintf("Best=%d\nAUC=%.4f", best_nrounds, cv_auc),
     pos = 4, cex = 0.8)

dev.off()
cat("  Model comparison plot saved\n")

# ============================================================
# 14. Summary Report
# ============================================================
cat("\n=== [14] ML/SHAP Summary ===\n\n")

summary_lines <- c(
  "==============================================================================",
  "ML + SHAP Molecular Fingerprint - Summary",
  "XELOX Resistance Study",
  paste("Date:", Sys.Date()),
  "==============================================================================",
  "",
  sprintf("Training samples: %d (GSE39582, XELOX-like)", nrow(X)),
  sprintf("  Resistant: %d (%.1f%%)", sum(y == 1), mean(y) * 100),
  sprintf("  Sensitive: %d (%.1f%%)", sum(y == 0), (1 - mean(y)) * 100),
  sprintf("Features: %d (gene pool probes)", length(selected_probes)),
  "",
  "--- Model Performance ---",
  sprintf("  XGBoost CV AUC: %.4f", cv_auc),
  sprintf("  XGBoost Test AUC: %.4f", xgb_auc),
  sprintf("  LightGBM CV AUC: %.4f", lgb_cv_auc),
  sprintf("  LightGBM Test AUC: %.4f", lgb_auc),
  "",
  "--- Top 10 Fingerprint Genes ---"
)

for (i in 1:min(10, nrow(fingerprint))) {
  summary_lines <- c(summary_lines,
    sprintf("  %d. %s |SHAP|=%.4f (rel.imp=%.1f%%)",
            i, fingerprint$GeneSymbol[i],
            fingerprint$MeanAbsSHAP[i],
            fingerprint$RelativeImportance[i]))
}

summary_lines <- c(summary_lines,
  "",
  "--- Output Files ---",
  sprintf("  Models: %s", file.path(ML_TAB_DIR, "XGBoost_model.xgb")),
  sprintf("  SHAP: %s", file.path(ML_TAB_DIR, "SHAP_summary.csv")),
  sprintf("  Fingerprint: %s", file.path(ML_TAB_DIR, "XELOX_molecular_fingerprint.csv")),
  sprintf("  Performance: %s", file.path(ML_TAB_DIR, "ML_model_performance.csv")),
  sprintf("  Figures: %s", ML_FIG_DIR),
  "",
  "=============================================================================="
)

writeLines(summary_lines, file.path(ML_TAB_DIR, "ML_SHAP_summary_report.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

cat("\n=== ML + SHAP Modeling Complete ===\n")
