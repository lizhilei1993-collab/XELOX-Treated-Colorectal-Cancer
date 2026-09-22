#!/usr/bin/env Rscript
# Compute SHAP values from XGBoost model and training data
# Output: CSV files for Python to use in generating Figure 2

cat("=== Computing SHAP Values ===\n")

# Paths
BASE <- "/path/to/xelox_project"
MODEL_PATH <- file.path(BASE, "results/tables/ml_phase2/XGBoost_model_phase2.xgb")
TRAIN_RDS <- file.path(BASE, "results/tables/ml_phase2/XELOX_training_data_prepared.rds")
OUT_DIR <- file.path(BASE, "results/tables/ml_phase2")
OUT_SHAP <- file.path(OUT_DIR, "SHAP_values_per_sample.csv")
OUT_BAR <- file.path(OUT_DIR, "SHAP_summary_full.csv")

# Load packages
library(xgboost)
library(SHAPforxgboost)

cat("Loading training data...\n")
train <- readRDS(TRAIN_RDS)
cat(sprintf("  Training data: %d samples\n", nrow(train)))

# Extract feature columns and target
# The fingerprint genes
fingerprint_genes <- c("CSNK1G2", "C5ORF42", "GGT1", "GGT2", "GGTLC1", "GGTLC2", "KAZN", "KLK6", "MGA", "MID2", "SOX11", "TAS2R40", "ZNF451")
target_col <- "group"  # resistant/sensitive

# Check which columns exist
cat("Available columns:\n")
cat(colnames(train)[1:20], "\n")
cat("...\n")

# Try to identify feature columns (numeric, non-target)
numeric_cols <- sapply(train, is.numeric)
feat_cols <- names(numeric_cols)[numeric_cols]
cat(sprintf("  Numeric feature columns: %d\n", length(feat_cols)))

# Check for target column
if (target_col %in% colnames(train)) {
  cat(sprintf("  Target column '%s' found\n", target_col))
  labels <- as.numeric(as.factor(train[[target_col]])) - 1  # 0=Sensitive, 1=Resistant
} else {
  cat(sprintf("  Target column '%s' not found, looking for alternatives...\n", target_col))
  # Check for group column
  for (alt in c("group", "status", "response", "label")) {
    if (alt %in% colnames(train)) {
      labels <- as.numeric(as.factor(train[[alt]])) - 1
      cat(sprintf("    Using '%s' as target\n", alt))
      break
    }
  }
}

# Check feature data
cat(sprintf("  Features: %d cols, Samples: %d\n", ncol(train), nrow(train)))

# Load pre-trained model
cat("Loading XGBoost model...\n")
model <- xgb.load(MODEL_PATH)
cat(sprintf("  Model loaded, %d trees\n", model$nrounds))

# Get feature names from model
model_features <- model$feature_names
cat(sprintf("  Model expects %d features: %s\n", length(model_features), paste(model_features[1:5], collapse=", ")))

# Prepare data matrix
if (all(model_features %in% colnames(train))) {
  x_data <- as.matrix(train[, model_features])
} else {
  # Try to find matching columns with partial names
  cat("  Looking for matching feature columns...\n")
  matched <- character()
  for (feat in model_features) {
    if (feat %in% colnames(train)) {
      matched <- c(matched, feat)
    } else {
      # Try partial matching
      idx <- grep(feat, colnames(train), ignore.case = TRUE)
      if (length(idx) > 0) {
        matched <- c(matched, colnames(train)[idx[1]])
        cat(sprintf("    Partial match: '%s' -> '%s'\n", feat, colnames(train)[idx[1]]))
      }
    }
  }
  if (length(matched) == length(model_features)) {
    x_data <- as.matrix(train[, matched])
    colnames(x_data) <- model_features
  } else {
    cat("ERROR: Cannot match all features!\n")
    cat(sprintf("  Matched: %d / %d\n", length(matched), length(model_features)))
    quit(status = 1)
  }
}

# Compute SHAP values using SHAPforxgboost
cat("Computing SHAP values...\n")

# Prepare for SHAP
dtrain <- xgb.DMatrix(data = x_data, label = labels)

# Compute SHAP values (using xgboost's built-in SHAP)
shap_contrib <- predict(model, dtrain, predcontrib = TRUE, predmargin = TRUE)

# Separate bias and features
bias <- shap_contrib[, 1]  # Last column is bias
shap_features <- shap_contrib[, -1]  # All except bias
colnames(shap_features) <- model_features

cat(sprintf("  SHAP matrix: %d samples x %d features\n", nrow(shap_features), ncol(shap_features)))

# Save per-sample SHAP values
shap_long <- reshape2::melt(shap_features, varnames = c("Sample", "Feature"), value.name = "SHAP_value")
shap_long$Sample_ID <- rownames(x_data)[shap_long$Sample]

# Also add gene expression values
for (i in 1:nrow(shap_long)) {
  feat <- shap_long$Feature[i]
  samp <- shap_long$Sample[i]
  if (feat %in% colnames(x_data)) {
    shap_long$Expression[i] <- x_data[samp, feat]
  }
}

# Save
write.csv(shap_long, OUT_SHAP, row.names = FALSE)
cat(sprintf("  Saved per-sample SHAP values: %s\n", OUT_SHAP))

# Compute summary statistics
shap_summary <- data.frame(
  Feature = model_features,
  MeanAbsSHAP = colMeans(abs(shap_features)),
  MedianSHAP = apply(abs(shap_features), 2, median),
  MeanSHAP = colMeans(shap_features),
  SD_SHAP = apply(shap_features, 2, sd)
)

# Add direction info
shap_summary$Direction <- ifelse(shap_summary$MeanSHAP > 0, "Risk", "Protective")

write.csv(shap_summary, OUT_BAR, row.names = FALSE)
cat(sprintf("  Saved SHAP summary: %s\n", OUT_BAR))

cat("\n=== SHAP Computation Complete ===\n")
cat("Output files:\n")
cat(sprintf("  Per-sample: %s\n", OUT_SHAP))
cat(sprintf("  Summary: %s\n", OUT_BAR))
