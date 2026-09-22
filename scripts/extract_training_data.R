#!/usr/bin/env Rscript
# Extract training data from RDS list structure
cat("=== Extracting Training Data ===\n")

BASE <- "/path/to/xelox_project"
RDS <- file.path(BASE, "results/tables/ml_phase2/XELOX_training_data_prepared.rds")
OUT_TRAIN <- file.path(BASE, "results/tables/ml_phase2/training_data_for_shap.csv")

# Load data
train_list <- readRDS(RDS)
cat("Loaded RDS list with", length(train_list), "items\n")

X <- train_list$X  # Matrix: 164 samples x 2528 features
y <- train_list$y  # Labels: 164 numeric
features <- train_list$features  # Feature names
groups_df <- train_list$groups_df

cat(sprintf("X: %d x %d\n", nrow(X), ncol(X)))
cat(sprintf("y: %d values\n", length(y)))
cat(sprintf("Features: %d\n", length(features)))

# Assign feature names
colnames(X) <- features

# Get the 10 fingerprint genes
fingerprint <- c("CSNK1G2", "C5ORF42", "GGT1 /// GGT2 /// GGTLC1 /// GGTLC2", 
                 "KAZN", "KLK6", "MGA", "MID2", "SOX11", "TAS2R40", "ZNF451")
cat("Looking for fingerprint features in data...\n")
for (f in fingerprint) {
  found <- f %in% features
  cat(sprintf("  %s: %s\n", f, ifelse(found, "FOUND", "NOT FOUND")))
}

# Extract just the 10 fingerprint features
x_sub <- X[, fingerprint, drop = FALSE]
cat(sprintf("\nSubset matrix: %d x %d\n", nrow(x_sub), ncol(x_sub)))

# Save to CSV
train_df <- as.data.frame(x_sub)
train_df$label <- y
train_df$group <- groups_df$group[1:nrow(train_df)]

write.csv(train_df, OUT_TRAIN, row.names = TRUE)
cat(sprintf("Saved: %s\n", OUT_TRAIN))
cat("=== Done ===\n")
