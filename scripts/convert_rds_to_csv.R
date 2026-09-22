#!/usr/bin/env Rscript
# Convert RDS training data to CSV for Python
cat("Converting RDS to CSV...\n")

BASE <- "/path/to/xelox_project"
RDS <- file.path(BASE, "results/tables/ml_phase2/XELOX_training_data_prepared.rds")
OUT <- file.path(BASE, "results/tables/ml_phase2/training_data.csv")

train <- readRDS(RDS)
cat(sprintf("Loaded: %d samples, %d columns\n", nrow(train), ncol(train)))
cat("Column names:\n")
print(colnames(train)[1:30])

# Save to CSV
write.csv(train, OUT, row.names = FALSE)
cat(sprintf("Saved: %s\n", OUT))
