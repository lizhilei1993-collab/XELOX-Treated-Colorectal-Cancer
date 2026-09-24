# Inspect GSE39582 characteristics structure
lib_path <- "/path/to/Rlibs"
.libPaths(c(lib_path, .libPaths()))
library(GEOquery)
gse <- readRDS("/path/to/xelox_project/data/geo/GSE39582_eset.rds")
pheno <- pData(gse)

# Find all characteristics columns
char_cols <- grep("characteristics_ch1", colnames(pheno), value = TRUE)
cat("Total characteristics columns:", length(char_cols), "\n")
cat("Column names:\n")
cat(paste(char_cols, collapse = "\n"), "\n")

# Check what's in characteristics_ch1 (the only column the old code uses)
cat("\n--- characteristics_ch1 first 10 values ---\n")
cat(paste(head(pheno[["characteristics_ch1"]], 10), collapse = "\n---\n"), "\n")

# Check column characteristics_ch1.1
if ("characteristics_ch1.1" %in% colnames(pheno)) {
  cat("\n--- characteristics_ch1.1 first 5 values ---\n")
  cat(paste(head(pheno[["characteristics_ch1.1"]], 5), collapse = "\n---\n"), "\n")
}

# Check column characteristics_ch1.2
if ("characteristics_ch1.2" %in% colnames(pheno)) {
  cat("\n--- characteristics_ch1.2 first 5 values ---\n")
  cat(paste(head(pheno[["characteristics_ch1.2"]], 5), collapse = "\n---\n"), "\n")
}

# Check last column
last_col <- char_cols[length(char_cols)]
cat("\n---", last_col, "first 5 values ---\n")
cat(paste(head(pheno[[last_col]], 5), collapse = "\n---\n"), "\n")

# How many columns have non-NA data
non_na_counts <- sapply(char_cols, function(col) sum(!is.na(pheno[[col]])))
cat("\nNon-NA counts per column:\n")
print(non_na_counts)

# Count rows
cat("\nTotal rows in pheno data:", nrow(pheno), "\n")
