# Compare :ch1 columns vs characteristics columns
lib_path <- "/path/to/Rlibs"
.libPaths(c(lib_path, .libPaths()))
library(GEOquery)

gse <- readRDS("/path/to/xelox_project/data/geo/GSE39582_eset.rds")
pheno <- pData(gse)

# Check what :ch1 columns contain (just values, no key prefix?)
ch1_cols <- grep(":ch1$", colnames(pheno), value = TRUE)
cat("Sample values from :ch1 columns:\n\n")

for (col in head(ch1_cols, 5)) {
  vals <- head(pheno[[col]], 5)
  cat(sprintf("  %s: [%s]\n", col, paste(vals, collapse = "], [")))
}

# Compare with raw characteristics
cat("\nCorresponding characteristics_ch1 values:\n\n")
# characteristics_ch1.X -> lookup which one matches each :ch1 column
for (col in head(ch1_cols, 5)) {
  # Find matching characteristics column
  key_name <- sub(":ch1$", "", col)
  # Search in raw chars
  match_val <- NA
  for (cc in grep("^characteristics_ch1", colnames(pheno), value = TRUE)) {
    v <- pheno[[cc]][1]
    if (grepl(paste0("^", key_name, ":"), v, ignore.case = TRUE)) {
      match_val <- v
      break
    }
  }
  vals <- head(pheno[[col]], 3)
  cat(sprintf("  :ch1 col '%s' values: [%s]\n", col, paste(vals, collapse = "], [")))
  cat(sprintf("  raw char sample: %s\n\n", match_val))
}

cat("\n\nVerification: are :ch1 values always the value part (no key: prefix)?\n")
verification_cols <- grep("rfs.event|rfs.delay|os.event|os.delay|chemotherapy|Sex|age", 
                          ch1_cols, value = TRUE, ignore.case = TRUE)
for (col in verification_cols) {
  r5 <- head(pheno[[col]], 5)
  cat(sprintf("  %s: [%s]\n", col, paste(r5, collapse = "], [")))
}
