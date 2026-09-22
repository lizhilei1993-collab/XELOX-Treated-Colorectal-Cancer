# Extract all key names from GSE39582 characteristics columns
lib_path <- "C:/Rlibs"
.libPaths(c(lib_path, .libPaths()))
library(GEOquery)

gse <- readRDS("C:/xelox_research/data/geo/GSE39582_eset.rds")
pheno <- pData(gse)

char_cols <- grep("characteristics_ch1", colnames(pheno), value = TRUE)

cat("All characteristic keys found:\n")
for (col in char_cols) {
  vals <- pheno[[col]]
  # Extract key (before colon)
  keys <- trimws(gsub(":.*", "", vals[1]))
  cat(sprintf("  %s -> %s (e.g. %s)\n", col, keys, vals[1]))
}

cat("\n\nAlso check if there are any ':ch1' style columns with clinical data:\n")
ch1_cols <- grep(":ch1$", colnames(pheno), value = TRUE)
cat(paste(ch1_cols, collapse = "\n"), "\n")

# Check for rfs.event etc in any column
cat("\nSearching for rfs, os, chemo in column names:\n")
target_cols <- grep("rfs|os|chemo|tnm|kras|braf|mmr|sex|age|tumor|stage|adjuvant", 
                     colnames(pheno), ignore.case = TRUE, value = TRUE)
cat(paste(target_cols, collapse = "\n"), "\n")
