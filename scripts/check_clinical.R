# Check clinical data
df <- read.csv("results/tables/GSE72970_clinical_data.csv", row.names=1)

cat("=== FOLFOX ===\n")
folfox <- df[grepl("FOLFOX", df[["regimen.ch1"]]), ]
cat("Samples:", nrow(folfox), "\n")
cat("PFS non-NA:", sum(!is.na(folfox[["pfs.ch1"]])), "\n")
pfs_vals <- as.numeric(folfox[["pfs.ch1"]])
cat("PFS (months):", paste(round(pfs_vals[1:10], 1), collapse=", "), "\n")
pfs_cens <- as.character(folfox[["pfs.censored.ch1"]])
cat("PFS censored:", paste(pfs_cens[1:10], collapse=", "), "\n")
cat("PFS censored (unique):", paste(unique(pfs_cens), collapse=", "), "\n")
cat("Response status:\n")
print(table(folfox[["response.status.ch1"]], useNA="ifany"))

cat("\n=== FOLFIRI ===\n")
folfiri <- df[grepl("FOLFIRI", df[["regimen.ch1"]]), ]
cat("Samples:", nrow(folfiri), "\n")
cat("PFS non-NA:", sum(!is.na(folfiri[["pfs.ch1"]])), "\n")
cat("Response status:\n")
print(table(folfiri[["response.status.ch1"]], useNA="ifany"))

cat("\n=== GSE69657 ===\n")
df2 <- read.csv("results/tables/GSE69657_clinical_data.csv", row.names=1)
cat("Samples:", nrow(df2), "\n")
cat("Columns:", paste(colnames(df2), collapse=", "), "\n")
cat("Response status:\n")
print(table(df2[["response.status.ch1"]], useNA="ifany"))

# Check for PFS in GSE69657
pfs_cols <- grep("pfs|survival|os|event|follow.up", colnames(df2), value=TRUE, ignore.case=TRUE)
cat("PFS/survival columns:", paste(pfs_cols, collapse=", "), "\n")
