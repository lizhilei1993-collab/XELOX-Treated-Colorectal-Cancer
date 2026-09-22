# 28b_parse_gse83129.R
# 解析 GSE83129 表达矩阵和临床数据
Sys.setenv(TMPDIR = "C:/temp", TMP = "C:/temp", TEMP = "C:/temp")
.libPaths(c("C:/Rlibs", .libPaths()))

PROJECT_ROOT <- "C:/xelox_research"
GEO_DIR <- file.path(PROJECT_ROOT, "data", "geo")

library(GEOquery)

# 从本地文件解析系列矩阵
gzfile <- file.path(GEO_DIR, "GSE83129_series_matrix.txt.gz")

gse <- getGEO(filename = gzfile)
expr <- exprs(gse)
pheno <- pData(phenoData(gse))

cat(sprintf("Expression matrix: %d probes x %d samples\n", nrow(expr), ncol(expr)))
cat("Sample names:", paste(colnames(expr), collapse=", "), "\n")

# 保存
saveRDS(expr, file.path(GEO_DIR, "GSE83129_expression.rds"))
saveRDS(pheno, file.path(GEO_DIR, "GSE83129_pheno.rds"))

# 输出所有列名
cat("\nAll column names:\n")
cat(paste(colnames(pheno), collapse = "\n"), "\n")

# 检查所有有用的列
cat("\n=== Checking columns with limited unique values ===\n")
for (col in colnames(pheno)) {
  vals <- as.character(pheno[[col]])
  n_unique <- length(unique(vals))
  if (n_unique <= 40 && n_unique > 1) {
    cat(sprintf("\n--- %s (N=%d) ---\n", col, n_unique))
    print(table(vals, useNA = "ifany"))
  }
}

# 特别关注 characteristics 列
cat("\n=== Characteristics ===\n")
for (col in grep("characteristics", colnames(pheno), value = TRUE, ignore.case = TRUE)) {
  cat(sprintf("\n--- %s ---\n", col))
  print(head(as.character(pheno[[col]]), 40))
}

cat("\n=== Sample titles ===\n")
print(pheno$title)

cat("\n=== description ===\n")
for (col in grep("description", colnames(pheno), value = TRUE, ignore.case = TRUE)) {
  cat(sprintf("\n--- %s ---\n", col))
  print(head(as.character(pheno[[col]]), 40))
}

cat("\nDone.\n")
