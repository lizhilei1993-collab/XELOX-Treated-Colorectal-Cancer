# 28_download_gse83129.R
# 下载 GSE83129 (GPL6244, 26例 oxaliplatin, 9R/17S)
Sys.setenv(TMPDIR = "C:/temp", TMP = "C:/temp", TEMP = "C:/temp")
.libPaths(c("C:/Rlibs", .libPaths()))

PROJECT_ROOT <- "C:/xelox_research"
GEO_DIR <- file.path(PROJECT_ROOT, "data", "geo")

library(GEOquery)

# 下载 GSE83129
gse <- getGEO("GSE83129", destdir = GEO_DIR, GSEMatrix = TRUE, AnnotGPL = TRUE)

# 检查平台
cat("Platforms:\n")
print(sapply(gse, annotation))

# 提取表达矩阵
gse_gpl6244 <- gse[[1]]  # GPL6244
expr <- exprs(gse_gpl6244)
pheno <- pData(phenoData(gse_gpl6244))

cat(sprintf("\nExpression matrix: %d genes x %d samples\n", nrow(expr), ncol(expr)))
cat("Sample names:", colnames(expr), "\n")

# 保存
saveRDS(expr, file.path(GEO_DIR, "GSE83129_expression.rds"))
saveRDS(pheno, file.path(GEO_DIR, "GSE83129_pheno.rds"))

# 输出临床信息
cat("\n--- Pheno columns ---\n")
print(head(pheno[, 1:min(10, ncol(pheno))], 2))
cat("\nAll column names:\n")
cat(paste(colnames(pheno), collapse = "\n"), "\n")

# 检查是否有 response/treatment 信息
for (col in colnames(pheno)) {
  vals <- as.character(pheno[[col]])
  if (length(unique(vals)) < 20 && length(unique(vals)) > 1) {
    cat(sprintf("\n--- %s ---\n", col))
    print(table(vals, useNA = "ifany"))
  }
}

cat("\nDone.\n")
