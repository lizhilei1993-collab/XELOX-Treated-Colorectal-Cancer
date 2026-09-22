.libPaths(c("C:/Rlibs", .libPaths()))
library(DOSE)

cat("=== enrichDO internal dispatch ===\n")
cat(deparse(body(enrichDO)), sep="\n")

cat("\n\n=== gseDisease function ===\n")
cat(deparse(body(gseDisease)), sep="\n")

cat("\n\n=== Check if HDO can be installed ===\n")
if (requireNamespace("BiocManager", quietly=TRUE)) {
  cat("BiocManager available\n")
  cat("Can install HDO.db:\n")
  cat(paste(tryCatch(BiocManager::available("HDO.db")[1], error=function(e) "NOT FOUND"), "\n"))
}

cat("\n\n=== Test enrichDO without ont parameter ===\n")
test_genes <- c("4312", "8318", "10874", "55143", "55388", "25805", "26060")
result <- tryCatch({
  enrichDO(gene = test_genes, pvalueCutoff = 0.05, readable = TRUE)
}, error = function(e) {
  cat("Error:", conditionMessage(e), "\n")
  NULL
})
cat("Result from enrichDO():", ifelse(is.null(result), "NULL", paste(nrow(result), "terms")), "\n")

cat("\n\n=== Test enrichDGN (DisGeNET, built-in) ===\n")
result_dgn <- tryCatch({
  enrichDGN(gene = test_genes, pvalueCutoff = 0.05, readable = TRUE)
}, error = function(e) {
  cat("Error:", conditionMessage(e), "\n")
  NULL
})
cat("Result from enrichDGN():", ifelse(is.null(result_dgn), "NULL", paste(nrow(result_dgn), "terms")), "\n")
if (!is.null(result_dgn)) {
  cat("Top 5 terms:\n")
  print(head(as.data.frame(result_dgn)[, c("Description", "p.adjust", "Count")], 5))
}

cat("\n\n=== Test gseDGN (DisGeNET, built-in) ===\n")
test_vec <- structure(seq(1, -1, length.out=100), names=as.character(1:100))
names(test_vec) <- test_genes
result_gse_dgn <- tryCatch({
  gseDGN(geneList = test_vec, pvalueCutoff = 0.05, seed=42)
}, error = function(e) {
  cat("Error:", conditionMessage(e), "\n")
  NULL
})
cat("Result from gseDGN():", ifelse(is.null(result_gse_dgn), "NULL", paste(nrow(result_gse_dgn), "terms")), "\n")
