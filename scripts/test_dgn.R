# Test DGN (DisGeNET) enrichment - uses built-in DOSE data
Sys.setenv(HOME = "/tmp", TMPDIR = "/tmp", TEMP = "/tmp", TMP = "/tmp")
.libPaths(c("/path/to/Rlibs", .libPaths()))

library(DOSE)
library(org.Hs.eg.db)

cat("DOSE version:", as.character(packageVersion("DOSE")), "\n\n")

# Test with a few cancer-related genes
test_genes <- c("4312", "8318", "10874", "55143", "55388", "25805", "26060", "7157", "1956")

cat("=== enrichDGN (DisGeNET) ===\n")
result_dgn <- tryCatch({
  enrichDGN(gene = test_genes, pvalueCutoff = 1, readable = TRUE)
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  NULL
})
if (!is.null(result_dgn)) {
  cat("Terms found:", nrow(result_dgn), "\n")
  print(head(as.data.frame(result_dgn)[, c("Description", "p.adjust", "Count")], 5))
} else {
  cat("FAILED\n")
}

cat("\n=== enrichNCG (Network of Cancer Genes) ===\n")
result_ncg <- tryCatch({
  enrichNCG(gene = test_genes, pvalueCutoff = 1, readable = TRUE)
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  NULL
})
if (!is.null(result_ncg)) {
  cat("Terms found:", nrow(result_ncg), "\n")
  print(head(as.data.frame(result_ncg)[, c("Description", "p.adjust", "Count")], 5))
} else {
  cat("FAILED\n")
}

cat("\n=== enrichDGNv (Variant DisGeNET) ===\n")
result_dgnv <- tryCatch({
  enrichDGNv(gene = test_genes, pvalueCutoff = 1, readable = TRUE)
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  NULL
})
if (!is.null(result_dgnv)) {
  cat("Terms found:", nrow(result_dgnv), "\n")
  print(head(as.data.frame(result_dgnv)[, c("Description", "p.adjust", "Count")], 5))
} else {
  cat("FAILED\n")
}

cat("\n=== Package search: HDO.db ===\n")
cat("HDO.db in /path/to/Rlibs:", file.exists("/path/to/Rlibs/HDO.db"), "\n")
cat("HDO.db in system lib:", file.exists(file.path(.Library, "HDO.db")), "\n")
