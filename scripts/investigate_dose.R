# Investigate DOSE v4.5.0 internal objects
.libPaths(c("/path/to/Rlibs", .libPaths()))
library(DOSE)
library(org.Hs.eg.db)
library(GO.db)

cat("=== DOSE version ===\n")
cat(as.character(packageVersion("DOSE")), "\n\n")

cat("=== enrichDO args ===\n")
print(args(enrichDO))

cat("\n=== gseDO args ===\n")
print(args(gseDO))

cat("\n=== gseDO body ===\n")
cat(deparse(body(gseDO)), sep="\n")

cat("\n\n=== Search for TERMID2EXTID / PATHID2EXTID in DOSE ===\n")
# Search all functions in DOSE
ns <- getNamespace("DOSE")
funcs <- ls(ns, all.names=TRUE)

for(f in funcs) {
  if (exists(f, ns) && is.function(get(f, ns))) {
    body_text <- deparse(body(get(f, ns)))
    if (any(grepl("TERMID2EXTID|PATHID2EXTID", body_text))) {
      cat("Found in function:", f, "\n")
      lines <- grep("TERMID2EXTID|PATHID2EXTID", body_text, value=TRUE)
      cat(paste("  ", lines), sep="\n")
    }
  }
}

cat("\n=== Available data objects in DOSE ===\n")
cat(paste(data(package="DOSE")$results[, "Item"], collapse="\n"), "\n")

cat("\n=== Check for HDO.db package ===\n")
cat("HDO.db available:", require("HDO.db", quietly=TRUE), "\n")
if (require("HDO.db", quietly=TRUE)) {
  cat("HDO.db version:", as.character(packageVersion("HDO.db")), "\n")
  cat("HDO.db location:", find.package("HDO.db"), "\n")
}
