# Test fixing locale and installing HDO.db
Sys.setenv(HOME = "C:/temp")
Sys.setenv(TMPDIR = "C:/temp")

# Fix locale - try several options
locales <- c(
  tryCatch(Sys.setlocale("LC_ALL", "English_United States.1252"), error=function(e) NULL),
  tryCatch(Sys.setlocale("LC_ALL", "en_US.UTF-8"), error=function(e) NULL),
  tryCatch(Sys.setlocale("LC_ALL", "C"), error=function(e) NULL)
)

cat("Locale:", Sys.getlocale(), "\n")
cat("HOME:", Sys.getenv("HOME"), "\n")
cat("TMPDIR:", Sys.getenv("TMPDIR"), "\n")
cat("tempdir():", tempdir(), "\n")

.libPaths(c("C:/Rlibs", .libPaths()))

# Try to install HDO.db
cat("\n=== Installing HDO.db ===\n")
if (requireNamespace("BiocManager", quietly=TRUE)) {
  result <- tryCatch({
    BiocManager::install("HDO.db", lib="C:/Rlibs", ask=FALSE, update=FALSE)
    "SUCCESS"
  }, error = function(e) {
    paste("FAILED:", conditionMessage(e))
  })
  cat("Install result:", result, "\n")
}

cat("\n=== Check if HDO.db is now available ===\n")
cat("HDO.db available:", require("HDO.db", quietly=TRUE), "\n")
if (require("HDO.db", quietly=TRUE)) {
  cat("Version:", as.character(packageVersion("HDO.db")), "\n")
}
