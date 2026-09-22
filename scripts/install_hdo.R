# Try installing HDO.db with full error capture
Sys.setenv(HOME = "C:/temp", TMPDIR = "C:/temp", TEMP = "C:/temp", TMP = "C:/temp")
.libPaths(c("C:/Rlibs", .libPaths()))

# Capture ALL output including errors
con <- file("C:/temp/hdo_install_log.txt", open = "wt")
sink(con, type = "output")
sink(con, type = "message")

result <- tryCatch({
  install.packages("C:/temp/HDO.db_1.0.0.tar.gz", 
                   lib = "C:/Rlibs", 
                   repos = NULL, 
                   type = "source",
                   INSTALL_opts = "--no-byte-compile",
                   verbose = TRUE)
  "SUCCESS"
}, error = function(e) {
  paste("ERROR:", conditionMessage(e))
}, warning = function(w) {
  paste("WARNING:", conditionMessage(w))
})

sink(type = "message")
sink(type = "output")
close(con)

cat("Result:", result, "\n")
cat("Full log saved to C:/temp/hdo_install_log.txt\n")
