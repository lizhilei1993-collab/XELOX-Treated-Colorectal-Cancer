# Try installing HDO.db with full error capture
Sys.setenv(HOME = "/tmp", TMPDIR = "/tmp", TEMP = "/tmp", TMP = "/tmp")
.libPaths(c("/path/to/Rlibs", .libPaths()))

# Capture ALL output including errors
con <- file("/tmp/hdo_install_log.txt", open = "wt")
sink(con, type = "output")
sink(con, type = "message")

result <- tryCatch({
  install.packages("/tmp/HDO.db_1.0.0.tar.gz", 
                   lib = "/path/to/Rlibs", 
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
cat("Full log saved to /tmp/hdo_install_log.txt\n")
