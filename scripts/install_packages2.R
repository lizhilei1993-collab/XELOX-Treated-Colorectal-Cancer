# Install oncoPredict with temp dir fix
Sys.setenv(TMPDIR = "/tmp", TMP = "/tmp", TEMP = "/tmp")
.libPaths(c("/path/to/Rlibs", .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(Ncpus = 4)

# Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", lib = "/path/to/Rlibs")
}

# Try Bioconductor install
cat("Trying BiocManager::install(oncoPredict)...\n")
tryCatch({
    BiocManager::install("oncoPredict", lib = "/path/to/Rlibs", 
                         ask = FALSE, update = FALSE, force = TRUE)
}, error = function(e) {
    cat("BiocManager failed:", e$message, "\n")
    # Fallback: try CRAN install of pRRophetic
    cat("Trying pRRophetic from CRAN...\n")
    install.packages("pRRophetic", lib = "/path/to/Rlibs")
})

cat("Done.\n")
cat("oncoPredict:", require("oncoPredict"), "\n")
cat("pRRophetic:", require("pRRophetic"), "\n")
