# Install required R packages for drug sensitivity prediction and scRNA-seq
.libPaths(c("/path/to/Rlibs", .libPaths()))

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(Ncpus = 4)

# First, install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", lib = "/path/to/Rlibs")
}

# Install oncoPredict from Bioconductor
cat("Installing oncoPredict...\n")
BiocManager::install("oncoPredict", lib = "/path/to/Rlibs", ask = FALSE, update = FALSE, force = TRUE)
cat("oncoPredict done.\n")

# Verify
cat("\nVerification:\n")
cat("oncoPredict:", require("oncoPredict"), "\n")

cat("\nAll done.\n")
