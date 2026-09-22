# Bioconductor + CRAN packages the reproduction needs.
#
# Why this file exists: at the time of deposit, none of limma, edgeR, Biobase or
# GSVA was present in the R library on the analysis machine (see
# verify/env_check_at_deposit.txt), even though the manuscripts results depend on
# them. The DEG layer therefore cannot be re-run straight after cloning.
#
# Run once, ideally into a project-local library so the shared R install is not
# touched:
#   Rscript install_deps.R
#
# Windows gets precompiled Bioc binaries, so no toolchain is required for these
# particular packages.

local <- file.path(getwd(), "r_libs_release")
dir.create(local, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(local, .libPaths()))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = local, repos = "https://cloud.r-project.org")
}

BiocManager::install(c("limma", "edgeR", "Biobase", "GSVA", "GSEABase"),
                     lib = local, update = FALSE, ask = FALSE)

install.packages(c("data.table", "clusterProfiler", "glmnet", "survival", "rms",
                   "WGCNA", "randomForest", "xgboost", "boot"),
                 lib = local, repos = "https://cloud.r-project.org")

need <- c("limma", "edgeR", "Biobase", "GSVA", "GSEABase", "data.table",
          "clusterProfiler", "glmnet", "survival", "rms", "WGCNA", "boot")
ok <- vapply(need, function(p) requireNamespace(p, quietly = TRUE), logical(1))
report <- file.path(getwd(), "verify", "session_after_install.txt")
dir.create(dirname(report), showWarnings = FALSE, recursive = TRUE)
writeLines(c(R.version.string,
             paste0("packages available: ", paste(names(ok)[ok], collapse = ", ")),
             paste0("MISSING: ", if (any(!ok)) paste(names(ok)[!ok], collapse = ", ") else "none"),
             "",
             capture.output(sessionInfo())),
           con = report)
cat("installed into", local, "\n")
cat("missing:", if (any(!ok)) paste(names(ok)[!ok], collapse = ", ") else "none", "\n")
