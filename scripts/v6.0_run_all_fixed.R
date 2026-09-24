# ============================================================
# V6.0 Run All Fixed Scripts - Wrapper with Chinese path fix
# ============================================================
# The key issue: R on Windows with C.UTF-8 locale cannot
# handle Chinese characters in file paths.
# Solution: Set locale to Chinese first, then source scripts.
# ============================================================

Sys.setenv(TMPDIR = "/tmp", TMP = "/tmp", TEMP = "/tmp")
.libPaths(c("/path/to/Rlibs", .libPaths()))

# Fix locale for Chinese paths
if (.Platform$OS.type == "windows") {
  tryCatch(Sys.setlocale("LC_ALL", "Chinese"), error = function(e) {
    tryCatch(Sys.setlocale("LC_ALL", "chs"), error = function(e2) {
      tryCatch(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.936"), error = function(e3) {})
    })
  })
}
cat(sprintf("Locale: %s\n", Sys.getlocale("LC_ALL")))

PROJECT_ROOT <- "/path/to/基于~1"
cat(sprintf("Project root: %s\n", PROJECT_ROOT))
cat(sprintf("Exists: %s\n", file.exists(PROJECT_ROOT)))

# Check key files
key_files <- c(
  "data/geo/GSE39582_eset.rds",
  "data/geo/GSE72970_expression.rds",
  "data/geo/GSE69657_expression.rds",
  "results/tables/GSE39582_xelox_groups.csv",
  "results/tables/GSE39582_clinical_data.csv",
  "results/tables/XELOX_final_gene_pool.csv"
)

cat("\nKey file check:\n")
all_ok <- TRUE
for (f in key_files) {
  full_path <- file.path(PROJECT_ROOT, f)
  exists <- file.exists(full_path)
  cat(sprintf("  %s: %s\n", ifelse(exists, "OK", "MISSING"), f))
  if (!exists) all_ok <- FALSE
}

if (!all_ok) {
  cat("\nERROR: Some key files are missing. Cannot proceed.\n")
  cat("Trying alternative: use list.files to check...\n")
  geo_dir <- file.path(PROJECT_ROOT, "data", "geo")
  if (file.exists(geo_dir)) {
    cat("  data/geo/ exists, contents:\n")
    cat(paste("  ", head(list.files(geo_dir), 10), collapse = "\n"), "\n")
  } else {
    cat("  data/geo/ does NOT exist\n")
    # Try listing what does exist
    if (file.exists(PROJECT_ROOT)) {
      cat("  Project root contents:\n")
      cat(paste("  ", list.files(PROJECT_ROOT), collapse = "\n"), "\n")
    } else {
      cat("  Project root does NOT exist either\n")
      # Try E: drive root
      cat("  /path/to/ contents:\n")
      cat(paste("  ", list.files("/path/to/"), collapse = "\n"), "\n")
    }
  }
  quit(status = 1)
}

cat("\nAll key files found. Ready to run scripts.\n")
cat("To run the full pipeline, uncomment the source() calls below.\n")
cat("NOTE: Each script takes 5-30 minutes depending on data size.\n")

# ============================================================
# Script 03: DEG Analysis (fixed filtering)
# ============================================================
cat("\n\n========== SCRIPT 03: DEG Analysis ==========\n")
cat("Starting DEG analysis with stricter filtering...\n")
source(file.path(PROJECT_ROOT, "scripts", "03_deg_analysis.R"))

# ============================================================
# Script 13: ComBat Batch Correction (fixed NA handling)
# ============================================================
cat("\n\n========== SCRIPT 13: ComBat Batch Correction ==========\n")
cat("Starting ComBat with NA removal...\n")
source(file.path(PROJECT_ROOT, "scripts", "13_combat_batch_correction.R"))

# ============================================================
# Script 10: Phase 2 ML + SHAP (fixed double-dipping)
# ============================================================
cat("\n\n========== SCRIPT 10: Phase 2 ML + SHAP ==========\n")
cat("Starting Phase 2 with LASSO CV pre-filter...\n")
source(file.path(PROJECT_ROOT, "scripts", "10_phase2_ml_shap.R"))

# ============================================================
# Script 30: GDSC Drug Sensitivity (fixed scoring)
# ============================================================
cat("\n\n========== SCRIPT 30: GDSC Drug Sensitivity ==========\n")
cat("Starting GDSC validation with weighted scoring...\n")
source(file.path(PROJECT_ROOT, "scripts", "30_oncoPredict_gdsc_validation.R"))

cat("\n\n========== ALL SCRIPTS COMPLETE ==========\n")
