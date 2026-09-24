# ============================================================
# 01_GEO_DATA_PIPELINE.R
# XELOX Resistance Study - GEO Data Download and Preprocessing
# ============================================================
# NOTE: Uses /path/to/xelox_project as project root (ASCII path
# to avoid Chinese character encoding issues in R on Windows)

# --- 0. Environment Setup ---
lib_path <- "/path/to/Rlibs"
.libPaths(c(lib_path, .libPaths()))
Sys.setenv(TMPDIR = "/tmp", TMP = "/tmp", TEMP = "/tmp")

library(GEOquery)
library(limma)
library(sva)
library(biomaRt)
library(WGCNA)
library(BiocManager)

# Project paths - use ASCII junction to avoid Chinese path issues
PROJECT_ROOT <- "/path/to/xelox_project"
DATA_GEO_DIR  <- file.path(PROJECT_ROOT, "data", "geo")
DATA_PROC_DIR <- file.path(PROJECT_ROOT, "data", "processed")
RESULTS_TAB_DIR <- file.path(PROJECT_ROOT, "results", "tables")

dir.create(DATA_GEO_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_PROC_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_TAB_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== XELOX Resistance Study: GEO Data Pipeline ===\n")
cat("Project root:", PROJECT_ROOT, "\n\n")

# ============================================================
# Part 1: Target GEO Datasets
# ============================================================
cat("[Step 1] Target GEO datasets for CRC + XELOX/Capecitabine/Oxaliplatin\n")

target_datasets <- c(
  "GSE39582",  # CRC adjuvant chemo (capecitabine+oxaliplatin subset available)
  "GSE14333",  # CRC prognosis with chemo info
  "GSE17538",  # CRC prognosis with adjuvant chemo
  "GSE72970",  # CRC capecitabine response
  "GSE104645", # CRC chemo response
  "GSE28702",  # CRC oxaliplatin sensitive/resistant
  "GSE42284",  # CRC 5-FU/oxaliplatin response
  "GSE19860",  # CRC capecitabine response
  "GSE69657",  # CRC XELOX neoadjuvant therapy
  "GSE144224"  # CRC chemo resistance
)

cat(paste("  Target datasets:", length(target_datasets), "\n"))
cat(paste("  ", paste(target_datasets, collapse = ", "), "\n\n"))

# ============================================================
# Part 2: Download Series Matrix Files
# ============================================================
cat("[Step 2] Downloading series matrix files\n")

download_geo_dataset <- function(gse_id, data_dir) {
  cat(paste("  Processing", gse_id, "...\n"))
  
  options("download.file.method.GEOquery" = "auto")
  options("GEOquery.inmemory.gpl" = FALSE)
  
  gse <- tryCatch({
    getGEO(gse_id, destdir = data_dir, GSEMatrix = TRUE, AnnotGPL = TRUE)
  }, error = function(e) {
    cat(paste("  ERROR for", gse_id, ":", conditionMessage(e), "\n"))
    return(NULL)
  })
  
  if (is.null(gse) || length(gse) == 0) {
    return(NULL)
  }
  
  return(gse[[1]])
}

geo_results <- list()
for (gse_id in target_datasets) {
  geo_obj <- download_geo_dataset(gse_id, DATA_GEO_DIR)
  if (!is.null(geo_obj)) {
    geo_results[[gse_id]] <- geo_obj
    cat(paste("  OK:", gse_id, "- Samples:", ncol(geo_obj), "\n"))
  } else {
    cat(paste("  FAILED:", gse_id, "\n"))
  }
}

cat(paste("\n  Successfully downloaded:", length(geo_results), "/", 
          length(target_datasets), "datasets\n"))

# ============================================================
# Part 3: Extract Clinical Information
# ============================================================
cat("\n[Step 3] Extracting clinical information\n")

extract_clinical_info <- function(gse_obj, gse_id) {
  pheno <- pData(gse_obj)
  col_names <- colnames(pheno)
  
  key_cols <- list(
    treatment = grep("treatment|therapy|chemotherap|regimen|drug|agent|chemo|adjuvant", 
                     col_names, ignore.case = TRUE, value = TRUE),
    survival_status = grep("dfs|pfs|survival|recurrence|progression|event|status|relapse", 
                           col_names, ignore.case = TRUE, value = TRUE),
    survival_time = grep("dfs_time|pfs_time|survival_time|time_to|months|days|follow_up", 
                         col_names, ignore.case = TRUE, value = TRUE),
    response = grep("response|responder|sensitive|resistant|resistance|complete_response", 
                    col_names, ignore.case = TRUE, value = TRUE),
    stage = grep("stage|tnm|tumor_node|dukes|classification", 
                 col_names, ignore.case = TRUE, value = TRUE),
    age = grep("age", col_names, ignore.case = TRUE, value = TRUE),
    sex = grep("sex|gender", col_names, ignore.case = TRUE, value = TRUE)
  )
  
  cat(paste("\n", gse_id, "clinical columns:\n"))
  for (cat_name in names(key_cols)) {
    if (length(key_cols[[cat_name]]) > 0) {
      cat(paste("   ", cat_name, ":", 
                paste(head(key_cols[[cat_name]], 5), collapse = ", "), "\n"))
    }
  }
  
  return(list(pheno = pheno, key_cols = key_cols))
}

clinical_info <- list()
for (gse_id in names(geo_results)) {
  clinical_info[[gse_id]] <- extract_clinical_info(geo_results[[gse_id]], gse_id)
}

# ============================================================
# Part 4: GSE39582 Deep Analysis - Primary Target
# ============================================================
cat("\n[Step 4] GSE39582 Deep Analysis - Extracting XELOX subset\n")

if ("GSE39582" %in% names(geo_results)) {
  gse39582 <- geo_results[["GSE39582"]]
  pheno39582 <- pData(gse39582)
  
  # Find chemo-related columns
  chemo_cols <- grep("chemotherapy|adjuvant|treatment|regimen|ct|chemo", 
                     colnames(pheno39582), ignore.case = TRUE, value = TRUE)
  cat(paste("  Chemo-related columns:", length(chemo_cols), "\n"))
  if (length(chemo_cols) > 0) {
    cat(paste("   ", paste(chemo_cols, collapse = ", "), "\n"))
  }
  
  # Show sample column names
  cat("\n  First 40 column names:\n")
  for (i in seq(1, min(40, length(colnames(pheno39582))), 5)) {
    cat(paste("    ", paste(colnames(pheno39582)[i:min(i+4, length(colnames(pheno39582)))], 
                           collapse = " | "), "\n"))
  }
  
  # Save full clinical data
  write.csv(pheno39582, 
            file = file.path(RESULTS_TAB_DIR, "GSE39582_clinical_data.csv"),
            row.names = TRUE)
  cat("\n  Clinical data saved to results/tables\n")
  
  # Examine characteristics/title columns for chemo info
  for (col in grep("characteristics_ch1|title", 
                   colnames(pheno39582), ignore.case = TRUE, value = TRUE)) {
    cat(paste("\n  Column:", col, "- first 8 values:\n"))
    print(head(as.character(pheno39582[[col]]), 8))
  }
  
  # Extract expression matrix
  expr_gse39582 <- exprs(gse39582)
  cat(paste("\n  Expression matrix:", nrow(expr_gse39582), "genes x", 
            ncol(expr_gse39582), "samples\n"))
  
  saveRDS(expr_gse39582, file = file.path(DATA_GEO_DIR, "GSE39582_expression.rds"))
  cat("  Expression matrix saved to data/geo\n")
  
  # ALSO save full ExpressionSet for clinical data extraction
  saveRDS(gse39582, file = file.path(DATA_GEO_DIR, "GSE39582_eset.rds"))
  cat("  Full ExpressionSet saved to data/geo (for clinical extraction)\n")
}

# ============================================================
# Part 5: Extract Expression Data for All Datasets
# ============================================================
cat("\n[Step 5] Extracting expression data for all datasets\n")

for (gse_id in names(geo_results)) {
  if (gse_id == "GSE39582") next
  
  tryCatch({
    expr_mat <- exprs(geo_results[[gse_id]])
    cat(paste(" ", gse_id, ":", nrow(expr_mat), "genes x", ncol(expr_mat), "samples\n"))
    
    saveRDS(expr_mat, file = file.path(DATA_GEO_DIR, paste0(gse_id, "_expression.rds")))
    
    pheno_df <- pData(geo_results[[gse_id]])
    write.csv(pheno_df, 
              file = file.path(RESULTS_TAB_DIR, paste0(gse_id, "_clinical_data.csv")),
              row.names = TRUE)
    
  }, error = function(e) {
    cat(paste("  ERROR extracting", gse_id, ":", conditionMessage(e), "\n"))
  })
}

# ============================================================
# Part 6: Summary Report
# ============================================================
cat("\n[Step 6] Generating data collection report\n")

if (length(geo_results) > 0) {
  report <- data.frame(
    Dataset = names(geo_results),
    Samples = sapply(geo_results, function(x) ncol(x)),
    Genes = sapply(geo_results, function(x) nrow(exprs(x))),
    Platform = sapply(geo_results, function(x) annotation(x)),
    stringsAsFactors = FALSE
  )
  
  cat("\n=== GEO Data Collection Summary ===\n")
  print(report)
  cat("===================================\n")
  
  write.csv(report, file = file.path(RESULTS_TAB_DIR, "geo_data_collection_report.csv"), 
            row.names = FALSE)
} else {
  cat("\n  WARNING: No datasets were downloaded successfully.\n")
  cat("  GEO FTP may be blocked by firewall or slow connection.\n")
  cat("  Consider downloading manually from NCBI GEO website.\n")
}

# Save workspace (using safe path)
save.image(file = file.path(DATA_PROC_DIR, "01_geo_pipeline_workspace.RData"))
cat("\nWorkspace saved\n")
cat("GEO data pipeline complete\n")
