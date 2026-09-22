# ============================================================
# 00_RUN_ALL.R - XELOX Resistance Master Pipeline
# ============================================================
# Usage: Rscript scripts/00_run_all.R
#
# Pipeline steps:
#   1. Initialize environment (setup.R)
#   2. GEO data pipeline (01_geo_data_pipeline.R) - SKIP if already done
#   3. Clinical group extraction - FIXED v3 (02_extract_xelox_groups.R)
#   4. DEG analysis - limma (03_deg_analysis.R)
#   5. WGCNA analysis (04_wgcna_analysis.R)
#   6. Gene pool integration (05_gene_pool_integration.R)
#   7. Probe-to-gene mapping (06_probe_to_gene_mapping.R)
#   8. Enrichment analysis (07_enrichment_analysis.R) - GO/KEGG/DO/GSEA
#   9. TCGA data extraction (08_tcga_data_extraction.R) - STAR - Counts
#      → DESeq2 DEGs + log2(CPM+1) for ML/SHAP
# ============================================================

Sys.setenv(TMPDIR = "C:/temp", TMP = "C:/temp", TEMP = "C:/temp")
.libPaths(c("C:/Rlibs", .libPaths()))

PROJECT_ROOT <- "C:/xelox_research"
SCRIPT_DIR <- file.path(PROJECT_ROOT, "scripts")
DATA_GEO_DIR <- file.path(PROJECT_ROOT, "data", "geo")
DATA_TCGA_DIR <- file.path(PROJECT_ROOT, "data", "tcga")
DATA_PROC_DIR <- file.path(PROJECT_ROOT, "data", "processed")
RESULTS_TAB_DIR <- file.path(PROJECT_ROOT, "results", "tables")
RESULTS_FIG_DIR <- file.path(PROJECT_ROOT, "results", "figures")

# Create all directories
dir.create(DATA_GEO_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_TCGA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_PROC_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_FIG_DIR, showWarnings = FALSE, recursive = TRUE)

cat("===========================================================\n")
cat("XELOX Resistance Study - Master Analysis Pipeline\n")
cat("===========================================================\n\n")

# ============================================================
# Step 0: Check environment
# ============================================================
cat("=== Step 0: Environment Check ===\n\n")

required_pkgs <- c("limma", "DESeq2", "WGCNA", "GEOquery", "biomaRt", "survival", "sva")
for (pkg in required_pkgs) {
  status <- tryCatch({
    library(pkg, character.only = TRUE, quietly = TRUE)
    cat(sprintf("  ✓ %s v%s\n", pkg, packageVersion(pkg)))
  }, error = function(e) {
    cat(sprintf("  ✗ %s - %s\n", pkg, conditionMessage(e)))
    return(NULL)
  })
}
cat("\n")

# ============================================================
# Step 1: Check GEO data exists
# ============================================================
cat("=== Step 1: Checking GEO data ===\n\n")

# Check if GEO data is already downloaded - skip if eset exists
step1_done <- file.exists(file.path(DATA_GEO_DIR, "GSE39582_eset.rds"))
if (step1_done) {
  cat("  GEO data already exists. Skipping Step 1.\n")
  cat("  GSE39582_eset.rds found\n")
  
  # List all downloaded datasets
  expr_files <- list.files(DATA_GEO_DIR, pattern = "_expression\\.rds$")
  cat(sprintf("  Expression matrices found: %d\n", length(expr_files)))
} else {
  cat("  GEO data not found. Running 01_geo_data_pipeline.R...\n")
  source(file.path(SCRIPT_DIR, "01_geo_data_pipeline.R"))
}
cat("\n")

# ============================================================
# Step 2: Clinical group extraction (FIXED v3)
# ============================================================
cat("=== Step 2: Clinical Group Extraction (Fixed GSE39582) ===\n\n")

# Check if clinical extraction results already exist
step2_done <- file.exists(file.path(RESULTS_TAB_DIR, "GSE39582_xelox_groups.csv"))
if (step2_done) {
  cat("  Clinical extraction results found. Running fresh to ensure updated extraction...\n")
}

source(file.path(SCRIPT_DIR, "02_extract_xelox_groups.R"))
cat("\n")

# ============================================================
# Step 3: DEG Analysis - Limma
# ============================================================
cat("=== Step 3: DEG Differential Expression Analysis ===\n\n")

step3_done <- file.exists(file.path(RESULTS_TAB_DIR, "DEG_meta_analysis.csv"))
if (step3_done) {
  cat("  DEG results found. Re-running to incorporate fixed clinical data...\n")
}

source(file.path(SCRIPT_DIR, "03_deg_analysis.R"))
cat("\n")

# ============================================================
# Step 4: WGCNA Co-expression Network
# ============================================================
cat("=== Step 4: WGCNA Co-expression Network Analysis ===\n\n")

step4_done <- file.exists(file.path(RESULTS_TAB_DIR, "WGCNA_all_gene_module_membership.csv"))
if (step4_done) {
  cat("  WGCNA results found. Re-running with fixed clinical data...\n")
}

source(file.path(SCRIPT_DIR, "04_wgcna_analysis.R"))
cat("\n")

# ============================================================
# Step 5: Gene Pool Integration
# ============================================================
cat("=== Step 5: Gene Pool Integration ===\n\n")

source(file.path(SCRIPT_DIR, "05_gene_pool_integration.R"))
cat("\n")

# ============================================================
# Step 6: Probe-to-gene mapping (for DEG mapped files)
# ============================================================
cat("=== Step 6: Probe-to-Gene Mapping ===\n\n")

step6_done <- file.exists(file.path(RESULTS_TAB_DIR, "DEG_GSE39582_mapped.csv"))
if (!step6_done) {
  cat("  Running probe-to-gene mapping...\n")
  source(file.path(SCRIPT_DIR, "06_probe_to_gene_mapping.R"))
} else {
  cat("  DEG mapped files found. Skipping Step 6.\n")
}
cat("\n")

# ============================================================
# Step 7: Enrichment Analysis (GO/KEGG/DO + GSEA)
# ============================================================
cat("=== Step 7: Enrichment Analysis (GO + KEGG + DO + GSEA) ===\n\n")

source(file.path(SCRIPT_DIR, "07_enrichment_analysis.R"))
cat("\n")

# ============================================================
# Step 8: TCGA-COAD + TCGA-READ Validation Data
#         (STAR - Counts → DESeq2 + log2(CPM+1) for ML)
# ============================================================
cat("=== Step 8: TCGA Data Extraction (COAD + READ) ===\n\n")
cat("  Data source: GDC STAR - Counts (no FPKM available)\n")
cat("  Path A: DESeq2 differential expression (resistant vs sensitive)\n")
cat("  Path B: log2(CPM+1) → gene pool expression for ML/SHAP\n\n")

step8_done <- file.exists(file.path(DATA_TCGA_DIR, "TCGA_COAD_READ_expression_counts.rds"))
if (step8_done) {
  cat("  TCGA data found. Re-running to ensure updated extraction...\n")
} else {
  cat("  No TCGA cache found. Starting fresh download...\n")
}

source(file.path(SCRIPT_DIR, "08_tcga_data_extraction.R"))
cat("\n")

# ============================================================
# Step 9: ML + SHAP Modeling (XGBoost + LightGBM)
# ============================================================
cat("=== Step 9: ML + SHAP Molecular Fingerprint ===\n\n")
cat("  XGBoost + LightGBM + SHAP on GSE39582\n")
cat("  External validation on TCGA (if available)\n\n")

source(file.path(SCRIPT_DIR, "09_ml_shap_modeling.R"))
cat("\n")

# ============================================================
# Final Summary
# ============================================================
cat("===========================================================\n")
cat("PIPELINE COMPLETE\n")
cat("===========================================================\n\n")

cat("Output files:\n")
cat(sprintf("  Clinical groups:         %s\n", 
            file.path(RESULTS_TAB_DIR, "GSE39582_xelox_groups.csv")))
cat(sprintf("  DEG per dataset:         %s\n",
            file.path(RESULTS_TAB_DIR, "DEG_GSE39582_limma.csv")))
cat(sprintf("  DEG meta-analysis:       %s\n",
            file.path(RESULTS_TAB_DIR, "DEG_meta_analysis.csv")))
cat(sprintf("  WGCNA modules:           %s\n",
            file.path(RESULTS_TAB_DIR, "WGCNA_all_gene_module_membership.csv")))
cat(sprintf("  WGCNA hub genes:         %s\n",
            file.path(RESULTS_TAB_DIR, "WGCNA_hub_genes.csv")))
cat(sprintf("  Final gene pool:         %s\n",
            file.path(RESULTS_TAB_DIR, "XELOX_final_gene_pool.csv")))
cat(sprintf("  Integration report:      %s\n",
            file.path(RESULTS_TAB_DIR, "gene_pool_integration_report.txt")))
cat(sprintf("  Enrichment tables:       %s\n",
            file.path(RESULTS_TAB_DIR, "enrichment")))
cat(sprintf("  Enrichment figures:      %s\n",
            file.path(RESULTS_FIG_DIR)))
cat(sprintf("  GSEA figures:            %s\n",
            file.path(RESULTS_FIG_DIR, "gsea")))
cat(sprintf("  Enrichment summary:      %s\n",
            file.path(RESULTS_TAB_DIR, "enrichment", "enrichment_summary_report.txt")))
cat(sprintf("  TCGA DESeq2:             %s\n",
            file.path(RESULTS_TAB_DIR, "TCGA_DESeq2_results.csv")))
cat(sprintf("  TCGA counts (RDS):       %s\n",
            file.path(DATA_TCGA_DIR, "TCGA_COAD_READ_expression_counts.rds")))
cat(sprintf("  TCGA logCPM (RDS):       %s\n",
            file.path(DATA_TCGA_DIR, "TCGA_COAD_READ_expression_logcpm.rds")))
cat(sprintf("  TCGA gene pool (CSV):    %s\n",
            file.path(DATA_TCGA_DIR, "TCGA_COAD_READ_gene_pool_expression.csv")))
cat(sprintf("  ML performance:          %s\n",
            file.path(RESULTS_TAB_DIR, "ml", "ML_model_performance.csv")))
cat(sprintf("  SHAP summary:            %s\n",
            file.path(RESULTS_TAB_DIR, "ml", "SHAP_summary.csv")))
cat(sprintf("  Molecular fingerprint:   %s\n",
            file.path(RESULTS_TAB_DIR, "ml", "XELOX_molecular_fingerprint.csv")))
cat(sprintf("  ML report:               %s\n",
            file.path(RESULTS_TAB_DIR, "ml", "ML_SHAP_summary_report.txt")))

cat("\nNext steps:\n")
cat("  1. Review enrichment_summary_report.txt for GO/KEGG/DO/GSEA results\n")
cat("  2. Review TCGA_DESeq2_results.csv for DEGs (resistant vs sensitive)\n")
cat("  3. Review ML_SHAP_summary_report.txt for molecular fingerprint\n")
cat("  4. Review XELOX_molecular_fingerprint.csv for top predictive genes\n")
cat("  5. Consider wet-lab validation of top fingerprint genes\n")
