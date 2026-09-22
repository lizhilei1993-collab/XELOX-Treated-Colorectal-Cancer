# ============================================================
# 08_TCGA_DATA_EXTRACTION.R
# TCGA-COAD + TCGA-READ data extraction for XELOX resistance
# XELOX Resistance Study
# ============================================================
# CRITICAL: GDC TCGA data ONLY provides STAR - Counts (raw counts).
# FPKM/FPKM-UQ are NO LONGER available from GDC as of 2024-2025.
#
# DATA STRATEGY (single source, two analysis paths):
#   STAR - Counts (raw integer counts)
#      ├── [Path A] DESeq2 differential expression analysis
#      │     DESeq2's negative binomial model REQUIRES raw counts.
#      │     FPKM would be invalid input — Counts are the CORRECT type.
#      │
#      └── [Path B] log2(CPM+1) for ML / SHAP modeling
#            CPM (Counts Per Million) is library-size normalization.
#            Equivalent to FPKM without gene-length correction —
#            for between-sample comparisons within the same assay,
#            gene-length bias cancels out. No information lost.
#
# BENEFIT: No FPKM download needed → ~50% less download time.
# ============================================================

Sys.setenv(TMPDIR = "C:/temp", TMP = "C:/temp", TEMP = "C:/temp")
.libPaths(c("C:/Rlibs", .libPaths()))

PROJECT_ROOT     <- "C:/xelox_research"
DATA_TCGA_DIR    <- file.path(PROJECT_ROOT, "data", "tcga")
RESULTS_TAB_DIR  <- file.path(PROJECT_ROOT, "results", "tables")

dir.create(DATA_TCGA_DIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("TCGA-COAD + TCGA-READ Data Extraction\n")
cat("  STAR - Counts (RNA-Seq) | Clinical + Drug Data\n")
cat("  log2(CPM+1) for ML | XELOX/FOLFOX Patient Grouping\n")
cat("============================================================\n\n")

# ============================================================
# 0. Load packages
# ============================================================
cat("=== [0] Loading packages ===\n\n")

required_pkgs <- c("TCGAbiolinks", "SummarizedExperiment", "DESeq2")
for (pkg in required_pkgs) {
  tryCatch({
    library(pkg, character.only = TRUE, quietly = TRUE)
    cat(sprintf("  \u2713 %s v%s\n", pkg, packageVersion(pkg)))
  }, error = function(e) {
    cat(sprintf("  \u2717 %s - %s\n", pkg, conditionMessage(e)))
    stop(pkg, " is required")
  })
}
cat("\n")

# ============================================================
# Helper: download and prepare RNA-Seq data with checkpoint
# ============================================================
download_rnaseq <- function(project, data_dir, label) {
  rds_file <- file.path(data_dir, sprintf("%s_counts.rds", gsub("-", "_", project)))
  
  if (file.exists(rds_file)) {
    cat(sprintf("  [CACHED] Loading %s from %s\n", label, basename(rds_file)))
    return(readRDS(rds_file))
  }
  
  cat(sprintf("  [DOWNLOAD] Querying %s | STAR - Counts ...\n", project))
  
  query <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = "Primary Tumor"
  )
  
  n_samples <- nrow(getResults(query))
  cat(sprintf("  Found %d Primary Tumor samples\n", n_samples))
  
  if (n_samples == 0) {
    cat(sprintf("  WARNING: No samples for %s\n", project))
    return(NULL)
  }
  
  cat(sprintf("  Downloading %s (%d files)...\n", label, n_samples))
  tryCatch({
    GDCdownload(query, method = "api", files.per.chunk = 20)
  }, error = function(e) {
    cat(sprintf("  API download failed: %s\n", conditionMessage(e)))
    cat("  Retrying with 'client' method...\n")
    tryCatch({
      GDCdownload(query, method = "client", files.per.chunk = 20)
    }, error = function(e2) {
      cat(sprintf("  Client download also failed: %s\n", conditionMessage(e2)))
      return(NULL)
    })
  })
  
  cat(sprintf("  Preparing %s (this may take a while)...\n", label))
  se <- GDCprepare(query, save = FALSE)
  
  cat(sprintf("  Saving to %s ...\n", basename(rds_file)))
  saveRDS(se, rds_file)
  
  cat(sprintf("  Done: %d genes x %d samples\n", nrow(se), ncol(se)))
  return(se)
}

# ============================================================
# 1. Download RNA-Seq data (STAR - Counts)
# ============================================================
cat("=== [1] Downloading RNA-Seq STAR - Counts data ===\n\n")

se_coad <- download_rnaseq("TCGA-COAD", DATA_TCGA_DIR, "COAD Counts")
se_read <- download_rnaseq("TCGA-READ", DATA_TCGA_DIR, "READ Counts")

# ============================================================
# 2. Download clinical data (XML format with drug info)
# ============================================================
cat("\n=== [2] Downloading clinical + drug data ===\n\n")

download_clinical <- function(project, data_dir) {
  clin_file <- file.path(data_dir, sprintf("%s_clinical_xml.rds", gsub("-", "_", project)))
  
  if (file.exists(clin_file)) {
    cat(sprintf("  [CACHED] Loading %s from %s\n", project, basename(clin_file)))
    return(readRDS(clin_file))
  }
  
  cat(sprintf("  [DOWNLOAD] Querying clinical XML for %s ...\n", project))
  
  query_clin <- GDCquery(
    project = project,
    data.category = "Clinical",
    data.format = "bcr xml"
  )
  
  n_clin_samples <- nrow(getResults(query_clin))
  cat(sprintf("  Found %d clinical XML files\n", n_clin_samples))
  
  if (n_clin_samples == 0) {
    cat(sprintf("  WARNING: No clinical data for %s\n", project))
    return(NULL)
  }
  
  GDCdownload(query_clin, method = "api")
  
  cat("  Parsing patient data...\n")
  clin_patient <- GDCprepare_clinic(query_clin, clinical.info = "patient")
  cat(sprintf("  Patient table: %d rows x %d cols\n", nrow(clin_patient), ncol(clin_patient)))
  
  cat("  Parsing drug data...\n")
  clin_drug <- GDCprepare_clinic(query_clin, clinical.info = "drug")
  cat(sprintf("  Drug table: %d rows x %d cols\n", nrow(clin_drug), ncol(clin_drug)))
  
  result <- list(patient = clin_patient, drug = clin_drug)
  saveRDS(result, clin_file)
  
  return(result)
}

clin_coad <- download_clinical("TCGA-COAD", DATA_TCGA_DIR)
clin_read <- download_clinical("TCGA-READ", DATA_TCGA_DIR)

# ============================================================
# 3. Merge COAD + READ expression data
# ============================================================
cat("\n=== [3] Merging COAD + READ expression data ===\n\n")

merge_expression <- function(se_coad, se_read) {
  if (is.null(se_coad) || is.null(se_read)) {
    cat("  WARNING: Missing data for merge\n")
    return(NULL)
  }
  
  mat_coad <- assay(se_coad)
  mat_read <- assay(se_read)
  
  common_genes <- intersect(rownames(mat_coad), rownames(mat_read))
  cat(sprintf("  COAD=%d genes, READ=%d genes, Common=%d genes\n",
              nrow(mat_coad), nrow(mat_read), length(common_genes)))
  
  mat_coad <- mat_coad[common_genes, ]
  mat_read <- mat_read[common_genes, ]
  mat_read <- mat_read[rownames(mat_coad), ]
  
  merged <- cbind(mat_coad, mat_read)
  cat(sprintf("  Merged: %d genes x %d samples (%d COAD + %d READ)\n",
              nrow(merged), ncol(merged), ncol(mat_coad), ncol(mat_read)))
  
  origins <- data.frame(
    sample_barcode = colnames(merged),
    cancer_type = c(rep("COAD", ncol(mat_coad)), rep("READ", ncol(mat_read))),
    stringsAsFactors = FALSE
  )
  
  gene_info <- NULL
  if (!is.null(rowData(se_coad))) {
    gene_info <- as.data.frame(rowData(se_coad))
    gene_info <- gene_info[common_genes, , drop = FALSE]
  }
  
  return(list(counts = merged, origins = origins, gene_info = gene_info))
}

merged_data <- merge_expression(se_coad, se_read)

# ============================================================
# 4. Calculate log2(CPM+1) for ML (replaces FPKM)
# ============================================================
cat("\n=== [4] Calculating log2(CPM+1) normalization ===\n\n")
cat("  WHY CPM instead of FPKM:\n")
cat("    - GDC stopped serving FPKM from mid 2024\n")
cat("    - CPM = Counts Per Million (library-size normalization only)\n")
cat("    - For between-sample comparison within same assay type:\n")
cat("      gene-length bias cancels out → CPM ≈ FPKM in rank order\n")
cat("    - log2(CPM+1) is standard input for ML/SHAP (XGBoost, RF, etc.)\n")
cat("    - Saves ~50% download time vs downloading FPKM separately\n\n")

if (!is.null(merged_data)) {
  counts_mat <- merged_data$counts
  
  # Calculate CPM: counts / (library_size_in_millions)
  col_libsizes <- colSums(counts_mat) / 1e6
  cpm_mat <- sweep(counts_mat, 2, col_libsizes, "/")
  
  # log2(CPM + 1) transformation (+1 avoids log(0))
  logcpm_mat <- log2(cpm_mat + 1)
  
  cat(sprintf("  CPM range: [%.2f, %.2f]\n", min(cpm_mat), max(cpm_mat)))
  cat(sprintf("  log2(CPM+1) range: [%.2f, %.2f]\n", min(logcpm_mat), max(logcpm_mat)))
  
  merged_data$logcpm <- logcpm_mat
}

# ============================================================
# 5. DESeq2 differential expression (Path A: Counts → DESeq2)
# ============================================================
cat("\n=== [5] DESeq2 differential expression analysis ===\n\n")
cat("  Using raw STAR - Counts directly (DESeq2 requires integer counts)\n")
cat("  Comparing resistant vs sensitive among XELOX/FOLFOX patients\n\n")

run_deseq2 <- function(counts_mat, response_labels, label_name = "TCGA") {
  if (is.null(counts_mat) || is.null(response_labels) || length(response_labels) < 10) {
    cat("  SKIP: insufficient samples for DESeq2\n")
    return(NULL)
  }
  
  # Build colData
  coldata <- data.frame(
    condition = factor(response_labels, levels = c("sensitive", "resistant")),
    row.names = names(response_labels)
  )
  
  # Match samples
  common_samps <- intersect(names(response_labels), colnames(counts_mat))
  cat(sprintf("  Samples with response labels: %d\n", length(common_samps)))
  
  if (length(common_samps) < 10) {
    cat("  SKIP: too few labeled samples\n")
    return(NULL)
  }
  
  cnt_sub <- round(counts_mat[, common_samps, drop = FALSE])
  cd_sub  <- coldata[common_samps, , drop = FALSE]
  
  cat(sprintf("  DESeq2 input: %d genes x %d samples\n", nrow(cnt_sub), ncol(cnt_sub)))
  cat(sprintf("    Sensitive=%d, Resistant=%d\n",
              sum(cd_sub$condition == "sensitive"),
              sum(cd_sub$condition == "resistant")))
  
  # Run DESeq2
  dds <- DESeqDataSetFromMatrix(
    countData = cnt_sub,
    colData   = cd_sub,
    design    = ~ condition
  )
  
  # Pre-filter low-count genes (at least 10 counts in >= 5 samples)
  keep <- rowSums(counts(dds) >= 10) >= 5
  dds <- dds[keep, ]
  cat(sprintf("  Genes after low-count filter: %d\n", nrow(dds)))
  
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", "resistant", "sensitive"))
  
  cat(sprintf("  DEGs (p.adj < 0.05): %d\n", sum(res$padj < 0.05, na.rm = TRUE)))
  cat(sprintf("  Up in resistant: %d | Down in resistant: %d\n",
              sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE),
              sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)))
  
  # Return full results table
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- res_df[order(res_df$padj), ]
  return(res_df)
}

# ============================================================
# 6. Save merged expression data
# ============================================================
cat("\n=== [6] Saving merged expression data ===\n\n")

if (!is.null(merged_data)) {
  counts_file <- file.path(DATA_TCGA_DIR, "TCGA_COAD_READ_expression_counts.rds")
  saveRDS(merged_data$counts, counts_file)
  cat(sprintf("  Counts matrix: %s (%s)\n", basename(counts_file),
              format(object.size(merged_data$counts), units = "MB")))
  
  logcpm_file <- file.path(DATA_TCGA_DIR, "TCGA_COAD_READ_expression_logcpm.rds")
  saveRDS(merged_data$logcpm, logcpm_file)
  cat(sprintf("  log2(CPM+1) matrix: %s (%s)\n", basename(logcpm_file),
              format(object.size(merged_data$logcpm), units = "MB")))
  
  origins_file <- file.path(DATA_TCGA_DIR, "TCGA_COAD_READ_sample_origins.csv")
  write.csv(merged_data$origins, origins_file, row.names = FALSE)
  cat(sprintf("  Sample origins: %s\n", basename(origins_file)))
  
  geneinfo_file <- file.path(DATA_TCGA_DIR, "TCGA_COAD_READ_gene_info.csv")
  write.csv(merged_data$gene_info, geneinfo_file, row.names = FALSE)
  cat(sprintf("  Gene info: %s\n", basename(geneinfo_file)))
}

# ============================================================
# 7. Merge and process clinical data
# ============================================================
cat("\n=== [7] Merging clinical data ===\n\n")

merge_clinical <- function(clin_list, project_names) {
  all_patient <- list()
  all_drug <- list()
  
  for (i in seq_along(clin_list)) {
    if (is.null(clin_list[[i]])) next
    
    pt <- clin_list[[i]]$patient
    dr <- clin_list[[i]]$drug
    
    if (!is.null(pt) && nrow(pt) > 0) {
      pt$cancer_type <- project_names[i]
      all_patient[[project_names[i]]] <- pt
    }
    if (!is.null(dr) && nrow(dr) > 0) {
      dr$cancer_type <- project_names[i]
      all_drug[[project_names[i]]] <- dr
    }
  }
  
  result <- list()
  if (length(all_patient) > 0) {
    result$patient <- do.call(rbind, all_patient)
    cat(sprintf("  Patient table: %d rows\n", nrow(result$patient)))
  }
  if (length(all_drug) > 0) {
    result$drug <- do.call(rbind, all_drug)
    cat(sprintf("  Drug table: %d rows\n", nrow(result$drug)))
  }
  return(result)
}

clin_merged <- merge_clinical(list(clin_coad, clin_read), c("COAD", "READ"))

# Save clinical data
if (!is.null(clin_merged$patient)) {
  write.csv(clin_merged$patient,
            file.path(DATA_TCGA_DIR, "TCGA_COAD_READ_clinical_data.csv"),
            row.names = FALSE, na = "")
  cat("  Clinical data saved\n")
}

if (!is.null(clin_merged$drug)) {
  write.csv(clin_merged$drug,
            file.path(DATA_TCGA_DIR, "TCGA_COAD_READ_drug_data.csv"),
            row.names = FALSE, na = "")
  cat("  Drug data saved\n")
}

# ============================================================
# 8. Identify XELOX/FOLFOX patients
# ============================================================
cat("\n=== [8] Identifying XELOX/FOLFOX-treated patients ===\n\n")

if (is.null(clin_merged$drug) || nrow(clin_merged$drug) == 0) {
  cat("  ERROR: No drug data available.\n")
  analysis_set <- NULL
} else {
  drug_df <- clin_merged$drug
  cat(sprintf("  Drug columns: %s\n", paste(colnames(drug_df), collapse = ", ")))
  
  # Drug name column
  drug_name_col <- NULL
  for (col in c("drug_name", "pharmaceutical_therapy_drug_name", "drug", "treatment_drug_name")) {
    if (col %in% colnames(drug_df)) { drug_name_col <- col; break }
  }
  if (is.null(drug_name_col)) {
    char_cols <- names(which(sapply(drug_df, is.character)))
    drug_name_col <- char_cols[1]
    cat(sprintf("  Using column: %s\n", drug_name_col))
  } else {
    cat(sprintf("  Drug name column: %s\n", drug_name_col))
  }
  
  # --- Fuzzy matching ---
  fluoro_patt <- paste(c("fluorouracil", "5.fu", "5fu", "capecitabine", "xeloda",
                          "folfox", "xelox", "capox", "tegafur", "uft",
                          "s.1", "ts.1", "doxifluridine", "adrucil",
                          "efudex", "carac"), collapse = "|")
  oxali_patt  <- paste(c("oxaliplatin", "eloxatin", "l.ohp", "lohp",
                          "folfox", "xelox", "capox"), collapse = "|")
  
  drug_df$is_fluoro <- grepl(fluoro_patt, drug_df[[drug_name_col]], ignore.case = TRUE, perl = TRUE)
  drug_df$is_oxali  <- grepl(oxali_patt,  drug_df[[drug_name_col]], ignore.case = TRUE, perl = TRUE)
  
  pt_fluoro <- unique(drug_df$bcr_patient_barcode[drug_df$is_fluoro])
  pt_oxali  <- unique(drug_df$bcr_patient_barcode[drug_df$is_oxali])
  
  cat(sprintf("  Fluoropyrimidine-treated: %d patients\n", length(pt_fluoro)))
  cat(sprintf("  Oxaliplatin-treated:      %d patients\n", length(pt_oxali)))
  
  xelox_pts <- intersect(pt_fluoro, pt_oxali)
  cat(sprintf("  BOTH (XELOX/FOLFOX):      %d patients\n", length(xelox_pts)))
  
  if (length(xelox_pts) == 0) {
    cat("\n  WARNING: No XELOX/FOLFOX patients found.\n")
    cat("  Checking drug names (top 30):\n")
    drug_tab <- table(drug_df[[drug_name_col]])
    print(head(sort(drug_tab, decreasing = TRUE), 30))
    analysis_set <- NULL
  } else {
    cat("  Sample patients:\n")
    cat(sprintf("    %s\n", paste(head(xelox_pts, 5), collapse = "\n    ")))
    
    # ============================================================
    # 9. Define response groups
    # ============================================================
    cat("\n=== [9] Defining response groups ===\n\n")
    
   # Outcome column  
    outcome_col <- NULL
    for (col in c("drug_treatment_outcome", "measure_of_response",
                   "treatment_response", "response", "clinical_response")) {
      if (col %in% colnames(drug_df)) { outcome_col <- col; break }
    }
    
    if (!is.null(outcome_col)) {
      cat(sprintf("  Outcome column: %s\n", outcome_col))
      cat("  Outcome distribution:\n")
      print(table(drug_df[[outcome_col]], useNA = "always"))
      
      xelox_drugs <- drug_df[drug_df$bcr_patient_barcode %in% xelox_pts &
                               (drug_df$is_fluoro | drug_df$is_oxali), ]
      
      # Classify each patient
      response_map <- list()
      for (pt in xelox_pts) {
        outcomes <- xelox_drugs[xelox_drugs$bcr_patient_barcode == pt, outcome_col]
        outcomes <- outcomes[!is.na(outcomes) & outcomes != ""]
        if (length(outcomes) == 0) { response_map[[pt]] <- "unknown"; next }
        
        sens_patt <- paste(c("complete response", "partial response", "respond", 
                              "CR", "PR", "regression"), collapse = "|")
        res_patt  <- paste(c("progressive disease", "stable disease", "resist",
                              "PD", "SD", "progression", "recurrence"), collapse = "|")
        
        is_sens <- any(grepl(sens_patt, outcomes, ignore.case = TRUE))
        is_res  <- any(grepl(res_patt,  outcomes, ignore.case = TRUE))
        
        if (is_sens && !is_res)      response_map[[pt]] <- "sensitive"
        else if (is_res && !is_sens) response_map[[pt]] <- "resistant"
        else if (is_sens && is_res)  response_map[[pt]] <- "mixed"
        else                         response_map[[pt]] <- "unknown"
      }
      
      resp_df <- data.frame(patient_barcode = names(response_map),
                             response_group = unlist(response_map),
                             stringsAsFactors = FALSE)
      
      cat("\n  Response group distribution:\n")
      print(table(resp_df$response_group))
      
      analysis_set <- resp_df[resp_df$response_group %in% c("sensitive", "resistant"), ]
      cat(sprintf("\n  Analysis set: %d (%d sensitive, %d resistant)\n",
                  nrow(analysis_set),
                  sum(analysis_set$response_group == "sensitive"),
                  sum(analysis_set$response_group == "resistant")))
      
      write.csv(resp_df, file.path(DATA_TCGA_DIR, "TCGA_COAD_READ_xelox_response.csv"),
                row.names = FALSE)
      cat("  Response data saved\n")
    } else {
      cat("  WARNING: No outcome column found\n")
      cat(sprintf("  Available: %s\n", paste(colnames(drug_df), collapse = ", ")))
      analysis_set <- NULL
    }
    
    # Save XELOX patient details
    xelox_out <- drug_df[drug_df$bcr_patient_barcode %in% xelox_pts,
                         intersect(c("bcr_patient_barcode", drug_name_col,
                                      outcome_col, "is_fluoro", "is_oxali", "cancer_type"),
                                   colnames(drug_df))]
    write.csv(xelox_out, file.path(DATA_TCGA_DIR, "TCGA_COAD_READ_xelox_patients.csv"),
              row.names = FALSE)
    cat("  XELOX patient details saved\n")
  }
}

# ============================================================
# 10. Extract gene pool expression for analysis set
# ============================================================
cat("\n=== [10] Extracting gene pool expression ===\n\n")

gp_file <- file.path(RESULTS_TAB_DIR, "XELOX_final_gene_pool.csv")
if (!file.exists(gp_file)) {
  cat("  WARNING: Gene pool file not found.\n")
} else if (is.null(merged_data) || is.null(analysis_set) || nrow(analysis_set) == 0) {
  cat("  WARNING: Missing expression or analysis set.\n")
} else {
  gene_pool <- read.csv(gp_file, stringsAsFactors = FALSE)
  gp_genes <- toupper(trimws(gene_pool$Gene))
  cat(sprintf("  Gene pool: %d genes\n", length(gp_genes)))
  
  # Get gene symbols from rowData
  gi <- merged_data$gene_info
  if (!is.null(gi) && "gene_name" %in% colnames(gi)) {
    gene_sym <- toupper(trimws(gi$gene_name))
  } else if (!is.null(gi) && "external_gene_name" %in% colnames(gi)) {
    gene_sym <- toupper(trimws(gi$external_gene_name))
  } else {
    gene_sym <- toupper(trimws(rownames(merged_data$logcpm)))
  }
  
  matched <- which(gene_sym %in% gp_genes)
  cat(sprintf("  Found in logCPM matrix: %d / %d genes\n", length(matched), length(gp_genes)))
  
  if (length(matched) >= 10) {
    logcpm_mat <- merged_data$logcpm[matched, , drop = FALSE]
    rownames(logcpm_mat) <- gene_sym[matched]
    
    # Match patients
    col_pts <- sapply(strsplit(colnames(logcpm_mat), "-"), function(x) paste(x[1:3], collapse = "-"))
    label_map <- stats::setNames(analysis_set$response_group, analysis_set$patient_barcode)
    resp_labs <- label_map[col_pts]
    
    analysis_cols <- which(!is.na(resp_labs))
    cat(sprintf("  Samples with response labels: %d\n", length(analysis_cols)))
    
    if (length(analysis_cols) >= 10) {
      sub_mat <- logcpm_mat[, analysis_cols, drop = FALSE]
      sub_resp <- resp_labs[analysis_cols]
      
      output_df <- as.data.frame(t(sub_mat))
      output_df$patient_barcode <- col_pts[analysis_cols]
      output_df$sample_barcode  <- colnames(sub_mat)
      output_df$response <- sub_resp
      
      origin_map <- stats::setNames(merged_data$origins$cancer_type, merged_data$origins$sample_barcode)
      output_df$cancer_type <- origin_map[output_df$sample_barcode]
      
      write.csv(output_df, file.path(DATA_TCGA_DIR, "TCGA_COAD_READ_gene_pool_expression.csv"),
                row.names = FALSE)
      
      cat(sprintf("  Final dataset: %d samples, %d genes\n", nrow(output_df), nrow(sub_mat)))
      cat(sprintf("    Sensitive: %d | Resistant: %d\n",
                  sum(sub_resp == "sensitive"), sum(sub_resp == "resistant")))
      cat(sprintf("    COAD: %d | READ: %d\n",
                  sum(output_df$cancer_type == "COAD"), sum(output_df$cancer_type == "READ")))
    } else {
      cat("  WARNING: Too few analysis samples\n")
    }
  } else {
    cat("  WARNING: Too few genes matched\n")
  }
}

# ============================================================
# 11. DESeq2 re-run with actual response labels
# ============================================================
cat("\n=== [11] DESeq2 re-run with response labels ===\n\n")

if (!is.null(merged_data) && exists("analysis_set") && !is.null(analysis_set) && nrow(analysis_set) > 0) {
  col_pts <- sapply(strsplit(colnames(merged_data$counts), "-"),
                    function(x) paste(x[1:3], collapse = "-"))
  names(col_pts) <- colnames(merged_data$counts)
  
  pt2resp <- stats::setNames(analysis_set$response_group, analysis_set$patient_barcode)
  resp_for_counts <- pt2resp[col_pts]
  names(resp_for_counts) <- names(col_pts)
  resp_for_counts <- resp_for_counts[!is.na(resp_for_counts)]
  
  cat(sprintf("  Samples with response labels: %d\n", length(resp_for_counts)))
  cat(sprintf("    Sensitive=%d, Resistant=%d\n",
              sum(resp_for_counts == "sensitive"),
              sum(resp_for_counts == "resistant")))
  
  deseq2_results <- run_deseq2(merged_data$counts, resp_for_counts, "TCGA COAD+READ")
  
  if (!is.null(deseq2_results)) {
    deseq2_file <- file.path(RESULTS_TAB_DIR, "TCGA_DESeq2_results.csv")
    write.csv(deseq2_results, deseq2_file, row.names = FALSE)
    cat(sprintf("  DESeq2 results saved: %s\n", basename(deseq2_file)))
    
    # Also save top DEGs for easy viewing
    top_degs <- head(deseq2_results[deseq2_results$padj < 0.05, ], 100)
    if (nrow(top_degs) > 0) {
      top_file <- file.path(RESULTS_TAB_DIR, "TCGA_DESeq2_top100_DEGs.csv")
      write.csv(top_degs, top_file, row.names = FALSE)
      cat(sprintf("  Top 100 DEGs saved: %s\n", basename(top_file)))
    }
  }
} else {
  cat("  SKIP: insufficient data for DESeq2\n")
}

# ============================================================
# 12. Summary
# ============================================================
cat("\n=== [12] TCGA Extraction Summary ===\n\n")

cat("Output in", DATA_TCGA_DIR, ":\n")
for (f in list.files(DATA_TCGA_DIR)) {
  sz <- file.info(file.path(DATA_TCGA_DIR, f))$size
  cat(sprintf("  %-45s %s\n", f, format(sz, units = "auto")))
}

cat("\n=== TCGA Data Extraction Complete ===\n\n")
cat("Data strategy: STAR - Counts → [DESeq2] + [log2(CPM+1) for ML]\n\n")
cat("Output files:\n")
cat("  Expression:\n")
cat("    TCGA_COAD_READ_expression_counts.rds   - raw counts for DESeq2\n")
cat("    TCGA_COAD_READ_expression_logcpm.rds    - log2(CPM+1) for ML/SHAP\n")
cat("    TCGA_COAD_READ_gene_pool_expression.csv - gene pool subset for ML\n")
cat("  Clinical:\n")
cat("    TCGA_COAD_READ_clinical_data.csv        - patient clinical data\n")
cat("    TCGA_COAD_READ_drug_data.csv             - drug treatment records\n")
cat("  Response:\n")
cat("    TCGA_COAD_READ_xelox_response.csv       - patient response groups\n")
cat("    TCGA_COAD_READ_xelox_patients.csv       - XELOX patient details\n")
cat("  DESeq2:\n")
cat("    TCGA_DESeq2_results.csv                 - full DESeq2 results\n")
cat("    TCGA_DESeq2_top100_DEGs.csv             - top 100 DEGs\n")
cat("\nNext steps:\n")
cat("  1. Review TCGA_COAD_READ_xelox_response.csv for patient grouping\n")
cat("  2. Review TCGA_DESeq2_results.csv for DEGs (resistant vs sensitive)\n")
cat("  3. Use TCGA_COAD_READ_gene_pool_expression.csv for ML/SHAP\n")
cat("  4. If batch effect significant: sva::ComBat() correction\n")
cat("  5. Run GEO enrichment: Rscript scripts/07_enrichment_analysis.R\n")
