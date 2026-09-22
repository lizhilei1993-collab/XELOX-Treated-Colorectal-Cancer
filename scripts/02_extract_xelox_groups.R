# ============================================================
# 02_EXTRACT_XELOX_GROUPS.R (v3)
# Extract XELOX-treated patients and define resistance/sensitive groups
# ============================================================

lib_path <- "C:/Rlibs"
.libPaths(c(lib_path, .libPaths()))
Sys.setenv(TMPDIR = "C:/temp", TMP = "C:/temp", TEMP = "C:/temp")

library(GEOquery)
library(survival)

PROJECT_ROOT <- "C:/xelox_research"
DATA_GEO_DIR  <- file.path(PROJECT_ROOT, "data", "geo")
DATA_PROC_DIR <- file.path(PROJECT_ROOT, "data", "processed")
RESULTS_TAB_DIR <- file.path(PROJECT_ROOT, "results", "tables")
RESULTS_FIG_DIR <- file.path(PROJECT_ROOT, "results", "figures")

dir.create(DATA_PROC_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_TAB_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== XELOX Resistance Study: Clinical Group Extraction v3 ===\n\n")

# ============================================================
# 1. GSE39582 - Primary CIT dataset (n=585)
# ============================================================
cat("[1] Processing GSE39582 (CIT dataset, n=585)\n")

gse39582 <- readRDS(file.path(DATA_GEO_DIR, "GSE39582_eset.rds"))
pheno39582 <- pData(gse39582)

cat("  Columns in pheno data:", ncol(pheno39582), "\n")

# Use :ch1 columns directly (pre-parsed values)
clin39582 <- data.frame(
  sample_id = rownames(pheno39582),
  title = pheno39582$title,
  dataset = pheno39582[["dataset:ch1"]],
  
  adjuvant_chemo = pheno39582[["chemotherapy.adjuvant:ch1"]],
  chemo_type = pheno39582[["chemotherapy.adjuvant.type:ch1"]],
  
  rfs_event = as.numeric(pheno39582[["rfs.event:ch1"]]),
  rfs_delay = as.numeric(pheno39582[["rfs.delay:ch1"]]),
  os_event = as.numeric(pheno39582[["os.event:ch1"]]),
  os_delay = as.numeric(pheno39582[["os.delay (months):ch1"]]),
  
  sex = pheno39582[["Sex:ch1"]],
  age = as.numeric(pheno39582[["age.at.diagnosis (year):ch1"]]),
  tnm_stage = pheno39582[["tnm.stage:ch1"]],
  tumor_location = pheno39582[["tumor.location:ch1"]],
  mmr_status = pheno39582[["mmr.status:ch1"]],
  kras = pheno39582[["kras.mutation:ch1"]],
  braf = pheno39582[["braf.mutation:ch1"]],
  cit_subtype = pheno39582[["cit.molecularsubtype:ch1"]],
  
  stringsAsFactors = FALSE
)

cat("  Samples extracted:", nrow(clin39582), "\n")
cat("  RFS delay range:", range(clin39582$rfs_delay, na.rm = TRUE), "\n")
cat("  OS delay range:", range(clin39582$os_delay, na.rm = TRUE), "\n")
cat("  OS delay available:", sum(!is.na(clin39582$os_delay)), "\n")
cat("  Age range:", range(clin39582$age, na.rm = TRUE), "\n")
cat("  Adjuvant chemo (Y):", sum(clin39582$adjuvant_chemo == "Y", na.rm = TRUE), "\n")

cat("  Chemo type distribution:\n")
print(table(clin39582$chemo_type, useNA = "always"))

# XELOX-like keyword matching
xelox_keywords <- "XELOX|CAPOX|Capecitabine|capecitabin|Xeloda|FOLFOX|5-FU|fluoropyrimidine|oxaliplatin|L-OHP"
xelox_patients <- grepl(xelox_keywords, clin39582$chemo_type, ignore.case = TRUE)
cat("\n  XELOX-like regimens:", sum(xelox_patients, na.rm = TRUE), "\n")

xelox_subset <- clin39582[xelox_patients & !is.na(clin39582$chemo_type), ]
cat("  Breakdown:\n")
print(table(xelox_subset$chemo_type))

# --- Resistance definition ---
cat("\n  Defining resistance/sensitive groups...\n")
chemo_idx <- which(clin39582$adjuvant_chemo == "Y" & !is.na(clin39582$adjuvant_chemo) &
                     !is.na(clin39582$rfs_event))
cat("  Chemo-treated with RFS data:", length(chemo_idx), "\n")

if (length(chemo_idx) > 0) {
  chemo_patients <- clin39582[chemo_idx, ]
  
  resistant_idx <- which(chemo_patients$rfs_event == 1 & chemo_patients$rfs_delay <= 12)
  sensitive_idx <- which(chemo_patients$rfs_event == 0 & chemo_patients$rfs_delay >= 36)
  
  chemo_patients$group <- "intermediate"
  chemo_patients$group[resistant_idx] <- "resistant"
  chemo_patients$group[sensitive_idx] <- "sensitive"
  
  cat("  Resistant (RFS<=12mo):", length(resistant_idx), "\n")
  cat("  Sensitive (FU>=36mo):", length(sensitive_idx), "\n")
  cat("  Intermediate:", nrow(chemo_patients) - length(resistant_idx) - length(sensitive_idx), "\n")
  
  xelox_chemo <- chemo_patients[grepl(xelox_keywords, chemo_patients$chemo_type, ignore.case = TRUE), ]
  cat("\n  XELOX-like Subgroup:\n")
  cat("  Total:", nrow(xelox_chemo),
      "| Resistant:", sum(xelox_chemo$group == "resistant"),
      "| Sensitive:", sum(xelox_chemo$group == "sensitive"),
      "| Intermediate:", sum(xelox_chemo$group == "intermediate"), "\n")
  
  write.csv(chemo_patients, file.path(RESULTS_TAB_DIR, "GSE39582_xelox_groups.csv"), row.names = FALSE)
}

write.csv(clin39582, file.path(RESULTS_TAB_DIR, "GSE39582_clinical_processed.csv"), row.names = FALSE)
saveRDS(clin39582, file.path(DATA_PROC_DIR, "GSE39582_clinical_processed.rds"))
cat("  Saved clinical data\n")

# ============================================================
# 2. GSE72970 - CRC Capecitabine Response (n=124)
# ============================================================
cat("\n[2] Processing GSE72970 (CRC capecitabine response, n=124)\n")

pheno72970 <- read.csv(file.path(RESULTS_TAB_DIR, "GSE72970_clinical_data.csv"),
                        row.names = 1, check.names = FALSE)

regimen_col <- "regimen:ch1"
response_cat <- "response category:ch1"

cat("  Regimen distribution:\n")
print(table(pheno72970[[regimen_col]]))
cat("  Response category:\n")
print(table(pheno72970[[response_cat]]))

# Group XELOX-related regimens
xelox_regimens <- c("XELOX", "XELIRI", "CAPOX", "FOLFOX")
xelox_flag <- grepl(paste(xelox_regimens, collapse = "|"), pheno72970[[regimen_col]], ignore.case = TRUE)

clin72970 <- data.frame(
  sample_id = rownames(pheno72970),
  regimen = pheno72970[[regimen_col]],
  response = pheno72970[[response_cat]],
  is_xelox_like = xelox_flag,
  stringsAsFactors = FALSE
)

n_r <- sum(xelox_flag & pheno72970[[response_cat]] %in% c("PD", "SD"), na.rm = TRUE)
n_s <- sum(xelox_flag & pheno72970[[response_cat]] %in% c("CR", "PR"), na.rm = TRUE)
cat("  XELOX-like (PD/SD=Resistant):", n_r, "| (CR/PR=Sensitive):", n_s, "\n")

write.csv(clin72970, file.path(RESULTS_TAB_DIR, "GSE72970_clinical_summary.csv"), row.names = FALSE)

# ============================================================
# 3. GSE104645 - CRC 1st-line Chemotherapy (n=193)
# ============================================================
cat("\n[3] Processing GSE104645 (CRC 1st-line chemo response, n=193)\n")

pheno104645 <- read.csv(file.path(RESULTS_TAB_DIR, "GSE104645_clinical_data.csv"),
                         row.names = 1, check.names = FALSE)

regimen_col <- "1st-line chemotherapy regimens:ch1"
response_col <- "best response of 1st-line chemotherapy:ch1"

cat("  Top regimens:\n")
top_reg <- head(sort(table(pheno104645[[regimen_col]]), decreasing = TRUE), 8)
print(top_reg)

cat("  Response:\n")
print(table(pheno104645[[response_col]]))

# XELOX/FOLFOX grouping
xelox_flag <- grepl("XELOX|FOLFOX|SOX|CAPOX", pheno104645[[regimen_col]], ignore.case = TRUE)
n_xr <- sum(xelox_flag & pheno104645[[response_col]] %in% c("Progressive disease", "Stable disease"), na.rm = TRUE)
n_xs <- sum(xelox_flag & pheno104645[[response_col]] %in% c("Complete response", "Partial response"), na.rm = TRUE)
cat("  XELOX/FOLFOX: Resistant (PD/SD):", n_xr, "| Sensitive (CR/PR):", n_xs, "\n")

# ============================================================
# 4. GSE28702 - CRC Oxaliplatin (MFOLFOX6) (n=83)
# ============================================================
cat("\n[4] Processing GSE28702 (CRC mFOLFOX6 response, n=83)\n")

pheno28702 <- read.csv(file.path(RESULTS_TAB_DIR, "GSE28702_clinical_data.csv"),
                        row.names = 1, check.names = FALSE)

cat("  mfolfox6 response:\n")
print(table(pheno28702[["mfolfox6:ch1"]]))

n_r <- sum(pheno28702[["mfolfox6:ch1"]] == "non-responder", na.rm = TRUE)
n_s <- sum(pheno28702[["mfolfox6:ch1"]] == "responder", na.rm = TRUE)
cat("  non-responder:", n_r, "| responder:", n_s, "\n")

# ============================================================
# 5. GSE19860 - CRC Capecitabine Response (n=40)
# ============================================================
cat("\n[5] Processing GSE19860 (CRC capecitabine response, n=40)\n")

pheno19860 <- read.csv(file.path(RESULTS_TAB_DIR, "GSE19860_clinical_data.csv"),
                        row.names = 1, check.names = FALSE)

cat("  Treatment response:\n")
print(table(pheno19860[["treatment response:ch1"]]))

# FL = Fluoropyrimidine (capecitabine/5-FU)
fl_responder <- grepl("FL_Responder", pheno19860[["treatment response:ch1"]])
fl_nonresponder <- grepl("FL_Non_responder", pheno19860[["treatment response:ch1"]])
cat("  FL (capecitabine) Responder:", sum(fl_responder), "| Non-responder:", sum(fl_nonresponder), "\n")

# ============================================================
# 6. GSE69657 - CRC XELOX Neoadjuvant (n=30) [PURE XELOX]
# ============================================================
cat("\n[6] Processing GSE69657 (CRC XELOX neoadjuvant, n=30)\n")

pheno69657 <- read.csv(file.path(RESULTS_TAB_DIR, "GSE69657_clinical_data.csv"),
                        row.names = 1, check.names = FALSE)

cat("  Chemoresponse:\n")
print(table(pheno69657[["chemoresponse:ch1"]]))

n_r <- sum(pheno69657[["chemoresponse:ch1"]] == "noresponder", na.rm = TRUE)
n_s <- sum(pheno69657[["chemoresponse:ch1"]] == "responder", na.rm = TRUE)
cat("  Resistant (NR/noresponder):", n_r, "| Sensitive (R/responder):", n_s, "\n")

# ============================================================
# 7. GSE42284 - CRC Molecular Characterization (n=188)
# ============================================================
cat("\n[7] Processing GSE42284 (CRC molecular, n=188 - no response data)\n")

pheno42284 <- read.csv(file.path(RESULTS_TAB_DIR, "GSE42284_clinical_data.csv"),
                        row.names = 1, check.names = FALSE)
cat("  Molecular data available: stage, KRAS, BRAF, MSI, PIK3CA\n")

# ============================================================
# 8-10. Other datasets
# ============================================================
cat("\n[8] GSE144224 (n=3 - too small, SKIP)\n")
cat("[9] GSE14333 expression dims:", paste(dim(readRDS(file.path(DATA_GEO_DIR, "GSE14333_expression.rds"))), collapse=" x "), "\n")
cat("[10] GSE17538 expression dims:", paste(dim(readRDS(file.path(DATA_GEO_DIR, "GSE17538_expression.rds"))), collapse=" x "), "\n")

# ============================================================
# CONSOLIDATED SUMMARY
# ============================================================
cat("\n"); cat(paste(rep("=", 60), collapse="")); cat("\n")
cat("CONSOLIDATED CLINICAL GROUP SUMMARY\n")
cat(paste(rep("=", 60), collapse="")); cat("\n\n")

# Build summary list
summaries <- list(
  data.frame(
    Dataset = "GSE39582",
    n = 585,
    Platform = "Array",
    Cohort = "CIT (adjuvant FOLFOX/5-FU)",
    Resistant_def = "RFS event <= 12mo",
    Sensitive_def = "No event, FU >= 36mo",
    n_Resistant = 45,
    n_Sensitive = 119,
    n_Intermediate = 75,
    XELOX_specific = "FOLFOX (n=23): 7R/1S/15I",
    Notes = "Primary CIT cohort. All chemo types pooled for main analysis."
  ),
  data.frame(
    Dataset = "GSE72970",
    n = 124,
    Platform = "Array",
    Cohort = "Metastatic CRC chemo response",
    Resistant_def = "PD/SD",
    Sensitive_def = "CR/PR",
    n_Resistant = sum(xelox_flag & pheno72970[[response_cat]] %in% c("PD", "SD"), na.rm = TRUE),
    n_Sensitive = sum(xelox_flag & pheno72970[[response_cat]] %in% c("CR", "PR"), na.rm = TRUE),
    n_Intermediate = NA,
    XELOX_specific = "FOLFOX/XELIRI (n=37): response-based",
    Notes = "Multiple regimens (FOLFOX/FOLFIRI/XELIRI). Subset to FOLFOX-like."
  ),
  data.frame(
    Dataset = "GSE104645",
    n = 193,
    Platform = "Array",
    Cohort = "1st-line metastatic CRC",
    Resistant_def = "PD/SD",
    Sensitive_def = "CR/PR",
    n_Resistant = n_xr,
    n_Sensitive = n_xs,
    n_Intermediate = NA,
    XELOX_specific = "FOLFOX/SOX/XELOX (n=117): response-based",
    Notes = "Largest XELOX/FOLFOX response dataset."
  ),
  data.frame(
    Dataset = "GSE28702",
    n = 83,
    Platform = "Array (GPL570)",
    Cohort = "mFOLFOX6 responder vs non-responder",
    Resistant_def = "non-responder (label)",
    Sensitive_def = "responder (label)",
    n_Resistant = n_r,
    n_Sensitive = n_s,
    n_Intermediate = NA,
    XELOX_specific = "mFOLFOX6 (n=83): clean binary labels",
    Notes = "Pre-treatment biopsies. Clean responder/non-responder labels."
  ),
  data.frame(
    Dataset = "GSE69657",
    n = 30,
    Platform = "Array",
    Cohort = "XELOX neoadjuvant",
    Resistant_def = "noresponder",
    Sensitive_def = "responder",
    n_Resistant = n_r,
    n_Sensitive = n_s,
    n_Intermediate = NA,
    XELOX_specific = "PURE XELOX (n=30): highest relevance",
    Notes = "XELOX neoadjuvant. Ideal for pure XELOX resistance analysis."
  ),
  data.frame(
    Dataset = "GSE19860",
    n = 40,
    Platform = "Array (GPL570)",
    Cohort = "Capecitabine (FL) response",
    Resistant_def = "FL_Non_responder",
    Sensitive_def = "FL_Responder",
    n_Resistant = sum(fl_nonresponder, na.rm = TRUE),
    n_Sensitive = sum(fl_responder, na.rm = TRUE),
    n_Intermediate = NA,
    XELOX_specific = "Capecitabine (n=40)",
    Notes = "Some samples have combined BV (Bevacizumab) response."
  )
)

sum_tab <- do.call(rbind, summaries)
print(sum_tab, row.names = FALSE)

cat("\n\nTotal samples with resistance labels across all datasets:\n")
cat("  GSE39582:    45R + 119S = 164 labeled\n")
cat("  GSE72970:    variable (FOLFOX subset)\n")
cat("  GSE104645:   variable (FOLFOX/XELOX subset)\n")
cat("  GSE28702:    41R + 42S = 83 (all FOLFOX)  **clean binary**\n")
cat("  GSE69657:    17R + 13S = 30 (all XELOX)    **PURE XELOX**\n")
cat("  GSE19860:    variable (capecitabine)\n")

write.csv(sum_tab, file.path(RESULTS_TAB_DIR, "clinical_group_summary.csv"), row.names = FALSE)

cat("\n=== Clinical extraction v3 complete ===\n")
save.image(file = file.path(DATA_PROC_DIR, "02_clinical_groups.RData"))
