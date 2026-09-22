# ============================================================
# 28d_validate_gse83129.R
# Script 28d: GSE83129 External Validation (Ultimate Blind Test)
#
# Purpose:
#   Blind validation of the 13-pathway LASSO-PRS model on an
#   independent cohort: GSE83129 (GPL6244, HuGene-1.0-ST).
#   26 mCRC patients treated with 5FU+oxaliplatin (XELOX-like),
#   17 Responder / 9 Non-Responder.
#
# Workflow:
#   1. Parse series matrix → response labels (tumor-only, 26 samples)
#   2. Parse GPL6244 annot → probe-to-gene mapping
#   3. Aggregate to gene level
#   4. ssGSEA (same 44 pathways, same GMT)
#   5. Z-score pathway standardization
#   6. Apply 13-pathway PRS formula (Script 27 coefficients)
#   7. AUC calculation
#
# Input:
#   data/geo/GSE83129_series_matrix.txt.gz
#   data/geo/GPL6244.annot.gz
#   results/tables/prs_merged/merged_prs_coefficients.csv
#   C:/temp/h.all.v2024.1.symbols.gmt
#   C:/temp/c2.cp.kegg_legacy.v2024.1.symbols.gmt
#
# Output:
#   results/tables/gse83129_validation/*.csv, *.rds
#   results/figures/gse83129_validation/*.pdf
# ============================================================

Sys.setenv(TMPDIR = "C:/temp", TMP = "C:/temp", TEMP = "C:/temp")
.libPaths(c("C:/Rlibs", .libPaths()))

PROJECT_ROOT <- "C:/xelox_research"
DATA_GEO_DIR <- file.path(PROJECT_ROOT, "data", "geo")
RESULTS_DIR <- file.path(PROJECT_ROOT, "results", "tables", "gse83129_validation")
FIG_DIR <- file.path(PROJECT_ROOT, "results", "figures", "gse83129_validation")

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

cat("============================================================\n")
cat("Script 28d: GSE83129 External Validation (Blind Test)\n")
cat("============================================================\n\n")

# ============================================================
# 0. Load packages
# ============================================================
cat("=== [0] Loading packages ===\n\n")
suppressPackageStartupMessages({
  library(GSVA)
  library(GSEABase)
  library(pROC)
  library(ggplot2)
  library(data.table)
})
cat("  All packages loaded\n\n")

# ============================================================
# 1. Parse series matrix → response labels
# ============================================================
cat("=== [1] Parsing GSE83129 clinical data ===\n\n")

lines <- readLines(gzfile(file.path(DATA_GEO_DIR, "GSE83129_series_matrix.txt.gz")))

extract_field <- function(lines, prefix) {
  for (line in lines) {
    if (grepl(paste0("^", prefix), line)) {
      parts <- strsplit(line, "\t")[[1]]
      return(gsub('"', '', parts[-1]))
    }
  }
  return(NULL)
}

sample_ids <- extract_field(lines, "!Sample_geo_accession")
titles    <- extract_field(lines, "!Sample_title")
sources   <- extract_field(lines, "!Sample_source_name_ch1")

# characteristics_ch1 appears 4 times per sample (patient no, cancer type,
# treatment, response). The response line is the LAST match.
all_char_lines <- list()
for (line in lines) {
  if (grepl("^!Sample_characteristics_ch1", line)) {
    all_char_lines[[length(all_char_lines) + 1]] <- line
  }
}
# Use the 4th appearance (response)
resp_line <- all_char_lines[[4]]
parts <- strsplit(resp_line, "\t")[[1]]
chars_all <- gsub('"', '', parts[-1])

# Build sample-level data
sample_df <- data.frame(
  sample_id = sample_ids,
  title = titles,
  source = sources,
  response_raw = chars_all,
  stringsAsFactors = FALSE
)

# Parse response: "1st line best response: Responder" → 0 (sensitive)
#               "1st line best response: Non-Responder" → 1 (resistant)
#                               "na" → NA
# Order matters: check Non-Responder FIRST to avoid "Responder" matching "Non-Responder"
sample_df$response <- ifelse(grepl("Non-Responder", sample_df$response_raw), 1,
                      ifelse(grepl("1st line best response: Responder", sample_df$response_raw, fixed = TRUE), 0, NA))

# Keep only tumor samples (Adenocarcinoma) with valid response
tumor_idx <- sample_df$source == "Adenocarcinoma" & !is.na(sample_df$response)
val_df <- sample_df[tumor_idx, ]
cat(sprintf("  Tumor samples with valid response: %d\n", nrow(val_df)))
cat(sprintf("    Responder (sensitive): %d\n", sum(val_df$response == 0)))
cat(sprintf("    Non-Responder (resistant): %d\n", sum(val_df$response == 1)))

write.csv(val_df, file.path(RESULTS_DIR, "gse83129_clinical.csv"), row.names = FALSE)

# ============================================================
# 2. Load expression data
# ============================================================
cat("\n=== [2] Loading expression data ===\n\n")

# The series matrix expression data was parsed by 28c_parse_gse83129.py
# and saved as RDS via Python pickle. Let's read from the series matrix directly.
# Find where the data table starts
data_start <- 0
for (i in seq_along(lines)) {
  if (grepl("^!series_matrix_table_begin", lines[i])) {
    data_start <- i + 1
    break
  }
}

# Read data lines
data_lines <- list()
for (i in data_start:length(lines)) {
  s <- trimws(lines[i])
  if (s == "!series_matrix_table_end") break
  data_lines[[length(data_lines) + 1]] <- s
}

cat(sprintf("  Data rows: %d\n", length(data_lines)))

# Parse header
header_parts <- strsplit(data_lines[[1]], "\t")[[1]]
col_names <- gsub('"', '', header_parts)
cat(sprintf("  Columns: %d (ID_REF + %d samples)\n", length(col_names), length(col_names) - 1))

# Parse expression matrix
probe_ids <- character(length(data_lines) - 1)
expr_probe <- matrix(NA, nrow = length(data_lines) - 1, ncol = length(col_names) - 1)

for (i in 2:length(data_lines)) {
  parts <- strsplit(data_lines[[i]], "\t")[[1]]
  probe_ids[i - 1] <- gsub('"', '', parts[1])
  expr_probe[i - 1, ] <- as.numeric(parts[-1])
}
colnames(expr_probe) <- col_names[-1]  # sample IDs
rownames(expr_probe) <- probe_ids

cat(sprintf("  Expression matrix: %d probes x %d samples\n",
            nrow(expr_probe), ncol(expr_probe)))

# Save probe-level expression
saveRDS(expr_probe, file.path(RESULTS_DIR, "gse83129_probe_expression.rds"))

# ============================================================
# 3. Map probes to genes using GPL6244 annotation
# ============================================================
cat("\n=== [3] Probe-to-gene mapping (GPL6244) ===\n\n")

# Parse GPL6244 annotation
cat("  Parsing GPL6244 annotation...\n")

annot_file <- file.path(DATA_GEO_DIR, "GPL6244.annot.gz")
annot_conn <- gzfile(annot_file, "rt")

while (TRUE) {
  line <- readLines(annot_conn, n = 1)
  if (length(line) == 0) break
  if (grepl("!platform_table_begin", line)) break
}
# Read column headers
annot_cols <- strsplit(readLines(annot_conn, n = 1), "\t")[[1]]
cat(sprintf("  Annotation columns: %s\n", paste(annot_cols, collapse = ", ")))

# Read the table
probe_annot <- read.table(annot_conn, sep = "\t", header = FALSE,
                           col.names = annot_cols, fill = TRUE,
                           quote = "\"", stringsAsFactors = FALSE)
close(annot_conn)

cat(sprintf("  GPL6244 annotation: %d rows\n", nrow(probe_annot)))

# Filter to probes present in our expression data
idx <- match(probe_ids, probe_annot$ID)
matched <- probe_annot[idx, ]
valid <- !is.na(idx) & !is.na(matched$Gene.symbol) & matched$Gene.symbol != ""
cat(sprintf("  Probes with valid gene symbols: %d / %d\n", sum(valid), length(probe_ids)))

expr_matched <- expr_probe[valid, , drop = FALSE]
genes <- matched$Gene.symbol[valid]

# Aggregate to gene level (mean of multiple probes)
uniq_genes <- unique(genes)
expr_gene <- matrix(NA, nrow = length(uniq_genes), ncol = ncol(expr_matched))
rownames(expr_gene) <- uniq_genes
colnames(expr_gene) <- colnames(expr_matched)

for (i in seq_along(uniq_genes)) {
  g <- uniq_genes[i]
  probe_set <- which(genes == g)
  if (length(probe_set) == 1) {
    expr_gene[i, ] <- expr_matched[probe_set, ]
  } else {
    expr_gene[i, ] <- colMeans(expr_matched[probe_set, , drop = FALSE])
  }
}

cat(sprintf("  Gene-level expression: %d genes x %d samples\n",
            nrow(expr_gene), ncol(expr_gene)))

saveRDS(expr_gene, file.path(RESULTS_DIR, "gse83129_gene_expression.rds"))

# ============================================================
# 4. ssGSEA for the 44 pathways
# ============================================================
cat("\n=== [4] ssGSEA pathway scoring ===\n\n")

# 4a. Load GMT files
read_gmt <- function(gmt_file) {
  lines <- readLines(gmt_file, warn = FALSE)
  gs_list <- list()
  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    gs_name <- parts[1]
    genes <- parts[-(1:2)]
    genes <- genes[nchar(genes) > 0]
    gs_list[[gs_name]] <- genes
  }
  return(gs_list)
}

gmt_h <- "C:/temp/h.all.v2024.1.symbols.gmt"
gmt_kegg <- "C:/temp/c2.cp.kegg_legacy.v2024.1.symbols.gmt"

if (!file.exists(gmt_h) || !file.exists(gmt_kegg)) {
  stop("GMT files not found at C:/temp/")
}

hallmark_all <- read_gmt(gmt_h)
kegg_all <- read_gmt(gmt_kegg)

# 4b. Select the same 44 pathways as Script 14/16/27
target_pathways <- list()

# KEGG targets
kegg_targets <- c(
  "KEGG_NUCLEOTIDE_EXCISION_REPAIR", "KEGG_BASE_EXCISION_REPAIR",
  "KEGG_MISMATCH_REPAIR", "KEGG_HOMOLOGOUS_RECOMBINATION",
  "KEGG_P53_SIGNALING_PATHWAY", "KEGG_APOPTOSIS", "KEGG_CELL_CYCLE",
  "KEGG_ABC_TRANSPORTERS", "KEGG_GLUTATHIONE_METABOLISM",
  "KEGG_DRUG_METABOLISM_CYTOCHROME_P450", "KEGG_DRUG_METABOLISM_OTHER_ENZYMES",
  "KEGG_METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450",
  "KEGG_WNT_SIGNALING_PATHWAY", "KEGG_TGF_BETA_SIGNALING_PATHWAY",
  "KEGG_NOTCH_SIGNALING_PATHWAY", "KEGG_FOCAL_ADHESION",
  "KEGG_ECM_RECEPTOR_INTERACTION",
  "KEGG_MAPK_SIGNALING_PATHWAY", "KEGG_MTOR_SIGNALING_PATHWAY",
  "KEGG_JAK_STAT_SIGNALING_PATHWAY",
  "KEGG_CHEMOKINE_SIGNALING_PATHWAY", "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY", "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_COLORECTAL_CANCER", "KEGG_PATHWAYS_IN_CANCER"
)

# Hallmark targets
hallmark_targets <- c(
  "HALLMARK_DNA_REPAIR", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_APOPTOSIS", "HALLMARK_ANGIOGENESIS", "HALLMARK_HYPOXIA",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_GLYCOLYSIS",
  "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_P53_PATHWAY",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING", "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2", "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MTORC1_SIGNALING", "HALLMARK_TGF_BETA_SIGNALING"
)

for (id in kegg_targets) {
  if (id %in% names(kegg_all)) {
    target_pathways[[id]] <- kegg_all[[id]]
  }
}
for (id in hallmark_targets) {
  if (id %in% names(hallmark_all)) {
    target_pathways[[id]] <- hallmark_all[[id]]
  }
}

cat(sprintf("  Target pathways: %d\n", length(target_pathways)))

# Filter by genes present in GSE83129
common_genes <- rownames(expr_gene)
filtered_pathways <- list()
for (nm in names(target_pathways)) {
  gs <- target_pathways[[nm]]
  overlapping <- intersect(gs, common_genes)
  if (length(overlapping) >= 10) {
    filtered_pathways[[nm]] <- overlapping
  } else {
    cat(sprintf("  SKIP %s: only %d/%d genes present\n",
                nm, length(overlapping), length(gs)))
  }
}
cat(sprintf("  Usable gene sets: %d / %d\n", length(filtered_pathways), length(target_pathways)))

# 4c. Build GeneSetCollection
gs_list <- lapply(names(filtered_pathways), function(nm) {
  GeneSet(setName = nm, geneIds = filtered_pathways[[nm]])
})
gs_collection <- GeneSetCollection(gs_list)

# 4d. Run ssGSEA
cat("\n  Running ssGSEA...\n")
expr_mat <- as.matrix(expr_gene)
ps_gse83129 <- gsva(
  ssgseaParam(expr_mat, gs_collection,
    minSize = 10, maxSize = 500,
    verbose = TRUE)
)
cat(sprintf("  ssGSEA complete: %d pathways x %d samples\n",
            nrow(ps_gse83129), ncol(ps_gse83129)))

saveRDS(ps_gse83129, file.path(RESULTS_DIR, "gse83129_pathway_scores.rds"))
cat("  Saved pathway scores\n")

# ============================================================
# 5. Apply 13-pathway PRS model
# ============================================================
cat("\n=== [5] Computing PRS using 13-pathway model ===\n\n")

# 5a. Load model coefficients
coef_table <- read.csv(file.path(PROJECT_ROOT, "results", "tables",
                                   "prs_merged", "merged_prs_coefficients.csv"),
                        stringsAsFactors = FALSE)
cat(sprintf("  PRS model: %d pathways\n", nrow(coef_table)))
for (i in seq_len(nrow(coef_table))) {
  cat(sprintf("    %s: %.4f\n", coef_table$Pathway[i], coef_table$Coefficient[i]))
}

# 5b. Check which model pathways are present
model_pws <- coef_table$Pathway
present_pws <- intersect(model_pws, rownames(ps_gse83129))
missing_pws <- setdiff(model_pws, rownames(ps_gse83129))

cat(sprintf("\n  Pathways present: %d / %d\n", length(present_pws), nrow(coef_table)))
if (length(missing_pws) > 0) {
  cat("  MISSING pathways:\n")
  for (pw in missing_pws) cat(sprintf("    - %s\n", pw))
}

# 5c. Z-score standardize pathway scores
ps_z <- t(scale(t(ps_gse83129)))
cat(sprintf("  Z-scored pathway scores: %d x %d\n", nrow(ps_z), ncol(ps_z)))

# 5d. Compute PRS for ALL samples
prs_all <- rep(0, ncol(ps_z))
names(prs_all) <- colnames(ps_z)

for (i in seq_len(nrow(coef_table))) {
  pw <- coef_table$Pathway[i]
  coef_val <- coef_table$Coefficient[i]
  if (pw %in% rownames(ps_z)) {
    prs_all <- prs_all + coef_val * ps_z[pw, ]
    cat(sprintf("  + %s (%.4f): contribution range [%.3f, %.3f]\n",
                pw, coef_val, min(coef_val * ps_z[pw, ]), max(coef_val * ps_z[pw, ])))
  }
}

cat(sprintf("\n  PRS range: [%.3f, %.3f]\n", min(prs_all), max(prs_all)))

# Save all-sample PRS
prs_df <- data.frame(
  Sample = names(prs_all),
  PRS = prs_all,
  Title = val_df$title[match(names(prs_all), val_df$sample_id)],
  Source = sample_df$source[match(names(prs_all), sample_df$sample_id)],
  Response = sample_df$response[match(names(prs_all), sample_df$sample_id)],
  stringsAsFactors = FALSE
)
prs_df$ResponseLabel <- ifelse(is.na(prs_df$Response), "Unknown",
                        ifelse(prs_df$Response == 0, "Sensitive", "Resistant"))
write.csv(prs_df, file.path(RESULTS_DIR, "gse83129_prs_all_samples.csv"), row.names = FALSE)
cat("  Saved all-sample PRS\n")

# ============================================================
# 6. AUC on tumor-only samples
# ============================================================
cat("\n=== [6] AUC calculation ===\n\n")

# 6a. Tumor samples only
tumor_samples <- intersect(val_df$sample_id, names(prs_all))
cat(sprintf("  Tumor samples with PRS: %d\n", length(tumor_samples)))

prs_tumor <- prs_all[tumor_samples]
y_tumor <- val_df$response[match(tumor_samples, val_df$sample_id)]

cat(sprintf("  Outcome: Resistant=%d, Sensitive=%d\n",
            sum(y_tumor == 1), sum(y_tumor == 0)))

# 6b. ROC + AUC
if (length(unique(y_tumor)) == 2 && sd(prs_tumor) > 0) {
  roc_obj <- roc(y_tumor, prs_tumor, direction = "<", quiet = TRUE)
  val_auc <- as.numeric(auc(roc_obj))
  ci_auc <- ci.auc(roc_obj)
  
  cat(sprintf("\n  *** GSE83129 VALIDATION AUC: %.3f (95%% CI: %.3f - %.3f) ***\n",
              val_auc, ci_auc[1], ci_auc[3]))
  
  # 6c. ROC plot
  roc_df <- data.frame(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities
  )
  
  p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
    geom_line(color = "#2166AC", linewidth = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    annotate("text", x = 0.75, y = 0.25,
             label = paste0("AUC = ", sprintf("%.3f", val_auc),
                        "\n95% CI: ", sprintf("%.3f", ci_auc[1]),
                        " - ", sprintf("%.3f", ci_auc[3])),
             hjust = 0, size = 4.5) +
    labs(title = "GSE83129: 13-Pathway PRS Validation",
         subtitle = paste0("Blind test: ", length(tumor_samples),
                          " mCRC (5FU+Oxaliplatin) | ",
                          sum(y_tumor == 0), "R / ", sum(y_tumor == 1), "NR"),
         x = "1 - Specificity (FPR)",
         y = "Sensitivity (TPR)") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 10, color = "gray40"))
  
  ggsave(file.path(FIG_DIR, "gse83129_roc_curve.pdf"), p_roc, width = 7, height = 6)
  ggsave(file.path(FIG_DIR, "gse83129_roc_curve.png"), p_roc, width = 7, height = 6, dpi = 150)
  cat("  ROC curve saved\n")
  
  # 6d. PRS boxplot by response
  box_df <- data.frame(
    PRS = prs_tumor,
    Response = ifelse(y_tumor == 0, "Sensitive", "Resistant")
  )
  
  p_box <- ggplot(box_df, aes(x = Response, y = PRS, fill = Response)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +
    scale_fill_manual(values = c("Sensitive" = "#4393C3", "Resistant" = "#D6604D")) +
    labs(title = "PRS by Response Status (GSE83129)",
         x = "", y = "Pathway-based Risk Score (PRS)") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))
  
  ggsave(file.path(FIG_DIR, "gse83129_prs_boxplot.pdf"), p_box, width = 5.5, height = 5)
  ggsave(file.path(FIG_DIR, "gse83129_prs_boxplot.png"), p_box, width = 5.5, height = 5, dpi = 150)
  cat("  Boxplot saved\n")
  
  # 6e. Wilcoxon test
  wt <- wilcox.test(prs_tumor ~ y_tumor)
  cat(sprintf("  Wilcoxon p-value: %.4f\n", wt$p.value))
  
} else {
  val_auc <- NA
  ci_auc <- c(NA, NA, NA)
  cat("  WARNING: Cannot compute AUC (insufficient data or variation)\n")
}

# ============================================================
# 7. Compare with other cohorts
# ============================================================
cat("\n=== [7] Validation summary comparison ===\n\n")

# Load previous validation results
prev_val <- read.csv(file.path(PROJECT_ROOT, "results", "tables",
                                 "prs_merged", "merged_validation_summary.csv"),
                      stringsAsFactors = FALSE)

combined <- rbind(
  prev_val,
  data.frame(Cohort = "GSE83129 (blind)",
              N = length(tumor_samples),
              N_resistant = sum(y_tumor == 1),
              AUC = round(val_auc, 4),
              stringsAsFactors = FALSE)
)

cat("\n  Combined validation summary:\n")
print(combined)
write.csv(combined, file.path(RESULTS_DIR, "validation_summary_with_gse83129.csv"), row.names = FALSE)

# ============================================================
# 8. Pathway activity profile comparison
# ============================================================
cat("\n=== [8] Pathway activity profile ===\n\n")

# Compare the 13 model pathways between Responders and Non-Responders
model_ps <- ps_z[intersect(model_pws, rownames(ps_z)), tumor_samples, drop = FALSE]
y_order <- order(y_tumor)

heatmap_data <- model_ps[, y_order]
colnames(heatmap_data) <- paste0(val_df$title[match(colnames(heatmap_data), val_df$sample_id)],
                                  "_", ifelse(y_tumor[y_order] == 0, "S", "R"))

# Save heatmap data
write.csv(cbind(Pathway = rownames(heatmap_data), as.data.frame(heatmap_data)),
          file.path(RESULTS_DIR, "gse83129_13pathway_heatmap.csv"), row.names = FALSE)

cat("  GSE83129 validation complete.\n\n")

# Final report
sink(file.path(RESULTS_DIR, "gse83129_validation_report.txt"))
cat("GSE83129 Blind Validation Report\n")
cat("================================\n\n")
cat(sprintf("Dataset: GSE83129 (GPL6244, HuGene-1.0-ST)\n"))
cat(sprintf("Samples: %d total (26 tumor + 10 normal laser-microdissected)\n", ncol(expr_probe)))
cat(sprintf("Tumor samples with oxaliplatin response: %d\n", nrow(val_df)))
cat(sprintf("  Responder (sensitive): %d\n", sum(val_df$response == 0)))
cat(sprintf("  Non-Responder (resistant): %d\n", sum(val_df$response == 1)))
cat(sprintf("Treatment: %s\n\n", "First-line 5FU + oxaliplatin (XELOX-like)"))
cat(sprintf("Genes mapped: %d\n", nrow(expr_gene)))
cat(sprintf("Pathways scored: %d\n", nrow(ps_gse83129)))
cat(sprintf("PRS model pathways: %d\n", nrow(coef_table)))
cat(sprintf("Model pathways present: %d / %d\n\n", length(present_pws), nrow(coef_table)))
cat(sprintf("Validation AUC: %.3f\n", val_auc))
if (!is.na(val_auc)) {
  cat(sprintf("95%% CI: %.3f - %.3f\n", ci_auc[1], ci_auc[3]))
}
cat(sprintf("Wilcoxon p-value: %.4f\n\n", ifelse(exists("wt"), wt$p.value, NA)))
cat("Comparison with training:\n")
cat(sprintf("  Training AUC (GSE39582+GSE19860 merged): 0.764\n"))
cat(sprintf("  GSE83129 (blind): %.3f\n", val_auc))
if (!is.na(val_auc)) {
  delta <- val_auc - 0.764
  cat(sprintf("  Delta: %+.3f\n", delta))
}
sink()

cat("  Report saved to:", file.path(RESULTS_DIR, "gse83129_validation_report.txt"), "\n")
cat("\n============================================================\n")
cat("Script 28d Complete\n")
cat("============================================================\n")
