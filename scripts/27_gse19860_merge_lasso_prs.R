# ============================================================
# 27_GSE19860_MERGE_LASSO_PRS.R
# Script 27: GSE19860 Gene-Level Merge + ComBat + ssGSEA + LASSO + Bootstrap
#
# Purpose:
#   Augment GSE39582 XELOX training with GSE19860 (40 mFOLFOX6, GPL570)
#   to increase training sample size and stabilize LASSO-PRS selection.
#   Key insight: merge at GENE level (not pathway level), then ComBat,
#   then ssGSEA, to preserve subtle co-expression signals.
#
# Step 1: Gene-level merge of GSE39582 (n=164 XELOX) + GSE19860 (n=40)
# Step 2: ComBat batch correction (2 batches)
# Step 3: ssGSEA on merged ComBat-corrected matrix (44 pathways)
# Step 4: LASSO-Logistic PRS (AUC pre-filter + cross-validated LASSO)
# Step 5: Bootstrap stability analysis (B=1000)
# Step 6: Validation in existing cohorts (cross-platform)
#
# Input:
#   data/geo/GSE19860_expression.rds            — probe-level (GPL570)
#   data/geo/GSE19860_series_matrix.txt.gz      — clinical labels
#   data/geo/GPL570_probe_gene_map.csv          — probe-to-gene map
#   results/tables/pathway_activity/GSE39582_gene_expression.rds
#   results/tables/GSE39582_xelox_groups.csv    — XELOX labels
#   C:/temp/c2.cp.kegg_legacy.v2024.1.symbols.gmt
#   C:/temp/h.all.v2024.1.symbols.gmt
#
# Output:
#   results/tables/prs_merged/*.csv             — model coefficients, PRS, AUC
#   results/figures/prs_merged/*.pdf            — ROC, stability, PCA
# ============================================================

Sys.setenv(TMPDIR = "C:/temp", TMP = "C:/temp", TEMP = "C:/temp")
.libPaths(c("C:/Rlibs", .libPaths()))

PROJECT_ROOT <- "C:/xelox_research"
DATA_GEO_DIR <- file.path(PROJECT_ROOT, "data", "geo")
RESULTS_TAB_DIR <- file.path(PROJECT_ROOT, "results", "tables")
RESULTS_FIG_DIR <- file.path(PROJECT_ROOT, "results", "figures")
PA_DIR <- file.path(RESULTS_TAB_DIR, "pathway_activity")

# Output directories
OUT_TAB_DIR <- file.path(RESULTS_TAB_DIR, "prs_merged")
OUT_FIG_DIR <- file.path(RESULTS_FIG_DIR, "prs_merged")
dir.create(OUT_TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_FIG_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

cat("============================================================\n")
cat("Script 27: GSE19860 Gene-Level Merge + ComBat + LASSO-PRS\n")
cat("============================================================\n\n")

# ============================================================
# 0. Load packages
# ============================================================
cat("=== [0] Loading packages ===\n\n")
suppressPackageStartupMessages({
  library(sva)         # ComBat
  library(GSVA)        # ssGSEA
  library(GSEABase)    # GeneSetCollection
  library(glmnet)      # LASSO
  library(pROC)        # AUC
  library(ggplot2)     # plotting
  library(data.table)  # fast I/O
})
cat("  All packages loaded\n\n")

# ============================================================
# 1. Map GSE19860 probes to genes
# ============================================================
cat("=== [1] Processing GSE19860 (GPL570) ===\n\n")

# 1a. Load clinical labels
lines <- readLines(gzfile(file.path(DATA_GEO_DIR, "GSE19860_series_matrix.txt.gz")))
sample_ids <- c()
response <- c()

for (i in seq_along(lines)) {
  if (grepl('^!Sample_geo_accession', lines[i])) {
    parts <- strsplit(lines[i], '\t')[[1]]
    sample_ids <- gsub('"', '', parts[-1])
  }
  if (grepl('^!Sample_characteristics_ch1', lines[i])) {
    parts <- strsplit(lines[i], '\t')[[1]][-1]
    if (any(grepl('treatment response', parts, fixed = TRUE))) {
      response <- parts
    }
  }
}

# Parse: "treatment response: FL_Responder" or "FL_Non_responder"
resp_clean <- sub('treatment response: ', '', response)
resp_clean <- sub(',.*', '', resp_clean)   # remove BV status
resp_clean <- gsub('"', '', resp_clean)
resp_clean <- trimws(resp_clean)

g19860_label <- ifelse(resp_clean == "FL_Responder", 0,
                       ifelse(resp_clean == "FL_Non_responder", 1, NA))
names(g19860_label) <- sample_ids
cat(sprintf("  GSE19860 labels: Responder=%d, Non_responder=%d\n",
            sum(g19860_label == 0, na.rm = TRUE),
            sum(g19860_label == 1, na.rm = TRUE)))
cat(sprintf("  Samples: %d\n\n", length(g19860_label)))

# 1b. Load probe-level expression
expr_19860_probe <- readRDS(file.path(DATA_GEO_DIR, "GSE19860_expression.rds"))
cat(sprintf("  GSE19860 probe-level dim: %d probes x %d samples\n",
            nrow(expr_19860_probe), ncol(expr_19860_probe)))

# 1c. Probe to gene mapping
probe_map <- read.csv(file.path(DATA_GEO_DIR, "GPL570_probe_gene_map.csv"),
                       stringsAsFactors = FALSE)
cat(sprintf("  GPL570 probe map: %d probes, %d unique genes\n",
            nrow(probe_map), length(unique(probe_map$gene_symbol))))

probe_ids <- rownames(expr_19860_probe)
idx <- match(probe_ids, probe_map$probe_id)
matched <- probe_map[idx, ]
valid <- !is.na(idx)
expr_matched <- expr_19860_probe[valid, , drop = FALSE]
genes <- matched$gene_symbol[valid]

has_gene <- !is.na(genes) & genes != ""
expr_matched <- expr_matched[has_gene, , drop = FALSE]
genes <- genes[has_gene]

# Aggregate to gene level (mean of multiple probes)
uniq_genes <- unique(genes)
expr_19860 <- matrix(NA, nrow = length(uniq_genes), ncol = ncol(expr_matched))
rownames(expr_19860) <- uniq_genes
colnames(expr_19860) <- colnames(expr_matched)

for (i in seq_along(uniq_genes)) {
  g <- uniq_genes[i]
  probe_set <- which(genes == g)
  if (length(probe_set) == 1) {
    expr_19860[i, ] <- expr_matched[probe_set, ]
  } else {
    expr_19860[i, ] <- colMeans(expr_matched[probe_set, , drop = FALSE])
  }
}

cat(sprintf("  GSE19860 gene-level dim: %d genes x %d samples\n\n",
            nrow(expr_19860), ncol(expr_19860)))

# Reorder to match label order
common_samples_19860 <- intersect(colnames(expr_19860), names(g19860_label))
expr_19860 <- expr_19860[, common_samples_19860, drop = FALSE]
g19860_label <- g19860_label[common_samples_19860]

# ============================================================
# 2. Load GSE39582 XELOX gene expression + labels
# ============================================================
cat("=== [2] Loading GSE39582 XELOX data ===\n\n")

expr_39582 <- readRDS(file.path(PA_DIR, "GSE39582_gene_expression.rds"))
cat(sprintf("  GSE39582 gene-level dim: %d genes x %d samples\n",
            nrow(expr_39582), ncol(expr_39582)))

groups_39582 <- read.csv(file.path(RESULTS_TAB_DIR, "GSE39582_xelox_groups.csv"),
                          stringsAsFactors = FALSE)

# Extract XELOX-like samples (binary labels only)
is_xelox <- groups_39582[["group"]] %in% c("sensitive", "resistant")
xelox_samples <- groups_39582$sample_id[is_xelox]
g39582_label <- ifelse(groups_39582[["group"]][is_xelox] == "resistant", 1, 0)
names(g39582_label) <- xelox_samples

cat(sprintf("  GSE39582 XELOX: Sensitive=%d, Resistant=%d\n",
            sum(g39582_label == 0), sum(g39582_label == 1)))

# Match to expression matrix
common_39582 <- intersect(xelox_samples, colnames(expr_39582))
cat(sprintf("  Samples matching expression: %d / %d\n",
            length(common_39582), length(xelox_samples)))

expr_39582_xelox <- expr_39582[, common_39582, drop = FALSE]
g39582_label <- g39582_label[common_39582]

# ============================================================
# 3. Merge at gene level + ComBat batch correction
# ============================================================
cat("\n=== [3] Gene-level merge + ComBat ===\n\n")

# 3a. Find common genes
common_genes <- intersect(rownames(expr_39582_xelox), rownames(expr_19860))
cat(sprintf("  Common genes: %d\n", length(common_genes)))

# 3b. Build merged expression matrix
merged_expr <- cbind(
  expr_39582_xelox[common_genes, , drop = FALSE],
  expr_19860[common_genes, , drop = FALSE]
)
cat(sprintf("  Merged expression: %d genes x %d samples\n",
            nrow(merged_expr), ncol(merged_expr)))

# 3c. Batch vector
batch <- c(rep("GSE39582", ncol(expr_39582_xelox)),
           rep("GSE19860", ncol(expr_19860)))
cat(sprintf("  Batches: GSE39582=%d, GSE19860=%d\n",
            sum(batch == "GSE39582"), sum(batch == "GSE19860")))

# 3d. Combined label vector
merged_y <- c(g39582_label, g19860_label)
cat(sprintf("  Combined: Resistant=%d, Sensitive=%d\n",
            sum(merged_y == 1), sum(merged_y == 0)))

# 3e. Apply ComBat
cat("  Running ComBat...\n")
merged_combat <- ComBat(dat = merged_expr, batch = batch, mod = NULL,
                         par.prior = TRUE, prior.plots = FALSE)
cat("  ComBat complete\n\n")

# 3f. PCA diagnostic (pre vs post)
cat("  Generating PCA diagnostics...\n")
pca_pre <- prcomp(t(merged_expr), center = TRUE, scale. = FALSE)
pca_post <- prcomp(t(merged_combat), center = TRUE, scale. = FALSE)
pca_df <- data.frame(
  PC1_pre = pca_pre$x[, 1], PC2_pre = pca_pre$x[, 2],
  PC1_post = pca_post$x[, 1], PC2_post = pca_post$x[, 2],
  Batch = batch,
  Response = ifelse(merged_y == 1, "Resistant", "Sensitive"),
  Dataset = ifelse(batch == "GSE39582", "GSE39582", "GSE19860")
)

p_pre <- ggplot(pca_df, aes(x = PC1_pre, y = PC2_pre, color = Dataset, shape = Response)) +
  geom_point(size = 2.5, alpha = 0.8) +
  labs(title = "PCA Before ComBat", x = "PC1", y = "PC2") +
  theme_minimal() +
  scale_color_manual(values = c("GSE39582" = "#2166AC", "GSE19860" = "#B2182B"))
ggsave(file.path(OUT_FIG_DIR, "pca_pre_combat.pdf"), p_pre, width = 7, height = 5.5)

p_post <- ggplot(pca_df, aes(x = PC1_post, y = PC2_post, color = Dataset, shape = Response)) +
  geom_point(size = 2.5, alpha = 0.8) +
  labs(title = "PCA After ComBat", x = "PC1", y = "PC2") +
  theme_minimal() +
  scale_color_manual(values = c("GSE39582" = "#2166AC", "GSE19860" = "#B2182B"))
ggsave(file.path(OUT_FIG_DIR, "pca_post_combat.pdf"), p_post, width = 7, height = 5.5)

# Save merged data
saveRDS(merged_combat, file.path(OUT_TAB_DIR, "merged_expression_combat.rds"))
write.csv(data.frame(Sample = colnames(merged_combat), Batch = batch,
                      Response = merged_y, stringsAsFactors = FALSE),
          file.path(OUT_TAB_DIR, "merged_sample_info.csv"), row.names = FALSE)
cat("  Saved merged ComBat expression + sample info\n\n")

# ============================================================
# 4. ssGSEA on merged ComBat-corrected matrix
# ============================================================
cat("=== [4] ssGSEA on merged ComBat-corrected matrix ===\n\n")

# 4a. Load GMT
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

# 4b. Select the same 44 pathways as Script 14/16
target_pathways <- list()

# KEGG: DNA repair
kegg_targets <- c(
  "KEGG_NUCLEOTIDE_EXCISION_REPAIR", "KEGG_BASE_EXCISION_REPAIR",
  "KEGG_MISMATCH_REPAIR", "KEGG_HOMOLOGOUS_RECOMBINATION",
  "KEGG_P53_SIGNALING_PATHWAY", "KEGG_APOPTOSIS", "KEGG_CELL_CYCLE"
)
# KEGG: Drug metabolism
kegg_targets <- c(kegg_targets,
  "KEGG_ABC_TRANSPORTERS", "KEGG_GLUTATHIONE_METABOLISM",
  "KEGG_DRUG_METABOLISM_CYTOCHROME_P450", "KEGG_DRUG_METABOLISM_OTHER_ENZYMES",
  "KEGG_METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450"
)
# KEGG: EMT/Metastasis
kegg_targets <- c(kegg_targets,
  "KEGG_WNT_SIGNALING_PATHWAY", "KEGG_TGF_BETA_SIGNALING_PATHWAY",
  "KEGG_NOTCH_SIGNALING_PATHWAY", "KEGG_FOCAL_ADHESION",
  "KEGG_ECM_RECEPTOR_INTERACTION"
)
# KEGG: Signaling
kegg_targets <- c(kegg_targets,
  "KEGG_MAPK_SIGNALING_PATHWAY", "KEGG_MTOR_SIGNALING_PATHWAY",
  "KEGG_JAK_STAT_SIGNALING_PATHWAY"
)
# KEGG: Immune
kegg_targets <- c(kegg_targets,
  "KEGG_CHEMOKINE_SIGNALING_PATHWAY", "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY", "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY"
)
# KEGG: Cancer
kegg_targets <- c(kegg_targets, "KEGG_COLORECTAL_CANCER", "KEGG_PATHWAYS_IN_CANCER")

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

# Collect matched gene sets
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

cat(sprintf("  Selected %d target pathways\n", length(target_pathways)))

# 4c. Prepare GeneSetCollection
gs_list <- lapply(names(target_pathways), function(nm) {
  GeneSet(setName = nm, geneIds = target_pathways[[nm]])
})
gs_collection <- GeneSetCollection(gs_list)

# 4d. Run ssGSEA (GSVA v2.x API: ssgseaParam + gsva)
cat("  Running ssGSEA...\n")
merged_es <- gsva(
  ssgseaParam(merged_combat, gs_collection,
    minSize = 10, maxSize = 500,
    verbose = TRUE)
)
cat(sprintf("  ssGSEA complete: %d pathways x %d samples\n",
            nrow(merged_es), ncol(merged_es)))

# 4e. Save merged pathway scores
saveRDS(merged_es, file.path(OUT_TAB_DIR, "merged_pathway_scores.rds"))
cat("  Saved merged pathway scores\n\n")

# ============================================================
# 5. LASSO-Logistic PRS modeling
# ============================================================
cat("=== [5] LASSO-Logistic PRS ===\n\n")

train_ps <- merged_es  # 44 pathways x 204 samples
train_y <- merged_y

# Ensure sample order matches
common_samples <- intersect(colnames(train_ps), names(train_y))
train_ps <- train_ps[, common_samples, drop = FALSE]
train_y <- train_y[common_samples]

cat(sprintf("  Training: %d pathways x %d samples\n",
            nrow(train_ps), ncol(train_ps)))
cat(sprintf("  Outcome: Resistant=%d, Sensitive=%d\n",
            sum(train_y == 1), sum(train_y == 0)))

# Transpose: samples x pathways
train_x <- t(train_ps)
pathway_names <- colnames(train_x)

# ============================================================
# 5a. Step 1: Univariate AUC pre-filter
# ============================================================
cat("\n  Step 1: Univariate AUC pre-filter\n")

compute_auc <- function(x, y) {
  tryCatch({
    roc_obj <- roc(y, x, direction = "<", quiet = TRUE)
    as.numeric(auc(roc_obj))
  }, error = function(e) NA)
}

auc_vec <- apply(train_ps, 1, function(pw) compute_auc(pw, train_y))
names(auc_vec) <- pathway_names

auc_filter <- abs(auc_vec - 0.5) > 0.05
auc_filter[is.na(auc_filter)] <- FALSE

cat(sprintf("  Passed univariate filter (|AUC-0.5| > 0.05): %d / %d pathways\n",
            sum(auc_filter), length(auc_filter)))

auc_sorted <- sort(auc_vec[auc_filter], decreasing = TRUE)
cat("\n  Top pro-resistance (AUC > 0.55):\n")
top_res <- head(auc_sorted[auc_sorted > 0.5], 10)
for (i in seq_along(top_res)) {
  cat(sprintf("    %2d. %s (AUC = %.3f)\n", i, names(top_res)[i], top_res[i]))
}
cat("\n  Top pro-sensitive (AUC < 0.45):\n")
top_sens <- head(auc_sorted[auc_sorted < 0.5], 10)
for (i in seq_along(top_sens)) {
  cat(sprintf("    %2d. %s (AUC = %.3f)\n", i, names(top_sens)[i], top_sens[i]))
}

train_x_filtered <- train_x[, auc_filter, drop = FALSE]
cat(sprintf("\n  Filtered X: %d samples x %d features\n",
            nrow(train_x_filtered), ncol(train_x_filtered)))

# Save AUC table
auc_df <- data.frame(Pathway = names(auc_vec), AUC = auc_vec,
                      PassFilter = auc_filter, stringsAsFactors = FALSE)
write.csv(auc_df, file.path(OUT_TAB_DIR, "merged_auc_filter.csv"), row.names = FALSE)

# ============================================================
# 5b. Step 2: LASSO-Logistic
# ============================================================
cat("\n  Step 2: LASSO-Logistic feature selection\n")

train_x_scaled <- scale(train_x_filtered)

set.seed(42)
cv_lasso <- cv.glmnet(
  x = train_x_scaled,
  y = train_y,
  family = "binomial",
  alpha = 1,
  nfolds = 10,
  type.measure = "deviance"
)

# Try lambda.1se first, fall back to lambda.min
coef_lasso <- as.matrix(coef(cv_lasso, s = "lambda.1se"))
selected_idx <- which(coef_lasso[-1, 1] != 0)
selected_pathways <- pathway_names[auc_filter][selected_idx]

cat(sprintf("  Selected at lambda.1se: %d pathways\n", length(selected_pathways)))

if (length(selected_pathways) == 0) {
  cat("  WARNING: No pathways at lambda.1se, falling back to lambda.min\n")
  coef_lasso <- as.matrix(coef(cv_lasso, s = "lambda.min"))
  selected_idx <- which(coef_lasso[-1, 1] != 0)
  selected_pathways <- pathway_names[auc_filter][selected_idx]
  cat(sprintf("  Selected at lambda.min: %d pathways\n", length(selected_pathways)))
}

# Coefficient table
if (length(selected_pathways) > 0) {
  coef_table <- data.frame(
    Pathway = selected_pathways,
    Coefficient = coef_lasso[selected_idx + 1, 1],
    stringsAsFactors = FALSE
  )
  cat("\n  PRS coefficients:\n")
  for (i in seq_len(nrow(coef_table))) {
    cat(sprintf("    %s: %.4f\n", coef_table$Pathway[i], coef_table$Coefficient[i]))
  }
} else {
  cat("\n  WARNING: No pathways selected! Creating empty model.\n")
  coef_table <- data.frame(Pathway = character(0), Coefficient = numeric(0),
                            stringsAsFactors = FALSE)
}

write.csv(coef_table, file.path(OUT_TAB_DIR, "merged_prs_coefficients.csv"), row.names = FALSE)

# Save CV curve
pdf(file.path(OUT_FIG_DIR, "lasso_cv_merged.pdf"), width = 7, height = 5.5)
plot(cv_lasso)
dev.off()
cat("  Saved LASSO CV curve\n")

# ============================================================
# 5c. Compute PRS for training set
# ============================================================
cat("\n  Computing PRS for training set\n")

if (nrow(coef_table) > 0) {
  train_ps_z <- t(scale(t(train_ps)))
  prs_train <- rep(0, ncol(train_ps_z))
  names(prs_train) <- colnames(train_ps_z)
  for (i in seq_len(nrow(coef_table))) {
    pw <- coef_table$Pathway[i]
    coef_val <- coef_table$Coefficient[i]
    if (pw %in% rownames(train_ps_z)) {
      prs_train <- prs_train + coef_val * train_ps_z[pw, ]
    }
  }
} else {
  prs_train <- rep(0, ncol(train_ps))
  names(prs_train) <- colnames(train_ps)
}

# Training AUC
if (length(unique(train_y)) == 2 && sd(prs_train) > 0) {
  roc_train <- roc(train_y, prs_train, direction = "<", quiet = TRUE)
  train_auc <- as.numeric(auc(roc_train))
  cat(sprintf("  Training AUC: %.3f\n", train_auc))
} else {
  train_auc <- NA
  cat("  Training AUC: NA (insufficient variation)\n")
}

prs_train_df <- data.frame(
  Sample = names(prs_train), PRS = prs_train,
  Response = ifelse(train_y == 1, "Resistant", "Sensitive"),
  Dataset = batch,
  stringsAsFactors = FALSE
)
write.csv(prs_train_df, file.path(OUT_TAB_DIR, "merged_prs_train.csv"), row.names = FALSE)

# ============================================================
# 6. Validate in existing cohorts
# ============================================================
cat("\n=== [6] Validation in existing cohorts ===\n\n")

# Load existing pathway scores
load_ps <- function(prefix) {
  f <- file.path(PA_DIR, paste0(prefix, "_pathway_scores.rds"))
  if (file.exists(f)) {
    ps <- readRDS(f)
    cat(sprintf("  %s: %d x %d\n", prefix, nrow(ps), ncol(ps)))
    return(ps)
  }
  cat(sprintf("  %s: NOT FOUND\n", prefix))
  return(NULL)
}

validate_cohorts <- list(
  GSE28702 = load_ps("GSE28702"),
  GSE69657 = load_ps("GSE69657"),
  GSE72970 = load_ps("GSE72970"),
  GSE104645 = load_ps("GSE104645")
)

# Load clinical labels for validation cohorts
get_val_labels <- function(prefix) {
  f <- file.path(RESULTS_TAB_DIR, paste0(prefix, "_clinical_data.csv"))
  if (!file.exists(f)) return(NULL)
  clin <- read.csv(f, row.names = 1, stringsAsFactors = FALSE)
  return(clin)
}

# Known label formats from Script 16
val_results <- data.frame(Cohort = character(), N = integer(),
                           N_resistant = integer(), AUC = numeric(),
                           stringsAsFactors = FALSE)

for (nm in names(validate_cohorts)) {
  ps_val <- validate_cohorts[[nm]]
  if (is.null(ps_val)) next
  
  clin <- get_val_labels(nm)
  if (is.null(clin)) {
    cat(sprintf("  %s: no clinical data, skipping\n", nm))
    next
  }
  
  # Extract labels (uses same logic as Script 16)
  if (nm == "GSE28702") {
    resp <- clin[["mfolfox6.ch1"]]
    names(resp) <- rownames(clin)
    val_y <- ifelse(resp == "responder", 0, ifelse(resp == "non-responder", 1, NA))
  } else if (nm == "GSE69657") {
    resp <- clin[["chemoresponse.ch1"]]
    names(resp) <- rownames(clin)
    val_y <- ifelse(resp == "responder", 0, ifelse(resp == "noresponder", 1, NA))
  } else if (nm == "GSE72970") {
    resp <- clin[["response.status.ch1"]]
    names(resp) <- rownames(clin)
    val_y <- ifelse(resp == "R", 0, ifelse(resp == "NR", 1, NA))
  } else if (nm == "GSE104645") {
    regimen_col <- "X1st.line.chemotherapy.regimens.ch1"
    resp_col <- "best.response.of.1st.line.chemotherapy.ch1"
    is_oxali <- grepl("XELOX|FOLFOX|SOX|CAPOX", clin[[regimen_col]], ignore.case = TRUE)
    oxali_samples <- rownames(clin)[is_oxali]
    resp <- clin[[resp_col]]
    names(resp) <- rownames(clin)
    val_y <- rep(NA, length(oxali_samples))
    names(val_y) <- oxali_samples
    val_y[resp[oxali_samples] %in% c("Complete response", "Partial response")] <- 0
    val_y[resp[oxali_samples] %in% c("Progressive disease", "Stable disease")] <- 1
  } else {
    next
  }
  
  if (nrow(coef_table) > 0) {
    ps_z <- t(scale(t(ps_val)))
    common_pw <- intersect(coef_table$Pathway, rownames(ps_z))
    
    if (length(common_pw) > 0) {
      prs_val <- rep(0, ncol(ps_z))
      names(prs_val) <- colnames(ps_z)
      for (i in seq_len(nrow(coef_table))) {
        pw <- coef_table$Pathway[i]
        coef_val <- coef_table$Coefficient[i]
        if (pw %in% rownames(ps_z)) {
          prs_val <- prs_val + coef_val * ps_z[pw, ]
        }
      }
      
      common_val <- intersect(names(prs_val), names(val_y[!is.na(val_y)]))
      if (length(common_val) >= 10 && length(unique(val_y[common_val])) == 2) {
        roc_val <- roc(val_y[common_val], prs_val[common_val], direction = "<", quiet = TRUE)
        val_auc <- as.numeric(auc(roc_val))
        cat(sprintf("  %s AUC: %.3f (n=%d, R=%d, S=%d)\n",
                    nm, val_auc, length(common_val),
                    sum(val_y[common_val] == 1), sum(val_y[common_val] == 0)))
      } else {
        val_auc <- NA
        cat(sprintf("  %s: insufficient data (n=%d)\n", nm, length(common_val)))
      }
    } else {
      val_auc <- NA
      cat(sprintf("  %s: no common pathways\n", nm))
    }
  } else {
    val_auc <- NA
    cat("  %s: no model (empty coefficients)\n", nm)
  }
  
  val_results <- rbind(val_results, data.frame(
    Cohort = nm,
    N = length(common_val),
    N_resistant = sum(val_y[common_val] == 1, na.rm = TRUE),
    AUC = round(val_auc, 4),
    stringsAsFactors = FALSE
  ))
}

# Also include training
val_results <- rbind(
  data.frame(Cohort = "Training", N = length(train_y),
              N_resistant = sum(train_y == 1),
              AUC = round(train_auc, 4), stringsAsFactors = FALSE),
  val_results
)

cat("\n  Validation summary:\n")
print(val_results)
write.csv(val_results, file.path(OUT_TAB_DIR, "merged_validation_summary.csv"), row.names = FALSE)

# ============================================================
# 7. Bootstrap stability analysis
# ============================================================
cat("\n=== [7] Bootstrap stability analysis ===\n\n")

B <- 1000
n_train <- length(train_y)
n_features <- ncol(train_x_filtered)

# Storage
boot_coef_matrix <- matrix(0, nrow = B, ncol = n_features)
colnames(boot_coef_matrix) <- colnames(train_x_filtered)
boot_selected_count <- 0
boot_lambda_min_count <- 0

set.seed(42)
for (b in 1:B) {
  # Stratified bootstrap (preserve class proportions)
  idx_res <- which(train_y == 1)
  idx_sens <- which(train_y == 0)
  boot_idx <- c(
    sample(idx_res, size = length(idx_res), replace = TRUE),
    sample(idx_sens, size = length(idx_sens), replace = TRUE)
  )
  
  x_boot <- train_x_scaled[boot_idx, , drop = FALSE]
  y_boot <- train_y[boot_idx]
  
  cv_boot <- tryCatch(
    cv.glmnet(x = x_boot, y = y_boot, family = "binomial",
               alpha = 1, nfolds = 5, type.measure = "deviance"),
    error = function(e) NULL
  )
  
  if (is.null(cv_boot)) next
  
  # Try lambda.1se
  coef_boot <- as.matrix(coef(cv_boot, s = "lambda.1se"))
  nz_1se <- sum(coef_boot[-1, 1] != 0)
  
  if (nz_1se == 0) {
    # Fall back to lambda.min
    coef_boot <- as.matrix(coef(cv_boot, s = "lambda.min"))
    nz_min <- sum(coef_boot[-1, 1] != 0)
    if (nz_min > 0) boot_lambda_min_count <- boot_lambda_min_count + 1
  }
  
  boot_coef_matrix[b, ] <- coef_boot[-1, 1]
  if (sum(coef_boot[-1, 1] != 0) > 0) boot_selected_count <- boot_selected_count + 1
  
  if (b %% 100 == 0) cat(sprintf("  Bootstrap %d / %d\n", b, B))
}

cat(sprintf("\n  Bootstraps with >=1 selected pathway: %d / %d\n",
            boot_selected_count, B))
cat(sprintf("  lambda.min fallbacks: %d / %d\n", boot_lambda_min_count, B))

# Frequency of each pathway being selected
selection_freq <- colMeans(boot_coef_matrix != 0)
freq_df <- data.frame(
  Pathway = names(selection_freq),
  SelectionFreq = round(selection_freq, 4),
  stringsAsFactors = FALSE
)
freq_df <- freq_df[order(freq_df$SelectionFreq, decreasing = TRUE), ]

cat("\n  Top 10 most selected pathways:\n")
for (i in 1:min(10, nrow(freq_df))) {
  cat(sprintf("    %2d. %s: %.1f%%\n", i, freq_df$Pathway[i],
              freq_df$SelectionFreq[i] * 100))
}

# Jaccard stability
compute_pairwise_jaccard <- function(coef_mat, n_sample = 200) {
  selected <- coef_mat != 0
  n_boot <- nrow(selected)
  if (n_boot < 2) return(NA)
  
  if (n_boot > n_sample) {
    idx <- sample(n_boot, n_sample)
    selected <- selected[idx, ]
    n_boot <- n_sample
  }
  
  jaccard_vals <- numeric(n_boot * (n_boot - 1) / 2)
  idx <- 1
  for (i in 1:(n_boot - 1)) {
    set_i <- which(selected[i, ])
    for (j in (i + 1):n_boot) {
      set_j <- which(selected[j, ])
      intersection <- length(intersect(set_i, set_j))
      union <- length(union(set_i, set_j))
      jaccard_vals[idx] <- ifelse(union > 0, intersection / union, 0)
      idx <- idx + 1
    }
  }
  return(mean(jaccard_vals, na.rm = TRUE))
}

# Only compute if there are enough non-zero bootstraps
valid_boot <- rowSums(boot_coef_matrix != 0) > 0
n_valid <- sum(valid_boot)

if (n_valid >= 10) {
  jaccard_stability <- compute_pairwise_jaccard(boot_coef_matrix[valid_boot, ])
  cat(sprintf("\n  Pairwise Jaccard stability: %.4f (based on %d valid bootstraps)\n",
              jaccard_stability, n_valid))
} else {
  jaccard_stability <- NA
  cat(sprintf("\n  Jaccard stability: NA (only %d valid bootstraps)\n", n_valid))
}

# Mean coefficient and sign consistency
mean_coef <- colMeans(boot_coef_matrix, na.rm = TRUE)
sign_consistency <- apply(boot_coef_matrix, 2, function(col) {
  if (all(col == 0)) return(NA)
  mean(sign(col[boot_coef_matrix[, which(colnames(boot_coef_matrix) == names(col))] != 0]) > 0, na.rm = TRUE)
})

# Fix sign consistency
sign_consistency <- apply(boot_coef_matrix, 2, function(col) {
  non_zero <- col != 0
  if (sum(non_zero) == 0) return(NA)
  pos_frac <- mean(col[non_zero] > 0, na.rm = TRUE)
  return(max(pos_frac, 1 - pos_frac))
})

# Mean |coefficient| among non-zero
mean_abs_coef_nonzero <- apply(boot_coef_matrix, 2, function(col) {
  nz <- col[col != 0]
  if (length(nz) == 0) return(NA)
  return(mean(abs(nz), na.rm = TRUE))
})

# Stability summary
stability_df <- data.frame(
  Pathway = colnames(boot_coef_matrix),
  SelectionFreq = round(selection_freq, 4),
  MeanCoef = round(mean_coef, 6),
  MeanAbsCoefNonzero = round(mean_abs_coef_nonzero, 6),
  SignConsistency = round(sign_consistency, 4),
  CV = apply(boot_coef_matrix, 2, function(col) {
    nz <- col[col != 0]
    if (length(nz) < 2) return(NA)
    sd(nz) / abs(mean(nz))
  }),
  stringsAsFactors = FALSE
)

stability_df <- stability_df[order(stability_df$SelectionFreq, decreasing = TRUE), ]
write.csv(stability_df, file.path(OUT_TAB_DIR, "merged_bootstrap_stability.csv"), row.names = FALSE)

# Stability report
sink(file.path(OUT_TAB_DIR, "merged_bootstrap_report.txt"))
cat("Bootstrap Stability Report (GSE19860 Merged Training)\n")
cat("===================================================\n\n")
cat(sprintf("Bootstraps: B = %d\n", B))
cat(sprintf("Valid bootstraps (>=1 nonzero): %d\n", n_valid))
cat(sprintf("lambda.min fallbacks: %d / %d (%.1f%%)\n\n",
            boot_lambda_min_count, B, boot_lambda_min_count / B * 100))
cat(sprintf("Pairwise Jaccard stability: %.4f\n", jaccard_stability))
cat(sprintf("Pathways with >50%% selection frequency: %d\n\n",
            sum(stability_df$SelectionFreq > 0.5, na.rm = TRUE)))

cat("Top 15 most selected pathways:\n")
for (i in 1:min(15, nrow(stability_df))) {
  cat(sprintf("  %2d. %s — Freq=%.1f%%, SignCons=%.1f%%, CV=%.2f\n",
              i, stability_df$Pathway[i],
              stability_df$SelectionFreq[i] * 100,
              stability_df$SignConsistency[i] * 100,
              stability_df$CV[i]))
}

cat("\n\nComparison with original training (GSE39582+GSE104645, Script 25):\n")
cat(sprintf("  Original Jaccard = 0.1002\n"))
cat(sprintf("  New Jaccard      = %.4f\n", jaccard_stability))
if (!is.na(jaccard_stability)) {
  delta <- jaccard_stability - 0.1002
  cat(sprintf("  Delta            = %+.4f\n", delta))
  if (delta > 0.10) {
    cat("  ✓ SIGNIFICANT IMPROVEMENT: Jaccard exceeded +0.10 threshold\n")
  } else if (delta > 0.05) {
    cat("  ✓ Moderate improvement\n")
  } else {
    cat("  Minimal change\n")
  }
}

cat("\n\nSelected Pathways in Final Model:\n")
if (nrow(coef_table) > 0) {
  for (i in seq_len(nrow(coef_table))) {
    pw <- coef_table$Pathway[i]
    freq <- stability_df$SelectionFreq[stability_df$Pathway == pw]
    cat(sprintf("  %s (coef=%.4f, freq=%.1f%%)\n",
                pw, coef_table$Coefficient[i], freq * 100))
  }
} else {
  cat("  (empty model)\n")
}
sink()

cat("\n  Stability report saved\n")

# Stability curve (cumulative Jaccard)
if (n_valid >= 20) {
  cat("\n  Generating cumulative stability curve...\n")
  selected_mat <- boot_coef_matrix[valid_boot, ] != 0
  n_valid_use <- min(n_valid, 500)
  
  cumulative_jaccard <- numeric(n_valid_use - 1)
  for (k in 2:n_valid_use) {
    set_a <- which(selected_mat[1, ])
    set_b <- which(colMeans(selected_mat[1:k, , drop = FALSE]) > 0.5)
    intersection <- length(intersect(set_a, set_b))
    union_len <- length(union(set_a, set_b))
    cumulative_jaccard[k - 1] <- ifelse(union_len > 0, intersection / union_len, 0)
  }
  
  tryCatch({
    pdf(file.path(OUT_FIG_DIR, "bootstrap_stability_curve_merged.pdf"),
        width = 7, height = 5.5)
    plot(2:n_valid_use, cumulative_jaccard, type = "l", lwd = 2,
         col = "#2166AC", xlab = "Bootstrap Replicates",
         ylab = "Cumulative Jaccard Index",
         main = "Stability Curve (Merged Training)",
         ylim = c(0, 1))
    abline(h = 0.5, lty = 2, col = "gray50")
    text(n_valid_use * 0.8, 0.52, "Jaccard=0.5", col = "gray50", cex = 0.8)
    dev.off()
    cat("  Stability curve saved\n")
  }, error = function(e) cat("  Stability curve failed:", e$message, "\n"))
}

# ============================================================
# 8. Final summary
# ============================================================
cat("\n============================================================\n")
cat("Script 27 Complete\n")
cat("============================================================\n")
cat(sprintf("  Training N: %d (GSE39582=%d + GSE19860=%d)\n",
            length(train_y), sum(batch == "GSE39582"), sum(batch == "GSE19860")))
cat(sprintf("  Training AUC: %.3f\n", train_auc))
cat(sprintf("  Pathways selected: %d\n", nrow(coef_table)))
cat(sprintf("  Jaccard stability: %.4f\n", jaccard_stability))
if (nrow(coef_table) > 0) {
  cat("  PRS pathways:", paste(coef_table$Pathway, collapse = ", "), "\n")
}
cat("\n  Output:\n")
cat(sprintf("    Tables: %s\n", OUT_TAB_DIR))
cat(sprintf("    Figures: %s\n", OUT_FIG_DIR))
cat("============================================================\n")
