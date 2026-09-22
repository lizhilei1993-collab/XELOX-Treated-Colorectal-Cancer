# ============================================================
# 14_ssgsea_pathway_scoring.R (v2 - GMT-based, no msigdbr)
# Phase I: ssGSEA Pathway Activity Scoring
# For: GSE104645 (FOLFOX, discovery) + GSE28702/GSE72970/GSE69657 (validation)
#      + GSE39582 (secondary validation, RFS endpoint)
# ============================================================

lib_path <- "C:/Rlibs"
.libPaths(c(lib_path, .libPaths()))
Sys.setenv(TMPDIR = "C:/temp", TMP = "C:/temp", TEMP = "C:/temp")

require(GSVA, lib.loc = lib_path, quietly = TRUE)
require(GSEABase, lib.loc = lib_path, quietly = TRUE)
require(Biobase, lib.loc = lib_path, quietly = TRUE)
require(data.table, lib.loc = lib_path, quietly = TRUE)

PROJECT_ROOT <- "/path/to/xelox_project"
DATA_GEO_DIR <- file.path(PROJECT_ROOT, "data", "geo")
DATA_PROC_DIR <- file.path(PROJECT_ROOT, "data", "processed")
RESULTS_TAB_DIR <- file.path(PROJECT_ROOT, "results", "tables")
RESULTS_FIG_DIR <- file.path(PROJECT_ROOT, "results", "figures")
OUT_DIR <- file.path(RESULTS_TAB_DIR, "pathway_activity")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Phase I: ssGSEA Pathway Activity Scoring ===\n\n")

# ============================================================
# Step 0: Read GMT files (replaces msigdbr)
# ============================================================
cat("[0] Reading MSigDB GMT files...\n")

read_gmt <- function(gmt_file) {
  lines <- readLines(gmt_file, warn = FALSE)
  gs_list <- list()
  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    gs_name <- parts[1]
    genes <- parts[-(1:2)]  # skip name + description
    genes <- genes[nchar(genes) > 0]
    gs_list[[gs_name]] <- genes
  }
  cat("  Loaded", length(gs_list), "gene sets from", basename(gmt_file), "\n")
  return(gs_list)
}

# Paths to downloaded GMT files
gmt_h <- "C:/temp/h.all.v2024.1.symbols.gmt"
gmt_kegg <- "C:/temp/c2.cp.kegg_legacy.v2024.1.symbols.gmt"

if (!file.exists(gmt_h) || !file.exists(gmt_kegg)) {
  stop("GMT files not found at C:/temp/. Download them first.")
}

hallmark_all <- read_gmt(gmt_h)
kegg_all <- read_gmt(gmt_kegg)

# ============================================================
# Step 1: Define oxaliplatin resistance-related gene sets
# ============================================================
cat("[1] Selecting oxaliplatin resistance pathway gene sets...\n")

target_pathways <- list()

# 1a. DNA Damage Repair (core oxaliplatin mechanism)
kegg_targets <- c(
  "KEGG_NUCLEOTIDE_EXCISION_REPAIR",    # NER - platinum crosslink repair
  "KEGG_BASE_EXCISION_REPAIR",           # BER - oxidative damage repair
  "KEGG_MISMATCH_REPAIR",               # MMR - platinum adduct repair  
  "KEGG_HOMOLOGOUS_RECOMBINATION",       # HR - DSB repair
  "KEGG_P53_SIGNALING_PATHWAY",         # p53 - apoptosis/DNA damage response
  "KEGG_APOPTOSIS",                     # Apoptosis
  "KEGG_CELL_CYCLE"                     # Cell cycle checkpoint
)

# 1b. Drug metabolism & transport
kegg_targets <- c(kegg_targets,
  "KEGG_ABC_TRANSPORTERS",              # Drug efflux
  "KEGG_GLUTATHIONE_METABOLISM",        # GST detoxification
  "KEGG_DRUG_METABOLISM_CYTOCHROME_P450",
  "KEGG_DRUG_METABOLISM_OTHER_ENZYMES",
  "KEGG_METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450"
)

# 1c. EMT, Metastasis, Stemness
kegg_targets <- c(kegg_targets,
  "KEGG_WNT_SIGNALING_PATHWAY",         # EMT driver
  "KEGG_TGF_BETA_SIGNALING_PATHWAY",    # EMT driver
  "KEGG_NOTCH_SIGNALING_PATHWAY",       # Stemness
  "KEGG_FOCAL_ADHESION",               # Metastasis
  "KEGG_ECM_RECEPTOR_INTERACTION"      # ECM remodeling
)

# 1d. Signaling pathways
kegg_targets <- c(kegg_targets,
  "KEGG_MAPK_SIGNALING_PATHWAY",
  "KEGG_MTOR_SIGNALING_PATHWAY",
  "KEGG_JAK_STAT_SIGNALING_PATHWAY"
)

# 1e. Immune / Inflammation
kegg_targets <- c(kegg_targets,
  "KEGG_CHEMOKINE_SIGNALING_PATHWAY",
  "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY"
)

# 1f. Additional cancer pathways
kegg_targets <- c(kegg_targets,
  "KEGG_COLORECTAL_CANCER",
  "KEGG_PATHWAYS_IN_CANCER"
)

# Match against available KEGG gene sets
for (id in kegg_targets) {
  if (id %in% names(kegg_all)) {
    target_pathways[[id]] <- kegg_all[[id]]
    cat("  +", id, "-", length(kegg_all[[id]]), "genes\n")
  } else {
    cat("  -", id, "NOT FOUND in KEGG Legacy, skipping\n")
  }
}

# Also add Hallmark gene sets that are relevant
hallmark_targets <- c(
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_HYPOXIA",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_TGF_BETA_SIGNALING"
)

for (id in hallmark_targets) {
  if (id %in% names(hallmark_all)) {
    target_pathways[[id]] <- hallmark_all[[id]]
    cat("  +", id, "-", length(hallmark_all[[id]]), "genes\n")
  } else {
    cat("  -", id, "NOT FOUND in Hallmark, skipping\n")
  }
}

cat("\n  Total target gene sets:", length(target_pathways), "\n")
sizes <- sapply(target_pathways, length)
cat("  Size range:", min(sizes), "-", max(sizes), "genes (median:", median(sizes), ")\n")

# ============================================================
# Step 2: Build gene-level expression for each dataset
# ============================================================
cat("\n[2] Building gene-level expression matrices...\n")

get_gene_expression <- function(expr_matrix, probe_map, probe_col, gene_col) {
  if (is.null(rownames(expr_matrix))) {
    stop("Expression matrix must have rownames (probe IDs)")
  }
  
  probe_ids <- rownames(expr_matrix)
  idx <- match(probe_ids, probe_map[[probe_col]])
  matched <- probe_map[idx, ]
  
  valid <- !is.na(idx)
  expr_matched <- expr_matrix[valid, , drop = FALSE]
  genes <- matched[[gene_col]][valid]
  
  has_gene <- !is.na(genes) & genes != ""
  expr_genes <- expr_matched[has_gene, , drop = FALSE]
  genes <- genes[has_gene]
  
  uniq_genes <- unique(genes)
  expr_gene <- matrix(NA, nrow = length(uniq_genes), ncol = ncol(expr_matrix))
  rownames(expr_gene) <- uniq_genes
  colnames(expr_gene) <- colnames(expr_matrix)
  
  for (i in seq_along(uniq_genes)) {
    g <- uniq_genes[i]
    probe_set <- which(genes == g)
    if (length(probe_set) == 1) {
      expr_gene[i, ] <- expr_genes[probe_set, ]
    } else {
      expr_gene[i, ] <- colMeans(expr_genes[probe_set, , drop = FALSE])
    }
  }
  
  return(expr_gene)
}

# --- Load GPL570 probe map ---
cat("  Loading GPL570 probe-to-gene map...\n")
probe_map_570 <- read.csv(file.path(DATA_GEO_DIR, "GPL570_probe_gene_map.csv"),
                           stringsAsFactors = FALSE)
cat("  GPL570 probe map:", nrow(probe_map_570), "probes\n")

# --- Process GSE39582 (GPL570, ExpressionSet) ---
cat("\n  Processing GSE39582...\n")
eset <- readRDS(file.path(DATA_GEO_DIR, "GSE39582_eset.rds"))
expr_39582_probe <- exprs(eset)
cat("  Probe-level dim:", nrow(expr_39582_probe), "x", ncol(expr_39582_probe), "\n")
expr_39582 <- get_gene_expression(expr_39582_probe, probe_map_570, "probe_id", "gene_symbol")
cat("  Gene-level dim:", nrow(expr_39582), "x", ncol(expr_39582), "\n")
saveRDS(expr_39582, file.path(OUT_DIR, "GSE39582_gene_expression.rds"))

# --- Process GSE28702 (GPL570, matrix) ---
cat("\n  Processing GSE28702...\n")
expr_28702_probe <- readRDS(file.path(DATA_GEO_DIR, "GSE28702_expression.rds"))
cat("  Probe-level dim:", nrow(expr_28702_probe), "x", ncol(expr_28702_probe), "\n")
expr_28702 <- get_gene_expression(expr_28702_probe, probe_map_570, "probe_id", "gene_symbol")
cat("  Gene-level dim:", nrow(expr_28702), "x", ncol(expr_28702), "\n")
saveRDS(expr_28702, file.path(OUT_DIR, "GSE28702_gene_expression.rds"))

# --- Process GSE72970 (GPL570, matrix) ---
cat("\n  Processing GSE72970...\n")
expr_72970_probe <- readRDS(file.path(DATA_GEO_DIR, "GSE72970_expression.rds"))
cat("  Probe-level dim:", nrow(expr_72970_probe), "x", ncol(expr_72970_probe), "\n")
expr_72970 <- get_gene_expression(expr_72970_probe, probe_map_570, "probe_id", "gene_symbol")
cat("  Gene-level dim:", nrow(expr_72970), "x", ncol(expr_72970), "\n")
saveRDS(expr_72970, file.path(OUT_DIR, "GSE72970_gene_expression.rds"))

# --- Process GSE69657 (GPL570, matrix) ---
cat("\n  Processing GSE69657...\n")
expr_69657_probe <- readRDS(file.path(DATA_GEO_DIR, "GSE69657_expression.rds"))
cat("  Probe-level dim:", nrow(expr_69657_probe), "x", ncol(expr_69657_probe), "\n")
expr_69657 <- get_gene_expression(expr_69657_probe, probe_map_570, "probe_id", "gene_symbol")
cat("  Gene-level dim:", nrow(expr_69657), "x", ncol(expr_69657), "\n")
saveRDS(expr_69657, file.path(OUT_DIR, "GSE69657_gene_expression.rds"))

# --- Process GSE104645 (GPL6480, Agilent, matrix) ---
cat("\n  Processing GSE104645 (GPL6480)...\n")
expr_104645 <- readRDS(file.path(DATA_GEO_DIR, "GSE104645_expression.rds"))
rns <- rownames(expr_104645)
non_ctrl <- !grepl("^\\(\\+\\)", rns)
expr_104645_nonctrl <- expr_104645[non_ctrl, ]
cat("  Agilent probes (non-control):", nrow(expr_104645_nonctrl), "\n")

# Build GPL6480 probe-to-gene map
cat("  Building GPL6480 probe-to-gene map...\n")
conn <- gzfile(file.path(DATA_GEO_DIR, "GPL6480.annot.gz"), "rt")
while (TRUE) {
  line <- readLines(conn, n = 1)
  if (length(line) == 0 || grepl("!platform_table_begin", line)) break
}
col_names <- strsplit(readLines(conn, n = 1), "\t")[[1]]
probe_map_6480 <- read.table(conn, sep = "\t", header = FALSE, 
  col.names = col_names, fill = TRUE, quote = "\"", 
  stringsAsFactors = FALSE)
close(conn)
cat("  GPL6480 annot loaded:", nrow(probe_map_6480), "rows\n")

# Filter to only Agilent probes in expression data
probe_map_6480 <- probe_map_6480[probe_map_6480$ID %in% rownames(expr_104645_nonctrl), ]
cat("  Matching probes in expr:", nrow(probe_map_6480), "\n")

# Aggregate to gene level
probe_genes <- probe_map_6480$Gene.symbol
has_gene <- !is.na(probe_genes) & probe_genes != ""
probe_map_6480 <- probe_map_6480[has_gene, ]
probe_genes <- probe_genes[has_gene]
cat("  Probes with gene symbols:", length(probe_genes), "\n")

uniq_genes <- unique(probe_genes)
gene_expr <- matrix(NA, nrow = length(uniq_genes), ncol = ncol(expr_104645_nonctrl))
rownames(gene_expr) <- uniq_genes
colnames(gene_expr) <- colnames(expr_104645_nonctrl)

row_idx <- match(probe_map_6480$ID, rownames(expr_104645_nonctrl))
for (i in seq_along(uniq_genes)) {
  g <- uniq_genes[i]
  probe_set <- which(probe_genes == g)
  if (length(probe_set) == 1) {
    gene_expr[i, ] <- expr_104645_nonctrl[row_idx[probe_set], ]
  } else {
    gene_expr[i, ] <- colMeans(expr_104645_nonctrl[row_idx[probe_set], , drop = FALSE])
  }
}

expr_104645 <- gene_expr
cat("  Gene-level dim:", nrow(expr_104645), "x", ncol(expr_104645), "\n")
saveRDS(expr_104645, file.path(OUT_DIR, "GSE104645_gene_expression.rds"))

# ============================================================
# Step 3: Filter gene sets to common genes
# ============================================================
cat("\n[3] Filtering gene sets to genes present across datasets...\n")

all_datasets <- list(
  GSE39582 = expr_39582,
  GSE28702 = expr_28702,
  GSE72970 = expr_72970,
  GSE69657 = expr_69657,
  GSE104645 = expr_104645
)

common_genes <- Reduce(intersect, lapply(all_datasets, rownames))
cat("  Genes common across all 5 datasets:", length(common_genes), "\n")

# Filter gene sets
filtered_pathways <- list()
for (nm in names(target_pathways)) {
  gs <- target_pathways[[nm]]
  overlapping <- intersect(gs, common_genes)
  if (length(overlapping) >= 10) {
    filtered_pathways[[nm]] <- overlapping
  }
}
cat("  Usable gene sets (>=10 common genes):", length(filtered_pathways), "\n")
cat("  Removed (too few genes):", length(target_pathways) - length(filtered_pathways), "\n")

# ============================================================
# Step 4: Run ssGSEA for each dataset
# ============================================================
cat("\n[4] Running ssGSEA pathway activity scoring...\n")

run_ssgsea <- function(expr_mat, gene_sets, dataset_name) {
  cat("  Computing ssGSEA for", dataset_name, "...\n")
  
  if (!is.matrix(expr_mat)) expr_mat <- as.matrix(expr_mat)
  
  common <- intersect(rownames(expr_mat), unique(unlist(gene_sets)))
  cat("    Genes in gene sets present in expression:", length(common), "/", 
      length(unique(unlist(gene_sets))), "\n")
  
  # GSVA v2.x API: use ssgseaParam object instead of direct matrix call
  param <- ssgseaParam(expr_mat, gene_sets,
                       minSize = 10, maxSize = 500,
                       alpha = 0.25, normalize = TRUE)
  result <- gsva(param)
  
  cat("    Pathway activity matrix:", nrow(result), "x", ncol(result), "\n")
  return(result)
}

pathway_scores <- list()
for (nm in names(all_datasets)) {
  pathway_scores[[nm]] <- run_ssgsea(all_datasets[[nm]], filtered_pathways, nm)
  saveRDS(pathway_scores[[nm]], file.path(OUT_DIR, paste0(nm, "_pathway_scores.rds")))
}

# ============================================================
# Step 5: Summary statistics
# ============================================================
cat("\n[5] Pathway activity summary...\n")

summary_df <- data.frame(
  Dataset = names(pathway_scores),
  Samples = sapply(pathway_scores, ncol),
  Pathways = sapply(pathway_scores, nrow),
  stringsAsFactors = FALSE
)
print(summary_df)
write.csv(summary_df, file.path(OUT_DIR, "pathway_scoring_summary.csv"), row.names = FALSE)

cat("\n=== ssGSEA Pathway Activity Scoring Complete ===\n")
cat("Output directory:", OUT_DIR, "\n")
cat("Files:\n")
cat("  *_gene_expression.rds - gene-level expression matrices\n")
cat("  *_pathway_scores.rds - ssGSEA pathway activity scores\n")
cat("  pathway_scoring_summary.csv - summary table\n")
