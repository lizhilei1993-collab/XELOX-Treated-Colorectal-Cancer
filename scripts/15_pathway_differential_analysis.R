# ============================================================
# 15_PATHWAY_DIFFERENTIAL_ANALYSIS.R
# Phase II: Differential Pathway Activity Analysis
#
# Purpose:
#   Compare ssGSEA pathway activity between Resistant and Sensitive
#   groups across 5 GEO cohorts. Find consensus resistance pathways.
#
# Input:
#   results/tables/pathway_activity/*_pathway_scores.rds
#   results/tables/GSE*_clinical_data.csv / GSE39582_xelox_groups.csv
#
# Output:
#   tables: differential_pathway_summary.csv, *_pathway_diff.csv,
#           pathway_consistency_matrix.csv, top_consensus_pathways.csv
#   figures: pathway_volcano_5panel.pdf, pathway_heatmap_consensus.pdf,
#            pathway_boxplot_top4.pdf
# ============================================================

Sys.setenv(TMPDIR = "/tmp", TMP = "/tmp", TEMP = "/tmp")
.libPaths(c("/path/to/Rlibs", .libPaths()))

# PROJECT_ROOT from environment variable to avoid Chinese path issues
# Set via PowerShell: [Environment]::SetEnvironmentVariable("XELOX_ROOT","/path/to/xelox_project","User")
PROJECT_ROOT <- Sys.getenv("XELOX_ROOT",
  "/path/to/xelox_project")
DATA_GEO_DIR     <- file.path(PROJECT_ROOT, "data", "geo")
RESULTS_TAB_DIR  <- file.path(PROJECT_ROOT, "results", "tables")
RESULTS_FIG_DIR  <- file.path(PROJECT_ROOT, "results", "figures")
PA_DIR           <- file.path(RESULTS_TAB_DIR, "pathway_activity")

suppressWarnings(dir.create(PA_DIR, showWarnings = FALSE, recursive = TRUE))

cat("=== Phase II: Differential Pathway Activity Analysis ===\n\n")

# ============================================================
# Dependencies
# ============================================================
cat("[0] Loading packages...\n")
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(gridExtra)
library(grid)

set.seed(42)

# Output subdirectory for figures
FIG_DIR <- file.path(RESULTS_FIG_DIR, "pathway_diff")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Step 1: Sample annotation functions
# ============================================================
cat("\n[1] Defining sample annotation functions...\n")

annotate_GSE104645 <- function(ps, clin_file) {
  clin <- read.csv(clin_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  
  regimen_col <- "1st-line chemotherapy regimens:ch1"
  response_col <- "best response of 1st-line chemotherapy:ch1"
  
  # Filter to XELOX/FOLFOX/SOX regimens
  is_oxali <- grepl("XELOX|FOLFOX|SOX|CAPOX", clin[[regimen_col]], ignore.case = TRUE)
  oxali_samples <- rownames(clin)[is_oxali]
  
  # Map response
  resp <- clin[[response_col]]
  names(resp) <- rownames(clin)
  
  group <- rep(NA, length(oxali_samples))
  names(group) <- oxali_samples
  group[resp[oxali_samples] %in% c("Complete response", "Partial response")] <- "Sensitive"
  group[resp[oxali_samples] %in% c("Progressive disease", "Stable disease")] <- "Resistant"
  
  # Match to pathway score samples
  common <- intersect(names(group), colnames(ps))
  group <- group[common]
  group <- group[!is.na(group)]
  
  cat("    GSE104645: Sensitive=", sum(group == "Sensitive"),
      " Resistant=", sum(group == "Resistant"),
      " (oxaliplatin-based,", length(group), "total)\n")
  
  list(pathway_matrix = ps[, names(group), drop = FALSE],
       group = group)
}

annotate_GSE28702 <- function(ps, clin_file) {
  clin <- read.csv(clin_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  
  resp_col <- "mfolfox6:ch1"
  resp <- clin[[resp_col]]
  names(resp) <- rownames(clin)
  
  group <- rep(NA, length(resp))
  names(group) <- rownames(clin)
  group[resp == "responder"] <- "Sensitive"
  group[resp == "non-responder"] <- "Resistant"
  
  common <- intersect(names(group), colnames(ps))
  group <- group[common]
  group <- group[!is.na(group)]
  
  cat("    GSE28702: Sensitive=", sum(group == "Sensitive"),
      " Resistant=", sum(group == "Resistant"),
      " (pure mFOLFOX6,", length(group), "total)\n")
  
  list(pathway_matrix = ps[, names(group), drop = FALSE],
       group = group)
}

annotate_GSE72970 <- function(ps, clin_file) {
  clin <- read.csv(clin_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  
  resp_col <- "response category:ch1"
  resp <- clin[[resp_col]]
  names(resp) <- rownames(clin)
  
  group <- rep(NA, length(resp))
  names(group) <- rownames(clin)
  group[resp %in% c("CR", "PR")] <- "Sensitive"
  group[resp %in% c("PD", "SD")] <- "Resistant"
  
  common <- intersect(names(group), colnames(ps))
  group <- group[common]
  group <- group[!is.na(group)]
  
  cat("    GSE72970: Sensitive=", sum(group == "Sensitive"),
      " Resistant=", sum(group == "Resistant"),
      " (all regimens,", length(group), "total)\n")
  
  list(pathway_matrix = ps[, names(group), drop = FALSE],
       group = group)
}

annotate_GSE69657 <- function(ps, clin_file) {
  clin <- read.csv(clin_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  
  resp_col <- "chemoresponse:ch1"
  resp <- clin[[resp_col]]
  names(resp) <- rownames(clin)
  
  group <- rep(NA, length(resp))
  names(group) <- rownames(clin)
  group[resp == "responder"] <- "Sensitive"
  group[resp == "noresponder"] <- "Resistant"
  
  common <- intersect(names(group), colnames(ps))
  group <- group[common]
  group <- group[!is.na(group)]
  
  cat("    GSE69657: Sensitive=", sum(group == "Sensitive"),
      " Resistant=", sum(group == "Resistant"),
      " (pure XELOX neoadjuvant,", length(group), "total)\n")
  
  list(pathway_matrix = ps[, names(group), drop = FALSE],
       group = group)
}

annotate_GSE39582 <- function(ps, group_file) {
  clin <- read.csv(group_file, stringsAsFactors = FALSE)
  
  # GSE39582_xelox_groups has sample_id, group, etc.
  group_raw <- clin$group
  names(group_raw) <- clin$sample_id
  
  # Normalize case (file uses lowercase: "resistant"/"sensitive")
  group_raw <- tolower(group_raw)
  
  common <- intersect(names(group_raw), colnames(ps))
  group <- group_raw[common]
  
  # Only keep resistant/sensitive (drop intermediate), convert to title case
  keep <- group %in% c("resistant", "sensitive")
  group <- group[keep]
  
  # downstream uses "Resistant" / "Sensitive"
  group <- ifelse(group == "resistant", "Resistant", "Sensitive")
  
  cat("    GSE39582: Sensitive=", sum(group == "Sensitive"),
      " Resistant=", sum(group == "Resistant"),
      " (RFS-based,", length(group), "used, intermediate excluded)\n")
  
  list(pathway_matrix = ps[, names(group), drop = FALSE],
       group = group)
}

# ============================================================
# Step 2: Differential analysis function
# ============================================================
cat("[2] Defining differential analysis function...\n")

run_differential_analysis <- function(pathway_matrix, group, dataset_name) {
  pathways <- rownames(pathway_matrix)
  
  sens_idx <- which(group == "Sensitive")
  res_idx <- which(group == "Resistant")
  
  sens_mat <- pathway_matrix[, sens_idx, drop = FALSE]
  res_mat <- pathway_matrix[, res_idx, drop = FALSE]
  
  results <- data.frame(
    pathway = pathways,
    dataset = dataset_name,
    n_sensitive = length(sens_idx),
    n_resistant = length(res_idx),
    median_sensitive = NA_real_,
    median_resistant = NA_real_,
    delta = NA_real_,
    p_value = NA_real_,
    fdr = NA_real_,
    mean_rank_sensitive = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(pathways)) {
    pw <- pathways[i]
    s_vals <- as.numeric(sens_mat[pw, ])
    r_vals <- as.numeric(res_mat[pw, ])
    
    med_s <- median(s_vals, na.rm = TRUE)
    med_r <- median(r_vals, na.rm = TRUE)
    
    # delta: Resistant - Sensitive (positive = up in resistant)
    delta <- med_r - med_s
    
    # Wilcoxon rank-sum test
    wt <- wilcox.test(r_vals, s_vals, alternative = "two.sided", exact = FALSE)
    pv <- wt$p.value
    
    # Mean rank for info
    all_vals <- c(r_vals, s_vals)
    ranks <- rank(all_vals)
    mean_rank_r <- mean(ranks[1:length(r_vals)])
    
    results$median_sensitive[i] <- med_s
    results$median_resistant[i] <- med_r
    results$delta[i] <- delta
    results$p_value[i] <- pv
    results$mean_rank_sensitive[i] <- mean_rank_r
  }
  
  # BH correction
  results$fdr <- p.adjust(results$p_value, method = "BH")
  results$significant <- results$fdr < 0.05
  
  # Sort by p-value
  results <- results[order(results$p_value), ]
  rownames(results) <- NULL
  
  cat("    Significant pathways (FDR<0.05):", sum(results$significant), "/", nrow(results), "\n")
  
  return(results)
}

# ============================================================
# Step 3: Main pipeline - Process all 5 datasets
# ============================================================
cat("\n[3] Running differential analysis across all datasets...\n")

dataset_config <- list(
  GSE104645 = list(
    ps_file = file.path(PA_DIR, "GSE104645_pathway_scores.rds"),
    clin_file = file.path(RESULTS_TAB_DIR, "GSE104645_clinical_data.csv"),
    annotator = "GSE104645"
  ),
  GSE28702 = list(
    ps_file = file.path(PA_DIR, "GSE28702_pathway_scores.rds"),
    clin_file = file.path(RESULTS_TAB_DIR, "GSE28702_clinical_data.csv"),
    annotator = "GSE28702"
  ),
  GSE72970 = list(
    ps_file = file.path(PA_DIR, "GSE72970_pathway_scores.rds"),
    clin_file = file.path(RESULTS_TAB_DIR, "GSE72970_clinical_data.csv"),
    annotator = "GSE72970"
  ),
  GSE69657 = list(
    ps_file = file.path(PA_DIR, "GSE69657_pathway_scores.rds"),
    clin_file = file.path(RESULTS_TAB_DIR, "GSE69657_clinical_data.csv"),
    annotator = "GSE69657"
  ),
  GSE39582 = list(
    ps_file = file.path(PA_DIR, "GSE39582_pathway_scores.rds"),
    clin_file = file.path(RESULTS_TAB_DIR, "GSE39582_xelox_groups.csv"),
    annotator = "GSE39582"
  )
)

annotator_map <- list(
  GSE104645 = annotate_GSE104645,
  GSE28702 = annotate_GSE28702,
  GSE72970 = annotate_GSE72970,
  GSE69657 = annotate_GSE69657,
  GSE39582 = annotate_GSE39582
)

all_results <- list()
all_delta <- list()
all_pval <- list()

for (ds_name in names(dataset_config)) {
  cfg <- dataset_config[[ds_name]]
  
  cat("\n  >>> Processing", ds_name, "...\n")
  
  # Load pathway scores
  ps <- readRDS(cfg$ps_file)
  cat("    Pathway scores dim:", nrow(ps), "x", ncol(ps), "\n")
  
  # Annotate samples
  annotated <- annotator_map[[ds_name]](ps, cfg$clin_file)
  
  # Run differential analysis
  res <- run_differential_analysis(annotated$pathway_matrix, annotated$group, ds_name)
  
  # Save per-dataset results
  write.csv(res, file.path(PA_DIR, paste0(ds_name, "_pathway_diff.csv")), row.names = FALSE)
  
  all_results[[ds_name]] <- res
  all_delta[[ds_name]] <- res$delta
  all_pval[[ds_name]] <- res$p_value
  names(all_delta[[ds_name]]) <- res$pathway
  names(all_pval[[ds_name]]) <- res$pathway
}

# Build consensus matrix
all_pathways <- unique(unlist(lapply(all_results, `[[`, "pathway")))

cat("\n[4] Cross-cohort consistency summary...\n")

consensus_list <- list()
for (pw in all_pathways) {
  deltas <- sapply(all_delta, function(x) ifelse(pw %in% names(x), x[pw], NA))
  pvals <- sapply(all_pval, function(x) ifelse(pw %in% names(x), x[pw], NA))
  
  n_available <- sum(!is.na(deltas))
  n_up <- sum(deltas > 0, na.rm = TRUE)
  n_down <- sum(deltas < 0, na.rm = TRUE)
  direction_consensus <- if (n_up >= n_down) n_up else n_down
  direction_sign <- if (n_up >= n_down) "Up_in_Resistant" else "Down_in_Resistant"
  
  # Voting: direction consistent + FDR < 0.05
  n_sig_up <- sum(deltas > 0 & pvals < 0.05, na.rm = TRUE)
  n_sig_down <- sum(deltas < 0 & pvals < 0.05, na.rm = TRUE)
  n_sig_consistent <- if (n_up >= n_down) n_sig_up else n_sig_down
  
  # Mean rank by abs(delta)
  abs_deltas <- abs(deltas)
  ranks <- rank(-abs_deltas, na.last = "keep", ties.method = "average")
  mean_rank <- mean(ranks, na.rm = TRUE)
  
  # Mean delta
  mean_delta <- mean(deltas, na.rm = TRUE)
  
  consensus_list[[pw]] <- data.frame(
    pathway = pw,
    n_cohorts = n_available,
    direction_consensus = paste0(direction_consensus, "/", n_available),
    direction_sign = direction_sign,
    n_sig_consistent = n_sig_consistent,
    mean_delta = round(mean_delta, 4),
    mean_rank = round(mean_rank, 1),
    stringsAsFactors = FALSE
  )
  
  # Add per-cohort delta
  for (ds in names(all_delta)) {
    val <- deltas[ds]
    consensus_list[[pw]][[paste0("delta_", ds)]] <- round(val, 4)
  }
  for (ds in names(all_pval)) {
    val <- pvals[ds]
    consensus_list[[pw]][[paste0("pval_", ds)]] <- ifelse(is.na(val), NA, format(val, scientific = TRUE, digits = 3))
    consensus_list[[pw]][[paste0("fdr_", ds)]] <- ifelse(is.na(val), NA,
      format(p.adjust(val, method = "BH"), scientific = TRUE, digits = 3))
  }
}

consensus_df <- do.call(rbind, consensus_list)
consensus_df <- consensus_df[order(consensus_df$mean_rank), ]
row.names(consensus_df) <- NULL

write.csv(consensus_df, file.path(PA_DIR, "pathway_consistency_matrix.csv"), row.names = FALSE)
cat("  Consistency matrix written:", nrow(consensus_df), "pathways\n")

# Top consensus pathways
top_consensus <- consensus_df[
  as.numeric(sub("/.*", "", consensus_df$direction_consensus)) >= 4 &
  consensus_df$n_sig_consistent >= 2, 
]
write.csv(top_consensus, file.path(PA_DIR, "top_consensus_pathways.csv"), row.names = FALSE)
cat("  Top consensus pathways (>=4/5 direction + >=2 sig):", nrow(top_consensus), "\n")

# Summary table
summary_rows <- list()
for (ds_name in names(all_results)) {
  r <- all_results[[ds_name]]
  summary_rows[[ds_name]] <- data.frame(
    Dataset = ds_name,
    Sensitive = r$n_sensitive[1],
    Resistant = r$n_resistant[1],
    n_pathways = nrow(r),
    sig_fdr05 = sum(r$significant),
    top_pathway = r$pathway[1],
    top_delta = round(r$delta[1], 4),
    top_pval = format(r$p_value[1], scientific = TRUE, digits = 3),
    stringsAsFactors = FALSE
  )
}
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(PA_DIR, "differential_pathway_summary.csv"), row.names = FALSE)
cat("\nSummary:\n")
print(summary_df[, c("Dataset", "Sensitive", "Resistant", "sig_fdr05", "top_pathway")])

# ============================================================
# Step 5: Visualization
# ============================================================
cat("\n[5] Generating figures...\n")

dataset_colors <- c(
  GSE104645 = "#E41A1C",
  GSE28702  = "#377EB8",
  GSE72970  = "#4DAF4A",
  GSE69657  = "#984EA3",
  GSE39582  = "#FF7F00"
)

# --- 5a. Volcano plots (5-panel) ---
cat("  5a. Volcano plot panel...\n")

plot_volcano <- function(res, ds_name, max_label = 5) {
  res$neg_log10p <- -log10(res$p_value)
  res$label <- ifelse(rank(res$p_value, ties.method = "first") <= max_label |
                      (res$significant & abs(res$delta) > 0.03), res$pathway, "")
  # Shorten pathway names for display
  res$label <- gsub("^KEGG_|^HALLMARK_", "", res$label)
  
  # Color by significance and direction
  res$color <- "NS"
  res$color[res$delta > 0 & res$p_value < 0.05] <- "Up_in_Resistant"
  res$color[res$delta < 0 & res$p_value < 0.05] <- "Down_in_Resistant"
  # FDR filter
  res$color_fdr <- res$color
  res$color_fdr[res$p_value >= 0.05] <- "NS"
  res$color_fdr[res$significant] <- res$color[res$significant]
  col_map <- c("Up_in_Resistant" = "#D73027", "Down_in_Resistant" = "#4575B4", "NS" = "grey60")
  
  max_p <- max(res$neg_log10p, na.rm = TRUE)
  max_d <- max(abs(res$delta), na.rm = TRUE)
  
  p <- ggplot(res, aes(x = delta, y = neg_log10p)) +
    geom_point(aes(color = color_fdr), alpha = 0.7, size = 2.5) +
    scale_color_manual(values = col_map) +
    geom_text_repel(aes(label = label), size = 2.8, max.overlaps = 12,
                    segment.color = "grey50", segment.size = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", alpha = 0.5) +
    labs(title = ds_name, x = "Delta (Resistant - Sensitive)", y = "-log10(p)") +
    xlim(-max_d * 1.1, max_d * 1.1) +
    ylim(0, max_p * 1.1) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  return(p)
}

volcano_plots <- list()
for (ds_name in names(all_results)) {
  volcano_plots[[ds_name]] <- plot_volcano(all_results[[ds_name]], ds_name, max_label = 4)
}

pdf(file.path(FIG_DIR, "pathway_volcano_5panel.pdf"), width = 14, height = 10)
grid.arrange(
  volcano_plots[[1]], volcano_plots[[2]],
  volcano_plots[[3]], volcano_plots[[4]],
  volcano_plots[[5]],
  ncol = 3, nrow = 2,
  top = textGrob("Differential Pathway Activity: Resistant vs Sensitive",
                 gp = gpar(fontsize = 14, fontface = "bold"))
)
dev.off()
cat("    Saved volcano panel\n")

# --- 5b. Heatmap of cross-cohort delta ---
cat("  5b. Heatmap...\n")

# Build a delta matrix (pathways x datasets)
delta_mat <- matrix(NA, nrow = length(all_pathways), ncol = length(all_results))
rownames(delta_mat) <- all_pathways
colnames(delta_mat) <- names(all_results)
for (ds_name in names(all_results)) {
  r <- all_results[[ds_name]]
  delta_mat[r$pathway, ds_name] <- r$delta
}

# Filter to pathways with at least 3 non-NA deltas
keep <- rowSums(!is.na(delta_mat)) >= 3
delta_mat <- delta_mat[keep, ]

# Order by mean absolute delta
delta_mat <- delta_mat[order(-rowMeans(abs(delta_mat), na.rm = TRUE)), ]

# Save original names for boxplot matching
pathway_original_names <- rownames(delta_mat)

# Simplify pathway names for display (append source prefix to avoid duplicates)
rownames_short <- gsub("^KEGG_|^HALLMARK_", "", rownames(delta_mat))
# Prefix with K: or H: to disambiguate when KEGG and HALLMARK share a name
prefix <- ifelse(grepl("^KEGG_", rownames(delta_mat)), "K:",
          ifelse(grepl("^HALLMARK_", rownames(delta_mat)), "H:", ""))
rownames_short <- paste0(prefix, rownames_short)
# Clean up signaling pathway suffix
rownames_short <- gsub("_SIGNALING_PATHWAY|_SIGNALING", "", rownames_short)
rownames_short <- gsub("_", " ", rownames_short)

# Build short→original name mapping for boxplot
short_to_orig <- setNames(pathway_original_names, rownames_short)

# Annotation: direction consensus
direction_labels <- rep("", nrow(delta_mat))
for (i in seq_len(nrow(delta_mat))) {
  pw <- rownames(delta_mat)[i]
  if (pw %in% consensus_df$pathway) {
    direction_labels[i] <- consensus_df$direction_sign[consensus_df$pathway == pw]
  }
}
annot_row <- data.frame(Direction = direction_labels, row.names = rownames(delta_mat))
annot_row$Direction <- gsub("Up_in_Resistant", "Up in Resistant", annot_row$Direction)
annot_row$Direction <- gsub("Down_in_Resistant", "Down in Resistant", annot_row$Direction)

# Truncate for display
row.names(delta_mat) <- rownames_short
rownames(annot_row) <- rownames_short

# Limit rows for readability (top 25)
n_show <- min(nrow(delta_mat), 25)
delta_mat_show <- delta_mat[1:n_show, ]
annot_row_show <- annot_row[1:n_show, , drop = FALSE]

# Color breaks
max_abs <- max(abs(delta_mat_show), na.rm = TRUE)

pdf(file.path(FIG_DIR, "pathway_heatmap_consensus.pdf"), width = 9, height = 8)
pheatmap(delta_mat_show,
         color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
         breaks = seq(-max_abs, max_abs, length.out = 101),
         cluster_rows = TRUE, cluster_cols = TRUE,
         annotation_row = annot_row_show,
         main = "Pathway Activity Delta (Resistant - Sensitive)\nAcross Cohorts",
         fontsize_row = 8, fontsize_col = 10,
         na_col = "grey90",
         display_numbers = TRUE, number_format = "%.3f", fontsize_number = 6)
dev.off()
cat("    Saved heatmap\n")

# --- 5c. Boxplot for top 4 consensus pathways ---
cat("  5c. Boxplot for top pathways...\n")

# Get top 4 consensus pathways (use ORIGINAL names for matrix matching)
if (nrow(top_consensus) > 0) {
  top4_orig <- head(top_consensus$pathway, min(4, nrow(top_consensus)))
  top4_display <- head(gsub("^KEGG_|^HALLMARK_", "", top_consensus$pathway), min(4, nrow(top_consensus)))
} else {
  # Fall back to pathways with highest mean |delta| (short names from delta_mat)
  top4_short <- head(names(sort(apply(abs(delta_mat), 1, mean, na.rm = TRUE), decreasing = TRUE)), 4)
  # Map back to original names using short_to_orig lookup
  top4_orig <- short_to_orig[top4_short]
  # Remove any NAs (shouldn't happen, but be safe)
  top4_orig <- top4_orig[!is.na(top4_orig)]
  top4_display <- top4_short
}

cat("    Top pathways for boxplot:", paste(top4_orig, collapse = ", "), "\n")

# Gather data for boxplot
box_data <- data.frame()
for (ds_name in names(all_results)) {
  ps <- readRDS(dataset_config[[ds_name]]$ps_file)
  annotated <- annotator_map[[ds_name]](ps, dataset_config[[ds_name]]$clin_file)
  mat <- annotated$pathway_matrix
  grp <- annotated$group
  
  for (pw in top4_orig) {
    if (pw %in% rownames(mat)) {
      vals <- as.numeric(mat[pw, ])
      box_data <- rbind(box_data, data.frame(
        pathway = gsub("^KEGG_|^HALLMARK_", "", pw),
        dataset = ds_name,
        group = grp,
        activity = vals,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Factor ordering
box_data$pathway <- factor(box_data$pathway, levels = rev(gsub("^KEGG_|^HALLMARK_", "", top4_orig)))
box_data$group <- factor(box_data$group, levels = c("Sensitive", "Resistant"))

p_box <- ggplot(box_data, aes(x = dataset, y = activity, fill = group)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3, alpha = 0.7) +
  facet_wrap(~ pathway, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = c("Sensitive" = "#4575B4", "Resistant" = "#D73027")) +
  labs(title = "Top Consensus Pathway Activity: Sensitive vs Resistant",
       x = "Cohort", y = "ssGSEA Score") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom")

ggsave(file.path(FIG_DIR, "pathway_boxplot_top4.pdf"), p_box, width = 10, height = 8)
cat("    Saved boxplot\n")

# ============================================================
# Done
# ============================================================
cat("\n=== Phase II: Differential Pathway Activity Analysis Complete ===\n")
cat("Output tables:", PA_DIR, "\n")
cat("Output figures:", FIG_DIR, "\n")
cat("Files:\n")
cat("  *_pathway_diff.csv - per-cohort differential analysis\n")
cat("  differential_pathway_summary.csv - summary statistics\n")
cat("  pathway_consistency_matrix.csv - cross-cohort consistency\n")
cat("  top_consensus_pathways.csv - high-confidence consensus pathways\n")
cat("  pathway_volcano_5panel.pdf - volcano plots\n")
cat("  pathway_heatmap_consensus.pdf - cross-cohort heatmap\n")
cat("  pathway_boxplot_top4.pdf - top pathway boxplots\n")
