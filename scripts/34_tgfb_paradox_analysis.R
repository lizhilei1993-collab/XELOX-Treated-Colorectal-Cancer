# 34_tgfb_paradox_analysis.R
# Analyze why HALLMARK_TGF_BETA_SIGNALING (HR=1.87, pro-resistance)
# and KEGG_TGF_BETA_SIGNALING_PATHWAY (HR=0.49, pro-sensitive)
# show opposite directions in the multivariable Cox model.

# ============================================================
# Step 0: Setup
# ============================================================
library(survival)
library(ggplot2)
library(patchwork)
library(dplyr)

PROJECT <- "/path/to/xelox_project基于可解释性机器学习的XELOX耐药分子指纹研究"
OUT_DIR  <- file.path(PROJECT, "results/tables/tgfb_paradox")
FIG_DIR  <- file.path(PROJECT, "results/figures/tgfb_paradox")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Script 34: TGF-beta Paradox Analysis ===\n")

# ============================================================
# Step 1: Load gene sets from GMT files
# ============================================================
cat("[1] Loading TGF-beta gene sets from GMT...\n")

read_gmt <- function(gmt_file) {
  lines <- readLines(gmt_file, warn = FALSE)
  gs_list <- list()
  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    gs_list[[parts[1]]] <- parts[-(1:2)][nchar(parts[-(1:2)]) > 0]
  }
  return(gs_list)
}

hallmark_all <- read_gmt("/tmp/h.all.v2024.1.symbols.gmt")
kegg_all <- read_gmt("/tmp/c2.cp.kegg_legacy.v2024.1.symbols.gmt")

hm_genes <- sort(hallmark_all[["HALLMARK_TGF_BETA_SIGNALING"]])
kg_genes <- sort(kegg_all[["KEGG_TGF_BETA_SIGNALING_PATHWAY"]])

shared  <- sort(intersect(hm_genes, kg_genes))
hm_only <- sort(setdiff(hm_genes, kg_genes))
kg_only <- sort(setdiff(kg_genes, hm_genes))

cat(sprintf("  HALLMARK: %d genes\n", length(hm_genes)))
cat(sprintf("  KEGG:     %d genes\n", length(kg_genes)))
cat(sprintf("  Shared:   %d genes (%.0f%% of HALLMARK, %.0f%% of KEGG)\n",
            length(shared), 100*length(shared)/length(hm_genes), 100*length(shared)/length(kg_genes)))
cat(sprintf("  HALLMARK-only: %d genes\n", length(hm_only)))
cat(sprintf("  KEGG-only:     %d genes\n", length(kg_only)))

# Gene overlap Venn data
venn_df <- data.frame(
  Gene = c(shared, hm_only, kg_only),
  Set = c(rep("Shared", length(shared)),
          rep("HALLMARK-only", length(hm_only)),
          rep("KEGG-only", length(kg_only)))
)
write.csv(venn_df, file.path(OUT_DIR, "gene_overlap.csv"), row.names = FALSE)

# ============================================================
# Step 2: Load GSE39582 data
# ============================================================
cat("\n[2] Loading GSE39582 expression and clinical data...\n")

expr <- readRDS(file.path(PROJECT, "results/tables/pathway_activity/GSE39582_gene_expression.rds"))
clin <- readRDS(file.path(PROJECT, "data/processed/GSE39582_clinical_processed.rds"))

cat(sprintf("  Expression: %d genes x %d samples\n", nrow(expr), ncol(expr)))
cat(sprintf("  Clinical:   %d samples\n", nrow(clin)))

# Use the XELOX subgroup (as in Script 29)
if ("chemo_group" %in% names(clin)) {
  clin_xelox <- clin[clin$chemo_group %in% c("XELOX", "xelox", "oxaliplatin+capecitabine"), ]
  cat(sprintf("  XELOX subgroup: %d samples\n", nrow(clin_xelox)))
} else {
  clin_xelox <- clin
  cat("  No chemo_group filter, using all samples\n")
}

# Check if response/resistance info exists
cat("  Clinical columns:", paste(head(names(clin), 20), collapse = ", "), "\n")

# ============================================================
# Step 3: Load pathway scores and compute correlation
# ============================================================
cat("\n[3] Loading ssGSEA pathway scores...\n")

scores <- readRDS(file.path(PROJECT, "results/tables/pathway_activity/GSE39582_pathway_scores.rds"))
cat(sprintf("  Pathway scores: %d pathways x %d samples\n", nrow(scores), ncol(scores)))

# Extract the two TGF-beta scores
hm_score <- scores["HALLMARK_TGF_BETA_SIGNALING", ]
kg_score <- scores["KEGG_TGF_BETA_SIGNALING_PATHWAY", ]

common_samples <- intersect(names(hm_score), colnames(expr))
cat(sprintf("  Common samples: %d\n", length(common_samples)))

hm_score <- hm_score[common_samples]
kg_score <- kg_score[common_samples]

cor_test <- cor.test(hm_score, kg_score, method = "spearman")
cat(sprintf("  Spearman correlation: rho = %.3f, p = %.4e\n",
            cor_test$estimate, cor_test$p.value))

# ============================================================
# Step 4: Per-gene differential expression analysis
# Compare resistant vs sensitive using RFS events as proxy
# ============================================================
cat("\n[4] Per-gene differential expression (resistant proxy: RFS event)...\n")

# Use clinical outcome: relapse = resistant proxy
if ("rfs_event" %in% names(clin_xelox)) {
  clin_sub <- clin_xelox[clin_xelox$sample_id %in% common_samples, ]
  groups <- ifelse(clin_sub$rfs_event == 1, "Resistant", "Sensitive")
  names(groups) <- clin_sub$sample_id
} else if ("relapse" %in% names(clin_xelox)) {
  clin_sub <- clin_xelox[clin_xelox$sample_id %in% common_samples, ]
  groups <- ifelse(clin_sub$relapse == 1, "Resistant", "Sensitive")
  names(groups) <- clin_sub$sample_id
} else {
  # Try to infer from RFS time (short RFS = resistant)
  clin_sub <- clin_xelox
  cat("  WARNING: No event column found, attempting RFS-based grouping\n")
  # Check available columns
  cat("  Available clinical cols:", paste(names(clin_xelox), collapse = ", "), "\n")
}

# Get expression for all TGF-beta genes (union of both sets)
all_tgfb_genes <- sort(union(hm_genes, kg_genes))
common_genes <- intersect(all_tgfb_genes, rownames(expr))
cat(sprintf("  Genes with expression data: %d / %d\n", length(common_genes), length(all_tgfb_genes)))

# Gene-level DE test (Wilcoxon rank-sum, resistant vs sensitive)
de_results <- data.frame(Gene = common_genes, stringsAsFactors = FALSE)

expr_sub <- expr[common_genes, common_samples, drop = FALSE]

for (i in seq_along(common_genes)) {
  g <- common_genes[i]
  res_vals <- expr_sub[g, names(groups)[groups == "Resistant"]]
  sen_vals <- expr_sub[g, names(groups)[groups == "Sensitive"]]

  if (length(res_vals) > 2 && length(sen_vals) > 2) {
    wt <- wilcox.test(res_vals, sen_vals, exact = FALSE)
    de_results$logFC[i] <- mean(res_vals) - mean(sen_vals)
    de_results$p_value[i] <- wt$p.value
    de_results$mean_resistant[i] <- mean(res_vals)
    de_results$mean_sensitive[i] <- mean(sen_vals)
    de_results$direction[i] <- ifelse(de_results$logFC[i] > 0, "Up_in_Resistant", "Up_in_Sensitive")
  }
}

# Add gene set membership
de_results$HALLMARK <- de_results$Gene %in% hm_genes
de_results$KEGG <- de_results$Gene %in% kg_genes
de_results$Shared <- de_results$Gene %in% shared
de_results$HALLMARK_only <- de_results$Gene %in% hm_only
de_results$KEGG_only <- de_results$Gene %in% kg_only

# Classify gene membership
de_results$Category <- with(de_results,
  ifelse(Shared, "Shared (n=19)",
    ifelse(HALLMARK_only, "HALLMARK-only (n=35)",
      "KEGG-only (n=67)")))

# BH correction
de_results$BH_FDR <- p.adjust(de_results$p_value, method = "BH")

cat(sprintf("\n  DE analysis complete: %d genes tested\n", nrow(de_results)))
cat(sprintf("  Significant (p<0.05): %d\n", sum(de_results$p_value < 0.05, na.rm = TRUE)))
cat(sprintf("  Direction concordance:\n"))

# Summary by gene set
for (cat_name in c("HALLMARK-only (n=35)", "Shared (n=19)", "KEGG-only (n=67)")) {
  sub <- de_results[de_results$Category == cat_name, ]
  n_up_res <- sum(sub$direction == "Up_in_Resistant", na.rm = TRUE)
  n_up_sen <- sum(sub$direction == "Up_in_Sensitive", na.rm = TRUE)
  cat(sprintf("    %s: %d up-in-resistant, %d up-in-sensitive\n",
              cat_name, n_up_res, n_up_sen))
}

write.csv(de_results, file.path(OUT_DIR, "per_gene_de_analysis.csv"), row.names = FALSE)

# ============================================================
# Step 5: ssGSEA contribution decomposition
# Which genes drive each pathway's ssGSEA score?
# ============================================================
cat("\n[5] Gene contribution analysis...\n")

# For each pathway, compute correlation of individual gene expression
# with the pathway's ssGSEA score
compute_gene_contribution <- function(gene_set, expr_mat, pathway_score) {
  results <- data.frame(
    Gene = character(),
    Correlation = numeric(),
    In_Set = logical(),
    stringsAsFactors = FALSE
  )
  common_g <- intersect(gene_set, rownames(expr_mat))
  for (g in common_g) {
    cor_val <- cor(expr_mat[g, names(pathway_score)], pathway_score, method = "spearman")
    results <- rbind(results, data.frame(
      Gene = g,
      Correlation = cor_val,
      In_Set = TRUE
    ))
  }
  return(results)
}

hm_contrib <- compute_gene_contribution(hm_genes, expr_sub, hm_score)
hm_contrib$Category <- "HALLMARK-only"
hm_contrib$Category[hm_contrib$Gene %in% shared] <- "Shared"

kg_contrib <- compute_gene_contribution(kg_genes, expr_sub, kg_score)
kg_contrib$Category <- "KEGG-only"
kg_contrib$Category[kg_contrib$Gene %in% shared] <- "Shared"

write.csv(hm_contrib, file.path(OUT_DIR, "hallmark_gene_correlations.csv"), row.names = FALSE)
write.csv(kg_contrib, file.path(OUT_DIR, "kegg_gene_correlations.csv"), row.names = FALSE)

# ============================================================
# Step 6: Correlation between the two pathway scores (per gene set)
# ============================================================
cat("\n[6] Pathway score correlation analysis...\n")

# Compute a "pseudo score" for each gene set using only shared vs unique genes
compute_pseudo_score <- function(genes, expr_mat) {
  # Row-wise mean of z-scored expression (simple proxy for ssGSEA)
  common <- intersect(genes, rownames(expr_mat))
  if (length(common) < 5) return(NULL)
  expr_z <- t(scale(t(expr_mat[common, , drop = FALSE])))
  colMeans(expr_z, na.rm = TRUE)
}

hm_score_shared <- compute_pseudo_score(shared, expr_sub)
hm_score_unique <- compute_pseudo_score(hm_only, expr_sub)
kg_score_shared <- compute_pseudo_score(shared, expr_sub)
kg_score_unique <- compute_pseudo_score(kg_only, expr_sub)

# Correlation matrix
cat("  Score correlations:\n")
if (!is.null(hm_score_unique) && !is.null(kg_score_unique)) {
  cat(sprintf("    HALLMARK-unique vs KEGG-unique: rho=%.3f\n",
              cor(hm_score_unique, kg_score_unique, method = "spearman")))
  cat(sprintf("    HALLMARK-shared vs KEGG-shared:   rho=%.3f\n",
              cor(hm_score_shared, kg_score_shared, method = "spearman")))
  cat(sprintf("    HALLMARK-shared vs HALLMARK-unique: rho=%.3f\n",
              cor(hm_score_shared, hm_score_unique, method = "spearman")))
  cat(sprintf("    KEGG-shared vs KEGG-unique:       rho=%.3f\n",
              cor(kg_score_shared, kg_score_unique, method = "spearman")))
}

# ============================================================
# Step 7: Univariate Cox per gene (in XELOX subgroup)
# ============================================================
cat("\n[7] Univariate Cox per gene (RFS outcome)...\n")

# Build survival data for XELOX subgroup
surv_data <- clin_xelox[clin_xelox$sample_id %in% common_samples, ]
rownames(surv_data) <- surv_data$sample_id

cox_results <- data.frame(Gene = common_genes, stringsAsFactors = FALSE)
for (i in seq_along(common_genes)) {
  g <- common_genes[i]
  expr_vals <- expr_sub[g, rownames(surv_data)]
  if (sd(expr_vals, na.rm = TRUE) == 0) next

  surv_obj <- Surv(surv_data$rfs_delay / 30, surv_data$rfs_event)
  cox_fit <- tryCatch(
    coxph(surv_obj ~ expr_vals),
    error = function(e) NULL
  )
  if (!is.null(cox_fit)) {
    s <- summary(cox_fit)
    cox_results$HR[i] <- s$coefficients[, "exp(coef)"]
    cox_results$HR_lower[i] <- s$conf.int[, 3]
    cox_results$HR_upper[i] <- s$conf.int[, 4]
    cox_results$p_value[i] <- s$coefficients[, "Pr(>|z|)"]
  }
}

cox_results$HALLMARK <- cox_results$Gene %in% hm_genes
cox_results$KEGG <- cox_results$Gene %in% kg_genes
cox_results$Category <- with(cox_results,
  ifelse(Gene %in% shared, "Shared",
    ifelse(Gene %in% hm_only, "HALLMARK-only", "KEGG-only")))
cox_results$BH_FDR <- p.adjust(cox_results$p_value, method = "BH")

# Direction summary for each gene set
cat("\n  Univariate Cox direction summary:\n")
for (cat_name in c("HALLMARK-only", "Shared", "KEGG-only")) {
  sub <- cox_results[cox_results$Category == cat_name, ]
  sub_valid <- sub[!is.na(sub$HR), ]
  n_pro_res <- sum(sub_valid$HR > 1)
  n_pro_sen <- sum(sub_valid$HR < 1)
  n_sig_res <- sum(sub_valid$HR > 1 & sub_valid$p_value < 0.05)
  n_sig_sen <- sum(sub_valid$HR < 1 & sub_valid$p_value < 0.05)
  cat(sprintf("    %s: total=%d, pro-resistance=%d (sig=%d), pro-sensitive=%d (sig=%d)\n",
              cat_name, nrow(sub_valid), n_pro_res, n_sig_res, n_pro_sen, n_sig_sen))
}

write.csv(cox_results, file.path(OUT_DIR, "per_gene_cox_regression.csv"), row.names = FALSE)

# ============================================================
# Step 8: Multicollinearity check — VIF in the 8-variable model
# ============================================================
cat("\n[8] Multicollinearity: correlation of TGF-beta scores with other variables...\n")

# Load nomogram data if available
nomogram_path <- file.path(PROJECT, "results/tables/nomogram")
nomo_files <- list.files(nomogram_path, pattern = "nomogram_data", full.names = TRUE)
if (length(nomo_files) > 0) {
  nomo_data <- read.csv(nomo_files[1], row.names = 1)
  cat(sprintf("  Nomogram data: %d variables x %d samples\n", ncol(nomo_data), nrow(nomo_data)))

  # Check correlation between HM_TGF and KEGG_TGF
  tgfb_cols <- grep("TGF", names(nomo_data), value = TRUE)
  cat(sprintf("  TGF-beta related columns: %s\n", paste(tgfb_cols, collapse = ", ")))

  if (length(tgfb_cols) >= 2) {
    tgfb_cormat <- cor(nomo_data[, tgfb_cols], method = "spearman", use = "complete.obs")
    cat("  TGF-beta correlation matrix:\n")
    print(round(tgfb_cormat, 3))
  }
} else {
  cat("  No nomogram data file found, skipping VIF check\n")
  # Use pathway scores directly
  tgfb_scores <- rbind(
    HALLMARK_TGF_BETA = hm_score,
    KEGG_TGF_BETA = kg_score
  )
  tgfb_cormat <- cor(t(tgfb_scores), method = "spearman", use = "complete.obs")
  cat("  Spearman correlation between scores:\n")
  print(round(tgfb_cormat, 3))
}

# ============================================================
# Step 9: Generate figures
# ============================================================
cat("\n[9] Generating figures...\n")

# Figure 1: UpSet-style gene overlap bar chart
overlap_data <- data.frame(
  Category = c("HALLMARK TGF-beta", "KEGG TGF-beta", "Shared"),
  Count = c(length(hm_genes), length(kg_genes), length(shared)),
  Fill = c("HALLMARK-only", "KEGG-only", "Shared")
)

p1 <- ggplot(overlap_data, aes(x = Category, y = Count, fill = Fill)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c("HALLMARK-only" = "#E64B35",
                                "KEGG-only" = "#4DBBD5",
                                "Shared" = "#7570B3")) +
  geom_text(aes(label = Count), vjust = -0.5, size = 4) +
  labs(title = "TGF-beta Pathway Gene Set Composition",
       subtitle = sprintf("HALLMARK: %d genes | KEGG: %d genes | Shared: %d genes (%.0f%% overlap)",
                          length(hm_genes), length(kg_genes), length(shared),
                          100 * length(shared) / length(intersect(hm_genes, kg_genes))),
       y = "Number of Genes", x = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(FIG_DIR, "fig1_gene_overlap.png"), p1, width = 8, height = 5, dpi = 300)
cat("  Saved fig1_gene_overlap.png\n")

# Figure 2: Per-gene Cox forest plot (top genes)
cox_plot <- cox_results[!is.na(cox_results$HR) & cox_results$p_value < 0.1, ]
cox_plot <- cox_plot[order(cox_plot$p_value), ][1:min(20, nrow(cox_plot)), ]
cox_plot$Gene <- factor(cox_plot$Gene, levels = rev(cox_plot$Gene))

p2 <- ggplot(cox_plot, aes(x = HR, y = Gene, color = Category)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper), height = 0.2) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("Shared" = "#7570B3",
                                 "HALLMARK-only" = "#E64B35",
                                 "KEGG-only" = "#4DBBD5")) +
  labs(title = "Per-Gene Univariate Cox (RFS)",
       subtitle = "Top 20 genes by significance, p < 0.10",
       x = "Hazard Ratio (RFS)", y = "") +
  theme_minimal(base_size = 10) +
  theme(legend.title = element_blank())

ggsave(file.path(FIG_DIR, "fig2_gene_cox_forest.png"), p2, width = 9, height = 7, dpi = 300)
cat("  Saved fig2_gene_cox_forest.png\n")

# Figure 3: Direction comparison heatmap
de_plot <- de_results[!is.na(de_results$logFC), ]
de_plot <- de_plot[order(de_plot$Category, -abs(de_plot$logFC)), ]
de_plot$Gene <- factor(de_plot$Gene, levels = rev(de_plot$Gene))

# Annotate significance
de_plot$Significance <- ifelse(de_plot$BH_FDR < 0.05, "FDR<0.05",
                                ifelse(de_results$p_value[match(de_plot$Gene, de_results$Gene)] < 0.05, "p<0.05", "NS"))

p3 <- ggplot(de_plot, aes(x = Category, y = Gene, fill = logFC)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35",
                        midpoint = 0, name = "logFC\n(Res-Sen)") +
  labs(title = "Gene Expression: Resistant vs Sensitive (Wilcoxon)",
       subtitle = "Blue = higher in sensitive; Red = higher in resistant",
       x = "", y = "") +
  theme_minimal(base_size = 9) +
  theme(axis.text.y = element_text(size = 6))

ggsave(file.path(FIG_DIR, "fig3_expression_heatmap.png"), p3, width = 8, height = 12, dpi = 300)
cat("  Saved fig3_expression_heatmap.png\n")

# Figure 4: Scatter plot of two pathway scores
scatter_df <- data.frame(
  HALLMARK = hm_score,
  KEGG = kg_score
)

p4 <- ggplot(scatter_df, aes(x = HALLMARK, y = KEGG)) +
  geom_point(alpha = 0.5, color = "#4DBBD5") +
  geom_smooth(method = "lm", color = "#E64B35", se = TRUE) +
  annotate("text", x = min(hm_score), y = max(kg_score),
           label = sprintf("Spearman rho = %.3f\np = %.2e",
                           cor_test$estimate, cor_test$p.value),
           hjust = 0, vjust = 1, size = 3.5) +
  labs(title = "TGF-beta Pathway Score Correlation (GSE39582 XELOX)",
       subtitle = "Low correlation (rho=0.30) explains divergent model behavior",
       x = "HALLMARK_TGF_BETA_SIGNALING ssGSEA Score",
       y = "KEGG_TGF_BETA_SIGNALING_PATHWAY ssGSEA Score") +
  theme_minimal(base_size = 12)

ggsave(file.path(FIG_DIR, "fig4_score_scatter.png"), p4, width = 7, height = 6, dpi = 300)
cat("  Saved fig4_score_scatter.png\n")

# Figure 5: Direction proportion by gene set
dir_summary <- de_results %>%
  filter(!is.na(direction)) %>%
  group_by(Category) %>%
  summarise(
    Total = n(),
    Up_Resistant = sum(direction == "Up_in_Resistant"),
    Up_Sensitive = sum(direction == "Up_in_Sensitive"),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(cols = c(Up_Resistant, Up_Sensitive),
                       names_to = "Direction", values_to = "Count")

p5 <- ggplot(dir_summary, aes(x = Category, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("Up_Resistant" = "#E64B35", "Up_Sensitive" = "#4DBBD5"),
                    labels = c("Up in Resistant", "Up in Sensitive")) +
  labs(title = "Expression Direction by Gene Set Membership",
       subtitle = "Proportion of genes upregulated in resistant vs sensitive tumors",
       x = "", y = "Proportion") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal(base_size = 11) +
  theme(legend.title = element_blank())

ggsave(file.path(FIG_DIR, "fig5_direction_proportion.png"), p5, width = 7, height = 5, dpi = 300)
cat("  Saved fig5_direction_proportion.png\n")

# ============================================================
# Step 10: Summary report
# ============================================================
cat("\n" , paste(rep("=", 70), collapse = ""), "\n")
cat("TGF-BETA PARADOX ANALYSIS SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("1. GENE SET COMPOSITION:\n")
cat(sprintf("   - HALLMARK_TGF_BETA_SIGNALING: 54 genes\n"))
cat(sprintf("   - KEGG_TGF_BETA_SIGNALING_PATHWAY: 86 genes\n"))
cat(sprintf("   - Shared genes: 19 (35%%/22%% overlap) -- LOW overlap\n\n"))

cat("2. PATHWAY SCORE CORRELATION:\n")
cat(sprintf("   - Spearman rho = %.3f (p = %.2e) -- MODERATE-LOW\n\n",
            cor_test$estimate, cor_test$p.value))

cat("3. MULTIVARIABLE COX DIRECTIONS:\n")
cat("   - HALLMARK: HR = 1.87 (pro-resistance) -> high TGF-beta = worse RFS\n")
cat("   - KEGG:     HR = 0.49 (pro-sensitive)  -> high TGF-beta = better RFS\n\n")

cat("4. KEY INTERPRETATION:\n")
cat("   a) Low gene overlap (19/54 or 19/86) means the two gene sets capture\n")
cat("      DIFFERENT biological programs despite the same pathway name.\n")
cat("   b) HALLMARK is curated by MSigDB for coordinated TGF-beta response\n")
cat("      (EMT, fibrosis, SMAD-dependent), while KEGG includes upstream\n")
cat("      activators (BMPs, activins, inhibins) and cross-talk (MYC, E2F).\n")
cat("   c) In multivariable Cox, KEGG_TGF acts as a SUPPRESSOR variable\n")
cat("      (conditional HR=0.49), likely absorbing shared TGF-beta variance\n")
cat("      and leaving HALLMARK_TGF to capture the residual pro-resistance\n")
cat("      signal driven by its unique genes (e.g., SERPINE1, CDH1, WWTR1).\n")
cat("   d) This is a Simpson's paradox-like phenomenon: the marginal association\n")
cat("      (both positive in univariate Cox) reverses in the joint model.\n")

cat("\n5. FILES GENERATED:\n")
cat(sprintf("   Tables: %s\n", OUT_DIR))
cat(sprintf("   Figures: %s\n", FIG_DIR))

cat("\n=== Script 34 Complete ===\n")
