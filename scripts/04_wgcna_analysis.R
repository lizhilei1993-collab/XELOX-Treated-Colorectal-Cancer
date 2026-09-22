# ============================================================
# 04_WGCNA_ANALYSIS.R
# WGCNA Co-expression network analysis
# XELOX Resistance Study
# ============================================================

Sys.setenv(TMPDIR = "C:/temp", TMP = "C:/temp", TEMP = "C:/temp")
.libPaths(c("C:/Rlibs", .libPaths()))

library(WGCNA)
library(GEOquery)
allowWGCNAThreads()

PROJECT_ROOT <- "C:/xelox_research"
DATA_GEO_DIR  <- file.path(PROJECT_ROOT, "data", "geo")
DATA_PROC_DIR <- file.path(PROJECT_ROOT, "data", "processed")
RESULTS_TAB_DIR <- file.path(PROJECT_ROOT, "results", "tables")
RESULTS_FIG_DIR <- file.path(PROJECT_ROOT, "results", "figures")

dir.create(RESULTS_FIG_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== XELOX Resistance Study: WGCNA Co-expression Network ===\n\n")

# ============================================================
# 1. Load and integrate expression data
# ============================================================
cat("[1] Loading GSE39582 as primary WGCNA dataset\n")

gse_eset <- readRDS(file.path(DATA_GEO_DIR, "GSE39582_eset.rds"))
expr_all <- exprs(gse_eset)
cat(sprintf("  Expression matrix: %d genes x %d samples\n", nrow(expr_all), ncol(expr_all)))

# Filter lowly expressed genes
keep <- rowSums(expr_all > 0, na.rm = TRUE) >= 10
expr_f <- expr_all[keep, ]
cat(sprintf("  After low-expression filtering: %d genes\n", nrow(expr_f)))

# Filter to top variable genes (WGCNA memory constraint: ~20GB RAM for 54K genes)
# Standard practice: use top 5000-8000 most variable genes
gene_vars <- apply(expr_f, 1, var, na.rm = TRUE)
n_top_genes <- min(8000, length(gene_vars))
top_genes <- order(gene_vars, decreasing = TRUE)[1:n_top_genes]
expr_f <- expr_f[top_genes, ]
cat(sprintf("  After variance filtering (top %d): %d genes\n", n_top_genes, nrow(expr_f)))

# Load group information
groups39582 <- read.csv(file.path(RESULTS_TAB_DIR, "GSE39582_xelox_groups.csv"),
                         stringsAsFactors = FALSE)
gs <- groups39582[match(colnames(expr_f), groups39582$sample_id), ]

# Create trait matrix
trait_data <- data.frame(
  resistant = ifelse(gs$group == "resistant", 1, 0),
  sensitive = ifelse(gs$group == "sensitive", 1, 0),
  rfs_event = gs$rfs_event,
  rfs_delay = gs$rfs_delay,
  age = gs$age,
  sex = ifelse(gs$sex == "M", 1, ifelse(gs$sex == "F", 0, NA)),
  stringsAsFactors = FALSE
)
rownames(trait_data) <- colnames(expr_f)
cat(sprintf("  Trait data: %d samples x %d traits\n", nrow(trait_data), ncol(trait_data)))

# Remove samples with missing traits
valid_samples <- complete.cases(trait_data[, c("resistant", "sensitive", "rfs_event", "rfs_delay")])
expr_f <- expr_f[, valid_samples]
trait_data <- trait_data[valid_samples, ]
cat(sprintf("  After trait filtering: %d samples\n", ncol(expr_f)))

# ============================================================
# 2. Pick soft threshold
# ============================================================
cat("\n[2] Picking soft threshold\n")

# Transpose: WGCNA expects samples as rows
datExpr <- as.data.frame(t(expr_f))

# Check for genes with too many missing values
gsg <- goodSamplesGenes(datExpr, verbose = 1)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  cat(sprintf("  After sample/gene QC: %d samples x %d genes\n", nrow(datExpr), ncol(datExpr)))
}

# Choose a set of soft-thresholding powers
powers <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 1)

# Plot scale independence
png(file.path(RESULTS_FIG_DIR, "WGCNA_scale_independence.png"), width = 800, height = 600)
par(mfrow = c(1, 2))
cex1 <- 0.8
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.85, col = "red", lty = 2)

# Mean connectivity
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5],
     labels = powers, cex = cex1, col = "red")
dev.off()
cat("  Scale independence plot saved\n")

# Determine optimal power
sft_threshold <- 0.85
r_squared <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
best_power <- sft$fitIndices$Power[which.max(r_squared[r_squared >= sft_threshold])]
if (length(best_power) == 0) {
  best_power <- sft$fitIndices$Power[which.max(r_squared)]
}
cat(sprintf("  Selected power: %d (R^2 = %.3f)\n", best_power, max(r_squared)))

# ============================================================
# 3. Build network and identify modules
# ============================================================
cat("\n[3] Building network and identifying modules\n")

net <- blockwiseModules(datExpr, power = best_power,
                        TOMType = "unsigned",
                        minModuleSize = 30,
                        reassignThreshold = 0,
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE,
                        pamRespectsDendro = FALSE,
                        verbose = 1)

module_labels <- net$colors
module_colors <- labels2colors(module_labels)

cat(sprintf("  Modules identified: %d\n", length(unique(module_labels)) - 1))
cat("  Module sizes:\n")
print(table(module_labels))

# ============================================================
# 4. Module-trait correlations
# ============================================================
cat("\n[4] Module-trait correlations\n")

# Calculate module eigengenes
MEs <- net$MEs
colnames(MEs) <- paste0("ME", colnames(MEs))

# Use only numeric trait columns
trait_cols <- c("resistant", "sensitive", "rfs_event", "rfs_delay")
trait_data_num <- as.data.frame(lapply(trait_data[, trait_cols], as.numeric))
rownames(trait_data_num) <- rownames(trait_data)

# Correlation
module_trait_cor <- cor(MEs, trait_data_num, use = "pairwise.complete.obs")
module_trait_pval <- corPvalueStudent(module_trait_cor, nrow(datExpr))

# Display correlation heatmap
png(file.path(RESULTS_FIG_DIR, "WGCNA_module_trait_heatmap.png"), width = 800, height = 600)
par(mar = c(6, 8, 3, 3))
labeledHeatmap(Matrix = module_trait_cor,
               xLabels = names(trait_data_num),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = sprintf("%.2f\n(%s)", 
                                    round(module_trait_cor, 2),
                                    ifelse(module_trait_pval < 0.05, "*", " ")),
               cex.text = 0.6,
               zlim = c(-1, 1),
               main = "Module-trait relationships")
dev.off()
cat("  Module-trait heatmap saved\n")

# Print significant correlations
cat("\n  Significant module-trait correlations (p < 0.05):\n")
sig_cors <- which(module_trait_pval < 0.05, arr.ind = TRUE)
if (nrow(sig_cors) > 0) {
  for (i in 1:nrow(sig_cors)) {
    cat(sprintf("    %s - %s: r=%.3f, p=%.4f\n",
                rownames(module_trait_cor)[sig_cors[i, 1]],
                colnames(module_trait_cor)[sig_cors[i, 2]],
                module_trait_cor[sig_cors[i, 1], sig_cors[i, 2]],
                module_trait_pval[sig_cors[i, 1], sig_cors[i, 2]]))
  }
}

# ============================================================
# 5. Gene module membership and significance
# ============================================================
cat("\n[5] Gene significance and module membership\n")

# Gene significance for resistance
gene_significance <- as.data.frame(cor(datExpr, trait_data_num$resistant, use = "pairwise.complete.obs"))
colnames(gene_significance) <- "GS_resistant"
gene_significance$pval <- corPvalueStudent(as.numeric(gene_significance$GS_resistant), nrow(datExpr))

# Module membership
module_membership <- as.data.frame(cor(datExpr, MEs, use = "pairwise.complete.obs"))
mm_pval <- as.data.frame(corPvalueStudent(as.matrix(module_membership), nrow(datExpr)))

# Combine results
wgcna_results <- data.frame(
  Gene = colnames(datExpr),
  Module = module_colors,
  ModuleLabel = module_labels,
  GS_resistant = gene_significance$GS_resistant,
  GS_pval = gene_significance$pval,
  stringsAsFactors = FALSE
)

# Add module membership for each module
for (me_name in names(MEs)) {
  module_name <- gsub("^ME", "", me_name)
  wgcna_results[[paste0("MM_", module_name)]] <- module_membership[[me_name]]
  wgcna_results[[paste0("MM_", module_name, "_pval")]] <- mm_pval[[me_name]]
}

# Identify resistance-correlated modules
r_modules <- which(module_trait_pval[, "resistant"] < 0.05)
cat(sprintf("  Modules correlated with resistance: %d\n", length(r_modules)))

if (length(r_modules) > 0) {
  r_module_names <- rownames(module_trait_cor)[r_modules]
  r_module_colors <- gsub("^ME", "", r_module_names)
  cat("  Resistance modules:", paste(r_module_colors, collapse = ", "), "\n")
  
  # Extract genes from resistance modules
  resistance_genes <- wgcna_results[wgcna_results$Module %in% r_module_colors, ]
  resistance_genes <- resistance_genes[order(resistance_genes$GS_resistant, decreasing = TRUE), ]
  
  write.csv(resistance_genes,
            file.path(RESULTS_TAB_DIR, "WGCNA_resistance_module_genes.csv"), row.names = FALSE)
  cat(sprintf("  %d genes in resistance modules saved\n", nrow(resistance_genes)))
}

# Save full results
write.csv(wgcna_results,
          file.path(RESULTS_TAB_DIR, "WGCNA_all_gene_module_membership.csv"), row.names = FALSE)

# ============================================================
# 6. Top hub genes in each module
# ============================================================
cat("\n[6] Extracting hub genes\n")

hub_genes <- list()

# Create color-to-label mapping
color_to_label <- setNames(module_labels, labels2colors(module_labels))

for (module in unique(wgcna_results$Module)) {
  module_genes <- wgcna_results[wgcna_results$Module == module, ]
  
  # Hub genes: high GS and high MM
  # Map color name back to numeric label for ME/MM column names
  module_label <- color_to_label[module]
  if (is.na(module_label)) next
  
  me_name <- paste0("ME", module_label)
  mm_col <- paste0("MM_", me_name)
  
  if (mm_col %in% colnames(module_genes)) {
    # Score = |GS| * |MM|
    module_genes$HubScore <- abs(module_genes$GS_resistant) * abs(module_genes[[mm_col]])
    module_genes <- module_genes[order(module_genes$HubScore, decreasing = TRUE), ]
    
    top_hubs <- head(module_genes, 20)
    hub_genes[[module]] <- top_hubs
    cat(sprintf("  Module %s (%s): %d hub genes extracted\n", module, me_name, nrow(top_hubs)))
  }
}

# Save top hub genes
hub_df <- do.call(rbind, hub_genes)
write.csv(hub_df, file.path(RESULTS_TAB_DIR, "WGCNA_hub_genes.csv"), row.names = FALSE)
cat(sprintf("  Hub genes saved: %d genes\n", nrow(hub_df)))

# ============================================================
# Summary
# ============================================================
cat("\n=== WGCNA Analysis Complete ===\n")
cat(sprintf("  Total modules (incl. grey): %d\n", length(unique(module_labels))))
cat(sprintf("  Modules correlated with resistance: %d\n", length(r_modules)))
cat(sprintf("  Hub genes saved: %d\n", nrow(hub_df)))

save.image(file.path(DATA_PROC_DIR, "04_wgcna_analysis.RData"))
cat("Workspace saved\n")
