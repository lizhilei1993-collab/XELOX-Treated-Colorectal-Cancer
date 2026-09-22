# ============================================================
# 07_ENRICHMENT_ANALYSIS.R (v2)
# GO/KEGG/DO/GSEA enrichment analysis for XELOX resistance
#   - Over-Representation Analysis (ORA): GO (BP/CC/MF), KEGG, DO
#   - Gene Set Enrichment Analysis (GSEA): preranked by DEG logFC
#   - All gene pool tiers + per-dataset GSEA
# XELOX Resistance Study
# ============================================================
# Usage:
#   Rscript scripts/07_enrichment_analysis.R
# ============================================================

Sys.setenv(HOME = "C:/temp", TMPDIR = "C:/temp", TMP = "C:/temp", TEMP = "C:/temp", R_LIBS = "C:/Rlibs")
Sys.setlocale("LC_ALL", "C")
.libPaths(c("C:/Rlibs", .libPaths()))

PROJECT_ROOT     <- "/path/to/xelox_project"
SCRIPT_DIR       <- file.path(PROJECT_ROOT, "scripts")
DATA_GEO_DIR     <- file.path(PROJECT_ROOT, "data", "geo")
DATA_PROC_DIR    <- file.path(PROJECT_ROOT, "data", "processed")
RESULTS_TAB_DIR  <- file.path(PROJECT_ROOT, "results", "tables")
RESULTS_FIG_DIR  <- file.path(PROJECT_ROOT, "results", "figures")
ENRICH_DIR       <- file.path(RESULTS_TAB_DIR, "enrichment")
GSEA_DIR         <- file.path(RESULTS_FIG_DIR, "gsea")

dir.create(ENRICH_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(GSEA_DIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("XELOX Resistance: Comprehensive Enrichment Analysis v2\n")
cat("  ORA (GO / KEGG / DO) + GSEA (Preranked)\n")
cat("============================================================\n\n")

# ============================================================
# 0. Load packages
# ============================================================
cat("=== [0] Loading packages ===\n\n")

required_pkgs <- c(
  "clusterProfiler", "enrichplot", "org.Hs.eg.db", "GO.db",
  "DOSE", "ggplot2", "ggridges", "pheatmap"
)
for (pkg in required_pkgs) {
  tryCatch({
    library(pkg, character.only = TRUE, quietly = TRUE)
    cat(sprintf("  \u2713 %s v%s\n", pkg, packageVersion(pkg)))
  }, error = function(e) {
    cat(sprintf("  \u2717 %s - %s\n", pkg, conditionMessage(e)))
  })
}
cat("\n")

# ============================================================
# 1. Load gene pool
# ============================================================
cat("=== [1] Loading gene pool ===\n\n")

gp_file <- file.path(RESULTS_TAB_DIR, "XELOX_final_gene_pool.csv")
if (!file.exists(gp_file)) {
  stop("Gene pool not found. Run 05_gene_pool_integration.R first: ", gp_file)
}

gp <- read.csv(gp_file, stringsAsFactors = FALSE)
cat(sprintf("  Total genes in pool: %d\n", nrow(gp)))
cat(sprintf("  Priority breakdown:\n"))
print(table(gp$Priority))

# Define gene sets by tier
tiers <- list(
  All                 = gp$Gene,
  Tier1_Cape_DEG_WGCNA = gp$Gene[gp$Priority == "Tier1_Capecitabine_DEG_WGCNA"],
  Tier2_DEG_WGCNA     = gp$Gene[gp$Priority == "Tier2_DEG_WGCNA"],
  Tier3_DEG_Metab     = gp$Gene[gp$Priority == "Tier3_DEG_Metabolism"],
  Tier4_MultiDEG      = gp$Gene[gp$Priority == "Tier4_MultiDataset_DEG"],
  Tier5_WGCNA_Hub     = gp$Gene[gp$Priority == "Tier5_WGCNA_Hub"]
)
# Composite sets
tiers[["HighPriority"]] <- gp$Gene[
  gp$Priority %in% c("Tier2_DEG_WGCNA", "Tier3_DEG_Metabolism", "Tier5_WGCNA_Hub")
]
tiers[["Tier1_2"]] <- gp$Gene[
  gp$Priority %in% c("Tier1_Capecitabine_DEG_WGCNA", "Tier2_DEG_WGCNA")
]

for (nm in names(tiers)) {
  cat(sprintf("  %-24s: %d genes\n", nm, length(tiers[[nm]])))
}
cat("\n")

# ============================================================
# 2. Gene symbol -> Entrez ID mapping
# ============================================================
cat("=== [2] Mapping gene symbols to Entrez IDs ===\n\n")

all_genes <- unique(gp$Gene)
entrez_map <- tryCatch({
  bitr(all_genes, fromType = "SYMBOL", toType = "ENTREZID",
       OrgDb = org.Hs.eg.db)
}, error = function(e) {
  cat(sprintf("  WARNING: bitr failed: %s\n", conditionMessage(e)))
  cat("  Trying with clean symbols only...\n")
  clean_genes <- grep("^[A-Za-z0-9.-]+$", all_genes, value = TRUE)
  bitr(clean_genes, fromType = "SYMBOL", toType = "ENTREZID",
       OrgDb = org.Hs.eg.db)
})

cat(sprintf("  Input symbols: %d\n", length(all_genes)))
cat(sprintf("  Mapped to Entrez: %d\n", nrow(entrez_map)))
cat(sprintf("  Unmapped: %d\n", length(all_genes) - nrow(entrez_map)))

if (length(all_genes) - nrow(entrez_map) > 0) {
  unmapped <- setdiff(all_genes, unique(entrez_map$SYMBOL))
  cat(sprintf("  Unmapped examples: %s\n",
              paste(head(unmapped, 10), collapse = ", ")))
}
cat("\n")

# Named list: Entrez IDs by tier
entrez_tiers <- list()
for (nm in names(tiers)) {
  entrez_tiers[[nm]] <- entrez_map$ENTREZID[entrez_map$SYMBOL %in% tiers[[nm]]]
  cat(sprintf("  %-24s: %d Entrez IDs\n", nm, length(entrez_tiers[[nm]])))
}
cat("\n")

# ============================================================
# 3. GO Enrichment (ORA) - ALL tiers
# ============================================================
cat("=== [3] GO ORA: All tiers × All ontologies ===\n\n")

go_results <- list()
go_ontologies <- c("BP", "CC", "MF")
ora_sets <- names(tiers)  # test all tiers

for (ont in go_ontologies) {
  cat(sprintf("--- GO-%s ---\n", ont))
  
  for (set_name in ora_sets) {
    genes <- entrez_tiers[[set_name]]
    if (length(genes) < 5) {
      cat(sprintf("  %-24s: SKIP (< 5 genes)\n", set_name))
      next
    }
    
    ego <- tryCatch({
      enrichGO(gene          = genes,
               OrgDb         = org.Hs.eg.db,
               keyType       = "ENTREZID",
               ont           = ont,
               pvalueCutoff  = 0.05,
               pAdjustMethod = "BH",
               minGSSize     = 5,
               maxGSSize     = 800,
               readable      = TRUE)
    }, error = function(e) {
      cat(sprintf("    ERROR: %s\n", conditionMessage(e)))
      return(NULL)
    })
    
    key <- paste(set_name, ont, sep = "_")
    if (!is.null(ego) && nrow(ego) > 0) {
      go_results[[key]] <- ego
      cat(sprintf("  %-24s: %d terms\n", set_name, nrow(ego)))
      
      outfile <- file.path(ENRICH_DIR, sprintf("GO_%s_%s.csv", set_name, ont))
      write.csv(as.data.frame(ego), outfile, row.names = FALSE)
    } else {
      go_results[[key]] <- NULL
      cat(sprintf("  %-24s: 0 terms\n", set_name))
    }
  }
  cat("\n")
}

# ============================================================
# 4. KEGG Pathway Enrichment (ORA) - ALL tiers
# ============================================================
cat("=== [4] KEGG ORA: All tiers ===\n\n")

kegg_results <- list()

for (set_name in ora_sets) {
  genes <- entrez_tiers[[set_name]]
  if (length(genes) < 5) {
    cat(sprintf("  %-24s: SKIP (< 5 genes)\n", set_name))
    next
  }
  
  kk <- tryCatch({
    enrichKEGG(gene          = genes,
               organism       = "hsa",
               keyType        = "ncbi-geneid",
               pvalueCutoff   = 0.05,
               pAdjustMethod  = "BH",
               minGSSize      = 5,
               maxGSSize      = 800)
  }, error = function(e) {
    cat(sprintf("    ERROR: %s\n", conditionMessage(e)))
    return(NULL)
  })
  
  if (!is.null(kk) && nrow(kk) > 0) {
    kegg_results[[set_name]] <- kk
    cat(sprintf("  %-24s: %d pathways\n", set_name, nrow(kk)))
    
    outfile <- file.path(ENRICH_DIR, sprintf("KEGG_%s.csv", set_name))
    write.csv(as.data.frame(kk), outfile, row.names = FALSE)
  } else {
    kegg_results[[set_name]] <- NULL
    cat(sprintf("  %-24s: 0 pathways\n", set_name))
  }
}
cat("\n")

# ============================================================
# 5. DisGeNET (DGN) Disease Enrichment
# ============================================================
cat("=== [5] DisGeNET (DGN) Disease Enrichment ===\n\n")
cat("  NOTE: Using DGN (DisGeNET) instead of DO/HDO because\n")
cat("  DOSE v4.5.0 delegates DO enrichment to HDO.db, which\n")
cat("  could not be installed due to Windows/C-locale path issues.\n")
cat("  DGN data is built into DOSE and has more comprehensive\n")
cat("  disease-gene associations, making it suitable for\n")
cat("  XELOX resistance phenotype enrichment.\n\n")

dgn_results <- list()

for (set_name in ora_sets) {
  genes <- entrez_tiers[[set_name]]
  if (length(genes) < 5) {
    cat(sprintf("  %-24s: SKIP (< 5 genes)\n", set_name))
    next
  }
  
  dgn_ego <- tryCatch({
    enrichDGN(gene          = genes,
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              minGSSize     = 5,
              maxGSSize     = 800,
              readable      = TRUE)
  }, error = function(e) {
    cat(sprintf("    ERROR: %s\n", conditionMessage(e)))
    return(NULL)
  })
  
  if (!is.null(dgn_ego) && nrow(dgn_ego) > 0) {
    dgn_results[[set_name]] <- dgn_ego
    cat(sprintf("  %-24s: %d DGN disease terms\n", set_name, nrow(dgn_ego)))
    
    outfile <- file.path(ENRICH_DIR, sprintf("DGN_%s.csv", set_name))
    write.csv(as.data.frame(dgn_ego), outfile, row.names = FALSE)
  } else {
    dgn_results[[set_name]] <- NULL
    cat(sprintf("  %-24s: 0 DGN terms\n", set_name))
  }
}
cat("\n")

# ============================================================
# 6. GSEA Preranked Analysis
# ============================================================
cat("=== [6] GSEA Preranked Analysis (DEG logFC-based) ===\n\n")

# -- 6a. Build ranking metric from cross-dataset DEG results --
cat("--- [6a] Building gene ranking from DEG mapped files ---\n\n")

deg_mapped_files <- list.files(RESULTS_TAB_DIR,
                                pattern = "^DEG_.*_mapped\\.csv$",
                                full.names = TRUE)
cat(sprintf("  Found %d DEG mapped files\n", length(deg_mapped_files)))

# Collect gene-level statistics across datasets
gene_stats_list <- list()

for (f in deg_mapped_files) {
  base_name <- basename(f)
  gse_id <- gsub("^DEG_|_mapped\\.csv$", "", base_name)
  
  df <- read.csv(f, stringsAsFactors = FALSE)
  
  # Remove rows without gene_symbol
  df <- df[!is.na(df$gene_symbol) & df$gene_symbol != "" &
             df$gene_symbol != "NA", ]
  
  # Collapse: for genes with multiple probes, keep the one with smallest P.Value
  df <- df[order(df$P.Value), ]
  df <- df[!duplicated(df$gene_symbol), ]
  
  # Ranking metric: sign(logFC) * -log10(P.Value)
  # Positive = upregulated in resistant, negative = downregulated
  df$RankMetric <- sign(df$logFC) * (-log10(df$P.Value))
  
  gene_stats_list[[gse_id]] <- df[, c("gene_symbol", "RankMetric", "logFC", "P.Value")]
  cat(sprintf("  %s: %d genes with ranking metric\n", gse_id, nrow(df)))
}

# Combine: meta-ranking = mean of RankMetric across datasets
# Only genes present in at least 2 datasets get a meta-rank
all_genes_across <- unique(unlist(lapply(gene_stats_list, `[[`, "gene_symbol")))

meta_rank <- data.frame(
  gene_symbol = all_genes_across,
  stringsAsFactors = FALSE
)

for (gse_id in names(gene_stats_list)) {
  tmp <- gene_stats_list[[gse_id]][, c("gene_symbol", "RankMetric")]
  colnames(tmp)[2] <- paste0("RankMetric_", gse_id)
  meta_rank <- merge(meta_rank, tmp, by = "gene_symbol", all = TRUE)
}

meta_rank$N_datasets <- rowSums(!is.na(meta_rank[, grep("^RankMetric_", names(meta_rank))]))
meta_rank$MetaRankMetric <- rowMeans(meta_rank[, grep("^RankMetric_", names(meta_rank))], na.rm = TRUE)

# Filter: keep genes present in >= 2 datasets
meta_rank_filt <- meta_rank[meta_rank$N_datasets >= 2, ]
meta_rank_filt <- meta_rank_filt[order(meta_rank_filt$MetaRankMetric, decreasing = TRUE), ]

cat(sprintf("\n  Genes in >= 2 datasets: %d\n", nrow(meta_rank_filt)))
cat(sprintf("  Top 5 (resistant-up):  %s\n",
            paste(head(meta_rank_filt$gene_symbol[meta_rank_filt$MetaRankMetric > 0], 5),
                  collapse = ", ")))
cat(sprintf("  Bottom 5 (resistant-down): %s\n",
            paste(tail(meta_rank_filt$gene_symbol[meta_rank_filt$MetaRankMetric < 0], 5),
                  collapse = ", ")))

# Save meta-ranking
write.csv(meta_rank_filt, file.path(ENRICH_DIR, "GSEA_meta_ranking.csv"), row.names = FALSE)

# -- 6b. Map to Entrez for GSEA --
cat("\n--- [6b] Mapping ranking to Entrez IDs ---\n")

rank_genes <- meta_rank_filt$gene_symbol
rank_map <- tryCatch({
  bitr(rank_genes, fromType = "SYMBOL", toType = "ENTREZID",
       OrgDb = org.Hs.eg.db)
}, error = function(e) NULL)

if (!is.null(rank_map)) {
  # Build sorted named vector for GSEA
  rank_df <- merge(meta_rank_filt, rank_map, by.x = "gene_symbol", by.y = "SYMBOL")
  rank_df <- rank_df[!duplicated(rank_df$ENTREZID), ]
  rank_vec <- rank_df$MetaRankMetric
  names(rank_vec) <- rank_df$ENTREZID
  rank_vec <- sort(rank_vec, decreasing = TRUE)
  
  cat(sprintf("  Ranked genes with Entrez: %d\n", length(rank_vec)))
  cat(sprintf("  Range: [%.3f, %.3f]\n", min(rank_vec), max(rank_vec)))
  cat("\n")
  
  # -- 6c. GSEA: GO --
  cat("--- [6c] GSEA: GO (BP / CC / MF) ---\n\n")
  
  gsea_go_results <- list()
  
  for (ont in go_ontologies) {
    cat(sprintf("  Running GSEA GO-%s ... ", ont))
    gse <- tryCatch({
      gseGO(geneList       = rank_vec,
            OrgDb          = org.Hs.eg.db,
            keyType        = "ENTREZID",
            ont            = ont,
            minGSSize      = 10,
            maxGSSize      = 800,
            pvalueCutoff   = 0.05,
            pAdjustMethod  = "BH",
            seed           = 42)
    }, error = function(e) {
      cat(sprintf("ERROR: %s\n", conditionMessage(e)))
      return(NULL)
    })
    
    if (!is.null(gse) && nrow(gse) > 0) {
      gsea_go_results[[ont]] <- gse
      n_up   <- sum(gse$NES > 0 & gse$p.adjust < 0.05)
      n_down <- sum(gse$NES < 0 & gse$p.adjust < 0.05)
      cat(sprintf("%d terms (%d up, %d down)\n", nrow(gse), n_up, n_down))
      
      outfile <- file.path(ENRICH_DIR, sprintf("GSEA_GO_%s.csv", ont))
      write.csv(as.data.frame(gse), outfile, row.names = FALSE)
    } else {
      gsea_go_results[[ont]] <- NULL
      cat("No significant terms\n")
    }
  }
  
  # -- 6d. GSEA: KEGG --
  cat("\n--- [6d] GSEA: KEGG ---\n")
  
  gsea_kegg <- tryCatch({
    gseKEGG(geneList       = rank_vec,
            organism        = "hsa",
            keyType         = "ncbi-geneid",
            minGSSize       = 10,
            maxGSSize       = 800,
            pvalueCutoff    = 0.05,
            pAdjustMethod   = "BH",
            seed            = 42)
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", conditionMessage(e)))
    return(NULL)
  })
  
  if (!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
    cat(sprintf("  GSEA KEGG: %d pathways (%d up, %d down)\n",
                nrow(gsea_kegg),
                sum(gsea_kegg$NES > 0 & gsea_kegg$p.adjust < 0.05),
                sum(gsea_kegg$NES < 0 & gsea_kegg$p.adjust < 0.05)))
    gsea_kegg_results <- gsea_kegg
    write.csv(as.data.frame(gsea_kegg),
              file.path(ENRICH_DIR, "GSEA_KEGG.csv"), row.names = FALSE)
  } else {
    gsea_kegg_results <- NULL
    cat("  GSEA KEGG: No significant pathways\n")
  }
  
  # -- 6e. GSEA: DisGeNET (DGN) --
  cat("\n--- [6e] GSEA: DisGeNET (DGN) ---\n")
  cat("  NOTE: Using gseDGN instead of gseDO (see Section 5 note)\n")
  
  gsea_dgn <- tryCatch({
    gseDGN(geneList        = rank_vec,
           minGSSize       = 10,
           maxGSSize       = 800,
           pvalueCutoff    = 0.05,
           pAdjustMethod   = "BH",
           seed            = 42)
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", conditionMessage(e)))
    return(NULL)
  })
  
  if (!is.null(gsea_dgn) && nrow(gsea_dgn) > 0) {
    cat(sprintf("  GSEA DGN: %d terms\n", nrow(gsea_dgn)))
    gsea_dgn_results <- gsea_dgn
    write.csv(as.data.frame(gsea_dgn),
              file.path(ENRICH_DIR, "GSEA_DGN.csv"), row.names = FALSE)
  } else {
    gsea_dgn_results <- NULL
    cat("  GSEA DGN: No significant terms\n")
  }
  
} else {
  cat("  WARNING: Entrez mapping failed. Skipping GSEA.\n")
  gsea_go_results <- list()
  gsea_kegg_results <- NULL
  gsea_dgn_results <- NULL
}

cat("\n")

# ============================================================
# 7. Generate Plots (PDF)
# ============================================================
cat("=== [7] Generating plots ===\n\n")

# ---- 7a. ORA plots ----
save_ora_plots <- function(ego, prefix, ont_or_kegg, type = "GO") {
  if (is.null(ego) || nrow(ego) == 0) return(invisible())
  
  n_show <- min(nrow(ego), 20)
  ego_plot <- ego[1:n_show, ]
  
  plot_label <- if (type == "GO") sprintf("%s GO-%s", prefix, ont_or_kegg)
                else if (type == "KEGG") sprintf("KEGG %s", prefix)
                else sprintf("DGN %s", prefix)
  
  # Barplot
  p <- tryCatch(barplot(ego_plot, showCategory = n_show, title = plot_label),
                error = function(e) NULL)
  if (!is.null(p)) {
    ggsave(file.path(RESULTS_FIG_DIR,
                     sprintf("%s_%s_barplot.pdf", prefix, ont_or_kegg)),
           p, width = 11, height = max(6, n_show * 0.35))
  }
  
  # Dotplot
  p <- tryCatch(dotplot(ego_plot, showCategory = n_show, title = plot_label),
                error = function(e) NULL)
  if (!is.null(p)) {
    ggsave(file.path(RESULTS_FIG_DIR,
                     sprintf("%s_%s_dotplot.pdf", prefix, ont_or_kegg)),
           p, width = 11, height = max(6, n_show * 0.35))
  }
  
  # Network plot (emapplot) - only if >= 3 terms
  if (n_show >= 3) {
    p <- tryCatch(emapplot(pairwise_termsim(ego_plot), showCategory = min(n_show, 15)),
                  error = function(e) NULL)
    if (!is.null(p)) {
      ggsave(file.path(RESULTS_FIG_DIR,
                       sprintf("%s_%s_emapplot.pdf", prefix, ont_or_kegg)),
             p, width = 12, height = 10)
    }
  }
  
  # Category cnetplot - only if <= 15 terms
  if (n_show <= 15) {
    p <- tryCatch(cnetplot(ego_plot, showCategory = n_show, node_label = "category"),
                  error = function(e) NULL)
    if (!is.null(p)) {
      ggsave(file.path(RESULTS_FIG_DIR,
                       sprintf("%s_%s_cnetplot.pdf", prefix, ont_or_kegg)),
             p, width = 16, height = 14)
    }
  }
  
  # Heatmap-like upset plot
  if (n_show >= 3) {
    p <- tryCatch(heatplot(ego_plot, showCategory = n_show),
                  error = function(e) NULL)
    if (!is.null(p)) {
      ggsave(file.path(RESULTS_FIG_DIR,
                       sprintf("%s_%s_heatplot.pdf", prefix, ont_or_kegg)),
             p, width = 14, height = max(6, n_show * 0.4))
    }
  }
}

# GO ORA plots
for (nm in names(go_results)) {
  parts <- strsplit(nm, "_")[[1]]
  # Heuristic: last element is ontology, everything before is set name
  ont <- parts[length(parts)]
  set_name <- paste(parts[-length(parts)], collapse = "_")
  save_ora_plots(go_results[[nm]], set_name, ont, "GO")
}

# KEGG ORA plots
for (nm in names(kegg_results)) {
  kk <- kegg_results[[nm]]
  if (is.null(kk) || nrow(kk) == 0) next
  n_show <- min(nrow(kk), 20)
  kk_plot <- kk[1:n_show, ]
  
  p_bar <- tryCatch(barplot(kk_plot, showCategory = n_show,
                             title = sprintf("KEGG %s", nm)),
                    error = function(e) NULL)
  if (!is.null(p_bar)) {
    ggsave(file.path(RESULTS_FIG_DIR, sprintf("KEGG_%s_barplot.pdf", nm)),
           p_bar, width = 11, height = max(6, n_show * 0.35))
  }
  
  p_dot <- tryCatch(dotplot(kk_plot, showCategory = n_show,
                             title = sprintf("KEGG %s", nm)),
                    error = function(e) NULL)
  if (!is.null(p_dot)) {
    ggsave(file.path(RESULTS_FIG_DIR, sprintf("KEGG_%s_dotplot.pdf", nm)),
           p_dot, width = 11, height = max(6, n_show * 0.35))
  }
  
  if (n_show >= 3) {
    p_cnet <- tryCatch(cnetplot(kk_plot, showCategory = min(n_show, 12),
                                 node_label = "category"),
                       error = function(e) NULL)
    if (!is.null(p_cnet)) {
      ggsave(file.path(RESULTS_FIG_DIR, sprintf("KEGG_%s_cnetplot.pdf", nm)),
             p_cnet, width = 16, height = 14)
    }
  }
}

# DGN plots
for (nm in names(dgn_results)) {
  dgn_ego <- dgn_results[[nm]]
  if (is.null(dgn_ego) || nrow(dgn_ego) == 0) next
  n_show <- min(nrow(dgn_ego), 20)
  dgn_plot <- dgn_ego[1:n_show, ]
  
  p_bar <- tryCatch(barplot(dgn_plot, showCategory = n_show,
                             title = sprintf("DGN %s", nm)),
                    error = function(e) NULL)
  if (!is.null(p_bar)) {
    ggsave(file.path(RESULTS_FIG_DIR, sprintf("DGN_%s_barplot.pdf", nm)),
           p_bar, width = 11, height = max(6, n_show * 0.35))
  }
  
  p_dot <- tryCatch(dotplot(dgn_plot, showCategory = n_show,
                             title = sprintf("DGN %s", nm)),
                    error = function(e) NULL)
  if (!is.null(p_dot)) {
    ggsave(file.path(RESULTS_FIG_DIR, sprintf("DGN_%s_dotplot.pdf", nm)),
           p_dot, width = 11, height = max(6, n_show * 0.35))
  }
}

# ---- 7b. GSEA plots ----
save_gsea_plots <- function(gse_obj, prefix, ont_or_kegg, type = "GO") {
  if (is.null(gse_obj) || nrow(gse_obj) == 0) return(invisible())
  
  n_show <- min(nrow(gse_obj), 20)
  gse_plot <- gse_obj[1:n_show, ]
  label <- if (type == "GO") sprintf("GSEA GO-%s", ont_or_kegg)
           else if (type == "KEGG") sprintf("GSEA KEGG %s", prefix)
           else sprintf("GSEA DGN %s", prefix)
  
  # Ridge plot (distribution of enrichment scores)
  p_ridge <- tryCatch(ridgeplot(gse_plot, showCategory = n_show) +
                        labs(title = label) +
                        theme(plot.title = element_text(hjust = 0.5)),
                      error = function(e) NULL)
  if (!is.null(p_ridge)) {
    ggsave(file.path(GSEA_DIR, sprintf("GSEA_%s_%s_ridgeplot.pdf", prefix, ont_or_kegg)),
           p_ridge, width = 12, height = max(6, n_show * 0.35))
  }
  
  # Dotplot (similar to ORA but for GSEA)
  p_dot <- tryCatch(dotplot(gse_plot, showCategory = n_show, title = label),
                    error = function(e) NULL)
  if (!is.null(p_dot)) {
    ggsave(file.path(GSEA_DIR, sprintf("GSEA_%s_%s_dotplot.pdf", prefix, ont_or_kegg)),
           p_dot, width = 11, height = max(6, n_show * 0.35))
  }
  
  # Top 4 running score plots
  top4 <- head(gse_plot, 4)
  for (i in seq_len(nrow(top4))) {
    desc <- gsub("[ /]", "_", top4$Description[i])
    p_gs <- tryCatch(gseaplot2(gse_obj, geneSetID = rownames(top4)[i],
                                title = top4$Description[i]),
                     error = function(e) NULL)
    if (!is.null(p_gs)) {
      ggsave(file.path(GSEA_DIR,
                       sprintf("GSEA_%s_%s_running_%s.pdf",
                               prefix, ont_or_kegg, substr(desc, 1, 40))),
             p_gs, width = 10, height = 6)
    }
  }
}

# GSEA GO plots
for (ont in names(gsea_go_results)) {
  save_gsea_plots(gsea_go_results[[ont]], "Meta", ont, "GO")
}

# GSEA KEGG plots
if (!is.null(gsea_kegg_results) && nrow(gsea_kegg_results) > 0) {
  n_show <- min(nrow(gsea_kegg_results), 20)
  gk <- gsea_kegg_results[1:n_show, ]
  
  p_ridge <- tryCatch(ridgeplot(gk, showCategory = n_show) +
                        labs(title = "GSEA KEGG Pathways") +
                        theme(plot.title = element_text(hjust = 0.5)),
                      error = function(e) NULL)
  if (!is.null(p_ridge)) {
    ggsave(file.path(GSEA_DIR, "GSEA_KEGG_ridgeplot.pdf"), p_ridge,
           width = 12, height = max(6, n_show * 0.35))
  }
  
  p_dot <- tryCatch(dotplot(gk, showCategory = n_show, title = "GSEA KEGG"),
                    error = function(e) NULL)
  if (!is.null(p_dot)) {
    ggsave(file.path(GSEA_DIR, "GSEA_KEGG_dotplot.pdf"), p_dot,
           width = 11, height = max(6, n_show * 0.35))
  }
  
  # Top 4 running score plots for KEGG
  top4 <- head(gk, 4)
  for (i in seq_len(nrow(top4))) {
    desc <- gsub("[ /]", "_", top4$Description[i])
    p_gs <- tryCatch(gseaplot2(gsea_kegg_results, geneSetID = rownames(top4)[i],
                                title = top4$Description[i]),
                     error = function(e) NULL)
    if (!is.null(p_gs)) {
      ggsave(file.path(GSEA_DIR,
                       sprintf("GSEA_KEGG_running_%s.pdf", substr(desc, 1, 40))),
             p_gs, width = 10, height = 6)
    }
  }
}

# GSEA DGN plots
if (!is.null(gsea_dgn_results) && nrow(gsea_dgn_results) > 0) {
  n_show <- min(nrow(gsea_dgn_results), 20)
  gd <- gsea_dgn_results[1:n_show, ]
  
  p_dot <- tryCatch(dotplot(gd, showCategory = n_show, title = "GSEA DisGeNET (DGN)"),
                    error = function(e) NULL)
  if (!is.null(p_dot)) {
    ggsave(file.path(GSEA_DIR, "GSEA_DGN_dotplot.pdf"), p_dot,
           width = 11, height = max(6, n_show * 0.35))
  }
}

# ---- 7c. Comparative multi-panel figure: top GO-BP across key tiers ----
cat("  Creating comparative multi-panel figure...\n")

key_tiers <- c("HighPriority", "Tier2_DEG_WGCNA", "All")
comp_go_list <- list()
for (kt in key_tiers) {
  key <- paste(kt, "BP", sep = "_")
  if (!is.null(go_results[[key]]) && nrow(go_results[[key]]) > 0) {
    comp_go_list[[kt]] <- go_results[[key]]
  }
}

if (length(comp_go_list) >= 2) {
  tryCatch({
    p_comp <- dotplot(comp_go_list, showCategory = 10, 
                      title = "GO-BP: Comparison across gene sets",
                      label_format = 60) +
      theme(legend.position = "right")
    ggsave(file.path(RESULTS_FIG_DIR, "GO_BP_comparison_dotplot.pdf"),
           p_comp, width = 14, height = 10)
    cat("  GO-BP comparison dotplot saved\n")
  }, error = function(e) {
    cat("  Comparison dotplot skipped:", conditionMessage(e), "\n")
  })
}

cat("\n")

# ============================================================
# 8. Cross-table: Top terms summary
# ============================================================
cat("=== [8] Cross-table summary ===\n\n")

# Collect top 5 GO-BP terms from each significant tier
cross_go <- data.frame(
  Tier = character(),
  Ontology = character(),
  Term = character(),
  p.adjust = numeric(),
  Count = numeric(),
  stringsAsFactors = FALSE
)

for (nm in names(go_results)) {
  ego <- go_results[[nm]]
  if (is.null(ego) || nrow(ego) == 0) next
  
  top5 <- head(as.data.frame(ego), 5)
  parts <- strsplit(nm, "_")[[1]]
  ont <- parts[length(parts)]
  set_name <- paste(parts[-length(parts)], collapse = "_")
  
  cross_go <- rbind(cross_go,
                     data.frame(Tier = set_name,
                                Ontology = ont,
                                Term = top5$Description,
                                p.adjust = top5$p.adjust,
                                Count = top5$Count,
                                stringsAsFactors = FALSE))
}

if (nrow(cross_go) > 0) {
  write.csv(cross_go, file.path(ENRICH_DIR, "GO_top_terms_cross_table.csv"),
            row.names = FALSE)
  cat(sprintf("  GO cross-table: %d rows\n", nrow(cross_go)))
}

# Collect top KEGG pathways
cross_kegg <- data.frame(
  Tier = character(),
  Pathway = character(),
  p.adjust = numeric(),
  Count = numeric(),
  stringsAsFactors = FALSE
)

for (nm in names(kegg_results)) {
  kk <- kegg_results[[nm]]
  if (is.null(kk) || nrow(kk) == 0) next
  
  top5 <- head(as.data.frame(kk), 5)
  cross_kegg <- rbind(cross_kegg,
                       data.frame(Tier = nm,
                                  Pathway = top5$Description,
                                  p.adjust = top5$p.adjust,
                                  Count = top5$Count,
                                  stringsAsFactors = FALSE))
}

if (nrow(cross_kegg) > 0) {
  write.csv(cross_kegg, file.path(ENRICH_DIR, "KEGG_top_pathways_cross_table.csv"),
            row.names = FALSE)
  cat(sprintf("  KEGG cross-table: %d rows\n", nrow(cross_kegg)))
}

cat("\n")

# ============================================================
# 9. Summary Report
# ============================================================
cat("=== [9] Generating summary report ===\n\n")

report_lines <- c(
  "========================================================================",
  "Enrichment Analysis Summary",
  "XELOX Resistance Study (GO + KEGG + DO + GSEA)",
  paste("Date:", Sys.Date()),
  "========================================================================",
  "",
  sprintf("Total genes in pool: %d", nrow(gp)),
  sprintf("Mapped to Entrez ID: %d / %d", nrow(entrez_map), length(all_genes)),
  sprintf("Genes with cross-dataset meta-ranking for GSEA: %d",
          if (exists("meta_rank_filt")) nrow(meta_rank_filt) else 0),
  "",
  "--- GO ORA Results (p.adjust < 0.05) ---"
)

for (nm in names(go_results)) {
  report_lines <- c(report_lines,
                    sprintf("  %-32s: %d terms", nm,
                            if (!is.null(go_results[[nm]])) nrow(go_results[[nm]]) else 0))
}

report_lines <- c(report_lines, "",
                  "--- KEGG ORA Results (p.adjust < 0.05) ---")

for (nm in names(kegg_results)) {
  report_lines <- c(report_lines,
                    sprintf("  %-32s: %d pathways", nm,
                            if (!is.null(kegg_results[[nm]])) nrow(kegg_results[[nm]]) else 0))
}

report_lines <- c(report_lines, "",
                  "--- DGN Disease Enrichment Results (p.adjust < 0.05) ---")

for (nm in names(dgn_results)) {
  report_lines <- c(report_lines,
                    sprintf("  %-32s: %d DGN disease terms", nm,
                            if (!is.null(dgn_results[[nm]])) nrow(dgn_results[[nm]]) else 0))
}

report_lines <- c(report_lines, "",
                  "--- GSEA GO Results (p.adjust < 0.05) ---")

if (exists("gsea_go_results")) {
  for (ont in names(gsea_go_results)) {
    gse <- gsea_go_results[[ont]]
    if (!is.null(gse) && nrow(gse) > 0) {
      n_up <- sum(gse$NES > 0, na.rm = TRUE)
      n_down <- sum(gse$NES < 0, na.rm = TRUE)
      report_lines <- c(report_lines,
                        sprintf("  GSEA GO-%s: %d terms (up: %d, down: %d)", ont, nrow(gse), n_up, n_down))
    } else {
      report_lines <- c(report_lines, sprintf("  GSEA GO-%s: 0 terms", ont))
    }
  }
}

report_lines <- c(report_lines, "",
                  "--- GSEA KEGG Results ---")

if (exists("gsea_kegg_results") && !is.null(gsea_kegg_results) && nrow(gsea_kegg_results) > 0) {
  n_up <- sum(gsea_kegg_results$NES > 0, na.rm = TRUE)
  n_down <- sum(gsea_kegg_results$NES < 0, na.rm = TRUE)
  report_lines <- c(report_lines,
                    sprintf("  GSEA KEGG: %d pathways (up: %d, down: %d)",
                            nrow(gsea_kegg_results), n_up, n_down))
} else {
  report_lines <- c(report_lines, "  GSEA KEGG: 0 pathways")
}

report_lines <- c(report_lines, "",
                  "--- GSEA DGN Results ---")

if (exists("gsea_dgn_results") && !is.null(gsea_dgn_results) && nrow(gsea_dgn_results) > 0) {
  report_lines <- c(report_lines,
                    sprintf("  GSEA DGN: %d terms", nrow(gsea_dgn_results)))
} else {
  report_lines <- c(report_lines, "  GSEA DGN: 0 terms")
}

# Top terms
report_lines <- c(report_lines, "",
                  "--- Top 10 GO-BP Terms (All genes) ---")

if (!is.null(go_results[["All_BP"]]) && nrow(go_results[["All_BP"]]) > 0) {
  top_go <- head(go_results[["All_BP"]], 10)
  for (i in seq_len(nrow(top_go))) {
    report_lines <- c(report_lines,
                      sprintf("  %s | p.adj=%.2e | Count=%d",
                              top_go$Description[i], top_go$p.adjust[i], top_go$Count[i]))
  }
}

report_lines <- c(report_lines, "",
                  "--- Top 10 KEGG Pathways (All) ---")

if (!is.null(kegg_results[["All"]]) && nrow(kegg_results[["All"]]) > 0) {
  top_k <- head(kegg_results[["All"]], 10)
  for (i in seq_len(nrow(top_k))) {
    report_lines <- c(report_lines,
                      sprintf("  %s | p.adj=%.2e | Count=%d",
                              top_k$Description[i], top_k$p.adjust[i], top_k$Count[i]))
  }
}

report_lines <- c(report_lines, "",
                  "--- Top 10 GSEA GO-BP Terms (by NES) ---")

if (exists("gsea_go_results") && !is.null(gsea_go_results[["BP"]]) &&
    nrow(gsea_go_results[["BP"]]) > 0) {
  gsea_bp <- as.data.frame(gsea_go_results[["BP"]])
  gsea_bp <- gsea_bp[order(gsea_bp$NES, decreasing = TRUE), ]
  top_gsea <- head(gsea_bp, 10)
  for (i in seq_len(nrow(top_gsea))) {
    report_lines <- c(report_lines,
                      sprintf("  [%s] %s | NES=%.2f | p.adj=%.2e",
                              ifelse(top_gsea$NES[i] > 0, "UP", "DOWN"),
                              top_gsea$Description[i],
                              top_gsea$NES[i], top_gsea$p.adjust[i]))
  }
}

report_lines <- c(report_lines, "",
                  "========================================================================")

writeLines(report_lines, file.path(ENRICH_DIR, "enrichment_summary_report.txt"))
cat(paste(report_lines, collapse = "\n"), "\n")

# ============================================================
# 10. Done
# ============================================================
cat("\n============================================================\n")
cat("ENRICHMENT ANALYSIS v2 COMPLETE\n")
cat("============================================================\n")
cat(sprintf("  ORA tables:  %s\n", ENRICH_DIR))
cat(sprintf("  ORA figures: %s\n", RESULTS_FIG_DIR))
cat(sprintf("  GSEA figures:%s\n", GSEA_DIR))
cat(sprintf("  Summary:     %s\n",
            file.path(ENRICH_DIR, "enrichment_summary_report.txt")))
cat("\n")

cat("Key outputs:\n")
cat(sprintf("  GO ORA:       %d tier-ontology combinations tested\n", length(go_results)))
cat(sprintf("  KEGG ORA:     %d tier-pathway results\n", length(kegg_results)))
cat(sprintf("  DGN ORA:       %d tier-disease results\n", length(dgn_results)))
if (exists("gsea_go_results")) {
  cat(sprintf("  GSEA GO:      %d ontologies\n", length(gsea_go_results)))
}
if (exists("gsea_kegg_results") && !is.null(gsea_kegg_results)) {
  cat(sprintf("  GSEA KEGG:    %d pathways\n", nrow(gsea_kegg_results)))
}
cat("\n")
