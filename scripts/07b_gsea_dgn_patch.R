# ============================================================
# 07b_GSEA_DGN_PATCH.R
# Patch: runs GSEA DGN + plots + summary for XELOX enrichment
# Skips already-completed ORA and GSEA (GO/KEGG) steps.
# ============================================================

Sys.setenv(HOME = "C:/temp", TMPDIR = "C:/temp", TMP = "C:/temp", TEMP = "C:/temp", R_LIBS = "C:/Rlibs")
Sys.setlocale("LC_ALL", "C")
.libPaths(c("C:/Rlibs", .libPaths()))

PROJECT_ROOT     <- "C:/xelox_research"
SCRIPT_DIR       <- file.path(PROJECT_ROOT, "scripts")
RESULTS_TAB_DIR  <- file.path(PROJECT_ROOT, "results", "tables")
RESULTS_FIG_DIR  <- file.path(PROJECT_ROOT, "results", "figures")
ENRICH_DIR       <- file.path(RESULTS_TAB_DIR, "enrichment")
GSEA_DIR         <- file.path(RESULTS_FIG_DIR, "gsea")

cat("=== Loading packages ===\n")
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(DOSE)
library(ggplot2)

cat("=== Loading existing GSEA meta-ranking ===\n")
meta_rank_filt <- read.csv(file.path(ENRICH_DIR, "GSEA_meta_ranking.csv"),
                           stringsAsFactors = FALSE)

rank_genes <- meta_rank_filt$gene_symbol
rank_map <- bitr(rank_genes, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)

rank_df <- merge(meta_rank_filt, rank_map, by.x = "gene_symbol", by.y = "SYMBOL")
rank_df <- rank_df[!duplicated(rank_df$ENTREZID), ]
rank_vec <- rank_df$MetaRankMetric
names(rank_vec) <- rank_df$ENTREZID
rank_vec <- sort(rank_vec, decreasing = TRUE)

# Store for later use
saveRDS(rank_vec, file.path(ENRICH_DIR, "gsea_rank_vec.rds"))

cat(sprintf("  Ranked genes with Entrez: %d\n", length(rank_vec)))
cat(sprintf("  Range: [%.3f, %.3f]\n", min(rank_vec), max(rank_vec)))

# ============================================================
# 6e. GSEA: DisGeNET (DGN) - THIS IS THE CRITICAL PART
# ============================================================
cat("\n=== [6e] GSEA: DisGeNET (DGN) ===\n")
cat("  Using gseDGN with R_LIBS fix for snow dependency\n")

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

# ============================================================
# Load existing results for plotting & summary
# ============================================================
cat("\n=== Loading existing enrichment results ===\n")

# GO results
go_files <- list.files(ENRICH_DIR, pattern = "^GO_.+\\.csv$")
go_results <- list()
for (f in go_files) {
  nm <- gsub("\\.csv$", "", f)
  nm <- gsub("^GO_", "", nm)
  df <- read.csv(file.path(ENRICH_DIR, f), stringsAsFactors = FALSE)
  if (nrow(df) > 0) go_results[[nm]] <- df
}
cat(sprintf("  Loaded %d GO result sets\n", length(go_results)))

# KEGG results
kegg_files <- list.files(ENRICH_DIR, pattern = "^KEGG_.+\\.csv$")
kegg_results <- list()
for (f in kegg_files) {
  nm <- gsub("\\.csv$", "", f)
  nm <- gsub("^KEGG_", "", nm)
  df <- read.csv(file.path(ENRICH_DIR, f), stringsAsFactors = FALSE)
  if (nrow(df) > 0) kegg_results[[nm]] <- df
}
cat(sprintf("  Loaded %d KEGG result sets\n", length(kegg_results)))

# DGN results
dgn_files <- list.files(ENRICH_DIR, pattern = "^DGN_.+\\.csv$")
dgn_results <- list()
for (f in dgn_files) {
  nm <- gsub("\\.csv$", "", f)
  nm <- gsub("^DGN_", "", nm)
  df <- read.csv(file.path(ENRICH_DIR, f), stringsAsFactors = FALSE)
  if (nrow(df) > 0) dgn_results[[nm]] <- df
}
cat(sprintf("  Loaded %d DGN result sets\n", length(dgn_results)))

# GSEA GO
gsea_go_results <- list()
for (ont in c("BP", "CC", "MF")) {
  f <- file.path(ENRICH_DIR, sprintf("GSEA_GO_%s.csv", ont))
  if (file.exists(f)) {
    df <- read.csv(f, stringsAsFactors = FALSE)
    if (nrow(df) > 0) gsea_go_results[[ont]] <- df
  }
}
cat(sprintf("  Loaded %d GSEA GO result sets\n", length(gsea_go_results)))

# GSEA KEGG
gsea_kegg_results <- NULL
f <- file.path(ENRICH_DIR, "GSEA_KEGG.csv")
if (file.exists(f)) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  if (nrow(df) > 0) gsea_kegg_results <- df
}
cat(sprintf("  GSEA KEGG loaded: %s\n", ifelse(is.null(gsea_kegg_results), "no", paste(nrow(gsea_kegg_results), "terms"))))

# ============================================================
# 7. Generate Plots (PDF)
# ============================================================
cat("\n=== [7] Generating plots ===\n\n")

# ---- 7a. ORA plots ----
save_ora_plots <- function(ego_df, prefix, ont_or_kegg, type = "GO") {
  if (is.null(ego_df) || nrow(ego_df) == 0) return(invisible())
  
  n_show <- min(nrow(ego_df), 20)
  ego_plot <- ego_df[1:n_show, ]
  
  plot_label <- if (type == "GO") sprintf("%s GO-%s", prefix, ont_or_kegg)
                else if (type == "KEGG") sprintf("KEGG %s", prefix)
                else sprintf("DGN %s", prefix)
  
  # Barplot
  p <- tryCatch({
    ego_plot$GeneRatio <- ego_plot$GeneRatio  # keep as is
    ggplot(ego_plot, aes(x = reorder(Description, Count), y = Count, fill = p.adjust)) +
      geom_bar(stat = "identity") + coord_flip() +
      labs(title = plot_label, x = "", y = "Count") +
      theme_minimal() + scale_fill_gradient(low = "red", high = "blue")
  }, error = function(e) NULL)
  if (!is.null(p)) {
    ggsave(file.path(RESULTS_FIG_DIR,
                     sprintf("%s_%s_barplot.pdf", prefix, ont_or_kegg)),
           p, width = 11, height = max(6, n_show * 0.35))
  }
  
  # Dotplot (manually)
  p <- tryCatch({
    ggplot(ego_plot, aes(x = Count, y = reorder(Description, Count),
                          size = Count, color = p.adjust)) +
      geom_point() + labs(title = plot_label, x = "Count", y = "") +
      theme_minimal() + scale_color_gradient(low = "red", high = "blue")
  }, error = function(e) NULL)
  if (!is.null(p)) {
    ggsave(file.path(RESULTS_FIG_DIR,
                     sprintf("%s_%s_dotplot.pdf", prefix, ont_or_kegg)),
           p, width = 11, height = max(6, n_show * 0.35))
  }
}

# GO ORA plots
for (nm in names(go_results)) {
  parts <- strsplit(nm, "_")[[1]]
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
  
  p_bar <- tryCatch({
    ggplot(kk_plot, aes(x = reorder(Description, Count), y = Count, fill = p.adjust)) +
      geom_bar(stat = "identity") + coord_flip() +
      labs(title = sprintf("KEGG %s", nm), x = "", y = "Count") +
      theme_minimal() + scale_fill_gradient(low = "red", high = "blue")
  }, error = function(e) NULL)
  if (!is.null(p_bar)) {
    ggsave(file.path(RESULTS_FIG_DIR, sprintf("KEGG_%s_barplot.pdf", nm)),
           p_bar, width = 11, height = max(6, n_show * 0.35))
  }
}

# DGN plots - using dgn_results variable
for (nm in names(dgn_results)) {
  dgn_df <- dgn_results[[nm]]
  if (is.null(dgn_df) || nrow(dgn_df) == 0) next
  n_show <- min(nrow(dgn_df), 20)
  dgn_plot <- dgn_df[1:n_show, ]
  
  p_bar <- tryCatch({
    ggplot(dgn_plot, aes(x = reorder(Description, Count), y = Count, fill = p.adjust)) +
      geom_bar(stat = "identity") + coord_flip() +
      labs(title = sprintf("DGN %s", nm), x = "", y = "Count") +
      theme_minimal() + scale_fill_gradient(low = "red", high = "blue")
  }, error = function(e) NULL)
  if (!is.null(p_bar)) {
    ggsave(file.path(RESULTS_FIG_DIR, sprintf("DGN_%s_barplot.pdf", nm)),
           p_bar, width = 11, height = max(6, n_show * 0.35))
  }
}

# ---- 7b. GSEA plots ----
# GSEA GO ridge plots
for (ont in names(gsea_go_results)) {
  gse_df <- gsea_go_results[[ont]]
  if (is.null(gse_df) || nrow(gse_df) == 0) next
  n_show <- min(nrow(gse_df), 20)
  gse_plot <- gse_df[1:n_show, ]
  
  p_ridge <- tryCatch({
    ggplot(gse_plot, aes(x = NES, y = reorder(Description, NES), fill = p.adjust)) +
      geom_density_ridges(stat = "identity", scale = 0.9) +
      labs(title = sprintf("GSEA GO-%s", ont), x = "NES", y = "") +
      theme_minimal() + scale_fill_gradient(low = "red", high = "blue")
  }, error = function(e) NULL)
  if (!is.null(p_ridge)) {
    ggsave(file.path(GSEA_DIR, sprintf("GSEA_Meta_%s_ridgeplot.pdf", ont)),
           p_ridge, width = 12, height = max(6, n_show * 0.35))
  }
}

# GSEA KEGG
if (!is.null(gsea_kegg_results) && nrow(gsea_kegg_results) > 0) {
  n_show <- min(nrow(gsea_kegg_results), 20)
  gk <- gsea_kegg_results[1:n_show, ]
  
  p_ridge <- tryCatch({
    ggplot(gk, aes(x = NES, y = reorder(Description, NES), fill = p.adjust)) +
      geom_density_ridges(stat = "identity", scale = 0.9) +
      labs(title = "GSEA KEGG Pathways", x = "NES", y = "") +
      theme_minimal() + scale_fill_gradient(low = "red", high = "blue")
  }, error = function(e) NULL)
  if (!is.null(p_ridge)) {
    ggsave(file.path(GSEA_DIR, "GSEA_KEGG_ridgeplot.pdf"),
           p_ridge, width = 12, height = max(6, n_show * 0.35))
  }
}

# GSEA DGN plots
if (!is.null(gsea_dgn_results) && nrow(gsea_dgn_results) > 0) {
  n_show <- min(nrow(gsea_dgn_results), 20)
  gd <- gsea_dgn_results[1:n_show, ]
  
  p_dot <- tryCatch({
    ggplot(gd, aes(x = NES, y = reorder(Description, NES), size = setSize, color = p.adjust)) +
      geom_point() + labs(title = "GSEA DisGeNET (DGN)", x = "NES", y = "") +
      theme_minimal() + scale_color_gradient(low = "red", high = "blue")
  }, error = function(e) NULL)
  if (!is.null(p_dot)) {
    ggsave(file.path(GSEA_DIR, "GSEA_DGN_dotplot.pdf"),
           p_dot, width = 11, height = max(6, n_show * 0.35))
  }
}

cat("  All plots generated\n")

# ============================================================
# 8. Cross-table: Top terms summary
# ============================================================
cat("\n=== [8] Cross-table summary ===\n")

cross_go <- data.frame(Tier = character(), Ontology = character(),
                        Term = character(), p.adjust = numeric(),
                        Count = numeric(), stringsAsFactors = FALSE)
for (nm in names(go_results)) {
  df <- go_results[[nm]]
  if (is.null(df) || nrow(df) == 0) next
  top5 <- head(df, 5)
  parts <- strsplit(nm, "_")[[1]]
  ont <- parts[length(parts)]
  set_name <- paste(parts[-length(parts)], collapse = "_")
  cross_go <- rbind(cross_go,
                     data.frame(Tier = set_name, Ontology = ont,
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

cross_kegg <- data.frame(Tier = character(), Pathway = character(),
                          p.adjust = numeric(), Count = numeric(),
                          stringsAsFactors = FALSE)
for (nm in names(kegg_results)) {
  df <- kegg_results[[nm]]
  if (is.null(df) || nrow(df) == 0) next
  top5 <- head(df, 5)
  cross_kegg <- rbind(cross_kegg,
                       data.frame(Tier = nm, Pathway = top5$Description,
                                  p.adjust = top5$p.adjust,
                                  Count = top5$Count,
                                  stringsAsFactors = FALSE))
}
if (nrow(cross_kegg) > 0) {
  write.csv(cross_kegg, file.path(ENRICH_DIR, "KEGG_top_pathways_cross_table.csv"),
            row.names = FALSE)
  cat(sprintf("  KEGG cross-table: %d rows\n", nrow(cross_kegg)))
}

# ============================================================
# 9. Summary Report
# ============================================================
cat("\n=== [9] Generating summary report ===\n")

# Count GO results
go_count <- sum(sapply(go_results, nrow))
kegg_count <- sum(sapply(kegg_results, nrow))
dgn_count <- sum(sapply(dgn_results, nrow))
gsea_go_count <- sum(sapply(gsea_go_results, function(x) if (is.data.frame(x)) nrow(x) else 0))
gsea_kegg_count <- if (!is.null(gsea_kegg_results)) nrow(gsea_kegg_results) else 0
gsea_dgn_count <- if (!is.null(gsea_dgn_results)) nrow(gsea_dgn_results) else 0

report_lines <- c(
  "========================================================================",
  "Enrichment Analysis Summary",
  "XELOX Resistance Study (GO + KEGG + DGN + GSEA)",
  paste("Date:", Sys.Date()),
  "========================================================================",
  "",
  sprintf("Total ORA results:"),
  sprintf("  GO (BP/CC/MF):      %d terms across all tiers", go_count),
  sprintf("  KEGG:               %d pathways across all tiers", kegg_count),
  sprintf("  DGN (DisGeNET):     %d disease terms across all tiers", dgn_count),
  "",
  sprintf("Total GSEA results (preranked meta-analysis):"),
  sprintf("  GSEA GO:            %d terms (BP/CC/MF)", gsea_go_count),
  sprintf("  GSEA KEGG:          %d pathways", gsea_kegg_count),
  sprintf("  GSEA DGN:           %d disease terms", gsea_dgn_count),
  "",
  "--- Key GO-BP Terms (All genes, top 10) ---"
)

if (!is.null(go_results[["All_BP"]]) && nrow(go_results[["All_BP"]]) > 0) {
  for (i in 1:min(10, nrow(go_results[["All_BP"]]))) {
    report_lines <- c(report_lines,
      sprintf("  %s | p.adj=%.2e | Count=%d",
              go_results[["All_BP"]]$Description[i],
              go_results[["All_BP"]]$p.adjust[i],
              go_results[["All_BP"]]$Count[i]))
  }
}

report_lines <- c(report_lines, "",
  "--- Key KEGG Pathways (All, top 10) ---")

if (!is.null(kegg_results[["All"]]) && nrow(kegg_results[["All"]]) > 0) {
  for (i in 1:min(10, nrow(kegg_results[["All"]]))) {
    report_lines <- c(report_lines,
      sprintf("  %s | p.adj=%.2e | Count=%d",
              kegg_results[["All"]]$Description[i],
              kegg_results[["All"]]$p.adjust[i],
              kegg_results[["All"]]$Count[i]))
  }
}

report_lines <- c(report_lines, "",
  "--- Key DGN Disease Terms (All, top 10) ---")

if (!is.null(dgn_results[["All"]]) && nrow(dgn_results[["All"]]) > 0) {
  for (i in 1:min(10, nrow(dgn_results[["All"]]))) {
    report_lines <- c(report_lines,
      sprintf("  %s | p.adj=%.2e | Count=%d",
              dgn_results[["All"]]$Description[i],
              dgn_results[["All"]]$p.adjust[i],
              dgn_results[["All"]]$Count[i]))
  }
}

report_lines <- c(report_lines, "",
  "========================================================================")

writeLines(report_lines, file.path(ENRICH_DIR, "enrichment_summary_report.txt"))
cat(paste(report_lines, collapse = "\n"), "\n")

cat("\n============================================================\n")
cat("ENRICHMENT ANALYSIS v2 COMPLETE (via patch)\n")
cat("============================================================\n")
cat(sprintf("  ORA tables:  %s\n", ENRICH_DIR))
cat(sprintf("  ORA figures: %s\n", RESULTS_FIG_DIR))
cat(sprintf("  GSEA figures:%s\n", GSEA_DIR))
cat(sprintf("  Summary:     %s\n", file.path(ENRICH_DIR, "enrichment_summary_report.txt")))
cat("\n")
