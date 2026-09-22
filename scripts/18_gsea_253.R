# ============================================================
# 18_GSEA_253.R
# Gene Set Enrichment Analysis on 253-meta entire gene ranking
# Uses ALL 25568 genes ranked by mean_logfc
# GO + KEGG + MSigDB Hallmark + C2
# XELOX Resistance Study
# ============================================================

Sys.setenv(HOME = "C:/temp", TMPDIR = "C:/temp", TMP = "C:/temp", TEMP = "C:/temp",
           R_LIBS_USER = "C:/Rlibs",
           R_USER_CACHE_DIR = "C:/temp/Rcache")
suppressWarnings(Sys.setlocale("LC_ALL", "C"))
.libPaths(c("C:/Rlibs", .libPaths()))

TEMP_INPUT      <- "C:/temp/xelox_input/deg_meta_results.csv"
ENRICH_DIR      <- "C:/temp/xelox_output/tables"
FIG_DIR         <- "C:/temp/xelox_output/figures"

dir.create(ENRICH_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("GSEA on 253-Meta Gene Ranking (ALL 25568 genes)\n")
cat("  GO (BP/CC/MF) + KEGG + MSigDB Hallmark + C2\n")
cat("============================================================\n\n")

# ============================================================
# 0. Load packages
# ============================================================
cat("=== [0] Loading packages ===\n\n")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(ggplot2)
})

cat(sprintf("  clusterProfiler v%s\n", packageVersion("clusterProfiler")))
cat(sprintf("  enrichplot     v%s\n\n", packageVersion("enrichplot")))

# ============================================================
# 1. Load meta-analysis results and build ranked gene list
# ============================================================
cat("=== [1] Loading meta-analysis results (ALL genes) ===\n\n")

dat <- read.csv(TEMP_INPUT, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("  Total genes loaded: %d\n", nrow(dat)))

# Remove rows with NA mean_logfc or stouffer_z
dat <- dat[!is.na(dat$mean_logfc) & !is.na(dat$stouffer_z), ]
cat(sprintf("  After NA removal: %d\n", nrow(dat)))

# --- 1a. Create ranking by mean_logfc (primary) ---
#   Positive = upregulated in resistant (risk)
#   Negative = downregulated in resistant (protective)
cat("  Building ranked gene list by mean_logfc ...\n")

# Remove ties: keep unique gene_symbol
dat <- dat[!duplicated(dat$gene_symbol), ]
cat(sprintf("  Unique gene symbols: %d\n", nrow(dat)))

# Build ranked vector
rank_vec <- dat$mean_logfc
names(rank_vec) <- dat$gene_symbol
rank_vec <- sort(rank_vec, decreasing = TRUE)

cat(sprintf("  Top 5 (resistant-up):    %s\n",
            paste(names(head(rank_vec[rank_vec > 0], 5)), collapse = ", ")))
cat(sprintf("  Bottom 5 (resistant-down): %s\n",
            paste(names(tail(rank_vec[rank_vec < 0], 5)), collapse = ", ")))
cat(sprintf("  Range: [%.4f, %.4f]\n", min(rank_vec), max(rank_vec)))
cat(sprintf("  Positive (up): %d, Negative (down): %d\n",
            sum(rank_vec > 0), sum(rank_vec < 0)))
cat("\n")

# ============================================================
# 2. Map gene symbols to Entrez IDs for GSEA
# ============================================================
cat("=== [2] Mapping to Entrez IDs ===\n\n")

map <- tryCatch({
  bitr(names(rank_vec), fromType = "SYMBOL", toType = "ENTREZID",
       OrgDb = org.Hs.eg.db)
}, error = function(e) {
  cat(sprintf("  WARNING bitr failed: %s\n", conditionMessage(e)))
  clean <- grep("^[A-Za-z0-9.-]+$", names(rank_vec), value = TRUE)
  bitr(clean, fromType = "SYMBOL", toType = "ENTREZID",
       OrgDb = org.Hs.eg.db)
})

cat(sprintf("  Symbols mapped to Entrez: %d\n", nrow(map)))
cat(sprintf("  Unmapped: %d\n", length(rank_vec) - nrow(map)))

# Build ranked vector with Entrez names
rank_df <- data.frame(
  gene_symbol = names(rank_vec),
  mean_logfc  = rank_vec,
  stringsAsFactors = FALSE,
  row.names = NULL
)
rank_df <- merge(rank_df, map, by.x = "gene_symbol", by.y = "SYMBOL", all.x = FALSE)
rank_df <- rank_df[!duplicated(rank_df$ENTREZID), ]

entrez_rank <- rank_df$mean_logfc
names(entrez_rank) <- rank_df$ENTREZID
entrez_rank <- sort(entrez_rank, decreasing = TRUE)

cat(sprintf("  Ranked Entrez vector: %d genes\n", length(entrez_rank)))
cat(sprintf("  Range: [%.4f, %.4f]\n\n", min(entrez_rank), max(entrez_rank)))

# Save pre-ranked list
write.csv(rank_df, file.path(ENRICH_DIR, "GSEA_253_ranking.csv"), row.names = FALSE)

# ============================================================
# 3. Ensure clusterProfiler KEGG cache directory exists
# ============================================================
cat("=== [3] Ensuring KEGG cache ===\n\n")

cp_cache <- file.path(tools::R_user_dir("clusterProfiler", "cache"))
suppressWarnings(dir.create(cp_cache, recursive = TRUE))
cat(sprintf("  KEGG cache dir: %s\n", cp_cache))
cat(sprintf("  kegg_category.rda: %s\n",
            file.exists(file.path(cp_cache, "kegg_category.rda"))))

# ============================================================
# 4. KEGG Cache Patch (same as 17_pathway_enrichment_253.R)
# ============================================================
cat("=== [4] Loading cached KEGG data ===\n\n")

KEGG_CACHE_DIR <- "C:/temp/tcga_review"
kegg_link_file <- file.path(KEGG_CACHE_DIR, "kegg_link_hsa_pathway.txt")
kegg_list_file <- file.path(KEGG_CACHE_DIR, "kegg_list_pathway_hsa.txt")
kegg_conv_file <- file.path(KEGG_CACHE_DIR, "kegg_conv_ncbi_hsa.txt")

if (file.exists(kegg_link_file) && file.exists(kegg_list_file)) {
  tryCatch({
    link_lines <- readLines(kegg_link_file, warn = FALSE)
    list_lines <- readLines(kegg_list_file, warn = FALSE)
    
    link_parts <- strsplit(link_lines, "\t")
    link_mat <- do.call(rbind, link_parts)
    kegg_extid <- data.frame(
      from = gsub("^[^:]+:", "", link_mat[, 1]),
      to   = gsub("^[^:]+:", "", link_mat[, 2]),
      stringsAsFactors = FALSE
    )
    
    list_parts <- strsplit(list_lines, "\t")
    list_mat <- do.call(rbind, list_parts)
    kegg_names <- data.frame(
      from = list_mat[, 1],
      to   = sub("\\s-\\s[^-]+$", "", list_mat[, 2]),
      stringsAsFactors = FALSE
    )
    
    KEGG_Env <- clusterProfiler:::get_KEGG_Env()
    assign("organism",    "hsa",       envir = KEGG_Env)
    assign("_type_",      "KEGG",      envir = KEGG_Env)
    assign("KEGGPATHID2EXTID", kegg_extid, envir = KEGG_Env)
    assign("KEGGPATHID2NAME",  kegg_names, envir = KEGG_Env)
    
    if (file.exists(kegg_conv_file)) {
      conv_lines <- readLines(kegg_conv_file, warn = FALSE)
      conv_parts <- strsplit(conv_lines, "\t")
      conv_mat <- do.call(rbind, conv_parts)
      idconv <- data.frame(
        from = gsub("^[^:]+:", "", conv_mat[, 1]),
        to   = gsub("^[^:]+:", "", conv_mat[, 2]),
        stringsAsFactors = FALSE
      )
      assign("key",    "ncbi-geneid", envir = KEGG_Env)
      assign("idconv", idconv,        envir = KEGG_Env)
      cat(sprintf("  KEGG idconv: %d rows\n", nrow(idconv)))
    }
    
    cat(sprintf("  KEGGPATHID2EXTID: %d rows\n", nrow(kegg_extid)))
    cat(sprintf("  KEGGPATHID2NAME:  %d rows\n", nrow(kegg_names)))
    cat("  [OK] KEGG cache populated\n\n")
  }, error = function(e) {
    cat(sprintf("  [WARN] KEGG cache patch failed: %s\n", conditionMessage(e)))
  })
} else {
  cat("  [WARN] KEGG cached files not found\n\n")
}

# ============================================================
# 5. GO GSEA (BP / CC / MF)
# ============================================================
cat("=== [5] GSEA: GO (BP / CC / MF) ===\n\n")

go_ontologies <- c("BP", "CC", "MF")
gsea_go_results <- list()

for (ont in go_ontologies) {
  cat(sprintf("  GO-%s ... ", ont))
  gse <- tryCatch({
    gseGO(geneList       = entrez_rank,
          OrgDb          = org.Hs.eg.db,
          keyType        = "ENTREZID",
          ont            = ont,
          minGSSize      = 10,
          maxGSSize      = 500,
          pvalueCutoff   = 0.05,
          pAdjustMethod  = "BH",
          seed           = 42,
          verbose        = FALSE)
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", conditionMessage(e)))
    return(NULL)
  })
  
  if (!is.null(gse) && nrow(gse) > 0) {
    gsea_go_results[[ont]] <- gse
    n_up   <- sum(gse@result$NES > 0 & gse@result$p.adjust < 0.05, na.rm = TRUE)
    n_down <- sum(gse@result$NES < 0 & gse@result$p.adjust < 0.05, na.rm = TRUE)
    cat(sprintf("%d terms (%d up, %d down)\n", nrow(gse), n_up, n_down))
    
    csv_out <- file.path(ENRICH_DIR, sprintf("GSEA_GO_%s.csv", ont))
    write.csv(as.data.frame(gse), csv_out, row.names = FALSE)
  } else {
    gsea_go_results[[ont]] <- NULL
    cat("0 terms\n")
  }
}
cat("\n")

# ============================================================
# 6. KEGG GSEA
# ============================================================
cat("=== [6] GSEA: KEGG ===\n\n")

gsea_kegg <- tryCatch({
  gseKEGG(geneList       = entrez_rank,
          organism        = "hsa",
          keyType         = "ncbi-geneid",
          minGSSize       = 10,
          maxGSSize       = 500,
          pvalueCutoff    = 0.05,
          pAdjustMethod   = "BH",
          seed            = 42,
          verbose         = FALSE)
}, error = function(e) {
  cat(sprintf("  ERROR: %s\n", conditionMessage(e)))
  return(NULL)
})

if (!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
  n_up   <- sum(gsea_kegg@result$NES > 0 & gsea_kegg@result$p.adjust < 0.05, na.rm = TRUE)
  n_down <- sum(gsea_kegg@result$NES < 0 & gsea_kegg@result$p.adjust < 0.05, na.rm = TRUE)
  cat(sprintf("  GSEA KEGG: %d pathways (%d up, %d down)\n",
              nrow(gsea_kegg), n_up, n_down))
  
  csv_out <- file.path(ENRICH_DIR, "GSEA_KEGG.csv")
  write.csv(as.data.frame(gsea_kegg), csv_out, row.names = FALSE)
} else {
  gsea_kegg <- NULL
  cat("  GSEA KEGG: 0 pathways\n")
}
cat("\n")

# ============================================================
# Helper: GMT file parser
# ============================================================
read_gmt_simple <- function(file) {
  lines <- readLines(file, warn = FALSE)
  parts <- strsplit(lines, "\t")
  term <- rep(vapply(parts, `[`, character(1), 1),
              times = vapply(parts, length, integer(1)) - 2)
  gene <- unlist(lapply(parts, function(p) p[-(1:2)]))
  data.frame(term = term, gene = gene, stringsAsFactors = FALSE)
}

# ============================================================
# 7. MSigDB Hallmark GSEA
# ============================================================
cat("=== [7] GSEA: MSigDB Hallmark ===\n\n")

GMT_HALLMARK <- "C:/temp/h.all.v2024.1.symbols.gmt"
h_gmt <- tryCatch({
  read_gmt_simple(GMT_HALLMARK)
}, error = function(e) {
  cat(sprintf("  ERROR loading Hallmark GMT: %s\n", conditionMessage(e)))
  return(NULL)
})

if (is.null(h_gmt) || nrow(h_gmt) == 0) {
  cat("  SKIP Hallmark: failed to load GMT\n")
  gsea_h_all <- NULL
} else {
  h_symbols <- unique(h_gmt$gene)
  h_map <- tryCatch({
    bitr(h_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    data.frame(SYMBOL = character(0), ENTREZID = character(0), stringsAsFactors = FALSE)
  })
  h_gmt_merged <- merge(h_gmt, h_map, by.x = "gene", by.y = "SYMBOL", all.y = FALSE)
  h_gmt_entrez <- data.frame(
    term = h_gmt_merged$term,
    gene = h_gmt_merged$ENTREZID,
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("  GMT: %d sets, %d symbols -> %d Entrez IDs\n",
              length(unique(h_gmt$term)), length(h_symbols),
              length(unique(h_gmt_entrez$gene))))
  
  # GSEA requires TERM2GENE but NOT TERM2NAME for GSEA()
  cat("  Running GSEA Hallmark ... ")
  gsea_h_all <- tryCatch({
    GSEA(geneList      = entrez_rank,
         TERM2GENE     = h_gmt_entrez,
         minGSSize     = 10,
         maxGSSize     = 500,
         pvalueCutoff  = 0.05,
         pAdjustMethod = "BH",
         seed          = 42,
         verbose       = FALSE)
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", conditionMessage(e)))
    return(NULL)
  })
  
  if (!is.null(gsea_h_all) && nrow(gsea_h_all) > 0) {
    n_up   <- sum(gsea_h_all@result$NES > 0 & gsea_h_all@result$p.adjust < 0.05, na.rm = TRUE)
    n_down <- sum(gsea_h_all@result$NES < 0 & gsea_h_all@result$p.adjust < 0.05, na.rm = TRUE)
    cat(sprintf("%d terms (%d up, %d down)\n", nrow(gsea_h_all), n_up, n_down))
    
    csv_out <- file.path(ENRICH_DIR, "GSEA_Hallmark.csv")
    write.csv(as.data.frame(gsea_h_all), csv_out, row.names = FALSE)
  } else {
    cat("0 terms\n")
  }
}
cat("\n")

# ============================================================
# 8. MSigDB C2 (KEGG Legacy) GSEA
# ============================================================
cat("=== [8] GSEA: MSigDB C2 (CP:KEGG_legacy) ===\n\n")

GMT_C2 <- "C:/temp/c2.cp.kegg_legacy.v2024.1.symbols.gmt"
c2_gmt <- tryCatch({
  read_gmt_simple(GMT_C2)
}, error = function(e) {
  cat(sprintf("  ERROR loading C2 GMT: %s\n", conditionMessage(e)))
  return(NULL)
})

if (is.null(c2_gmt) || nrow(c2_gmt) == 0) {
  cat("  SKIP C2: failed to load GMT\n")
  gsea_c2_all <- NULL
} else {
  c2_symbols <- unique(c2_gmt$gene)
  c2_map <- tryCatch({
    bitr(c2_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    data.frame(SYMBOL = character(0), ENTREZID = character(0), stringsAsFactors = FALSE)
  })
  c2_merged <- merge(c2_gmt, c2_map, by.x = "gene", by.y = "SYMBOL", all.y = FALSE)
  c2_entrez <- data.frame(
    term = c2_merged$term,
    gene = c2_merged$ENTREZID,
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("  GMT: %d sets, %d symbols -> %d Entrez IDs\n",
              length(unique(c2_gmt$term)), length(c2_symbols),
              length(unique(c2_entrez$gene))))
  
  cat("  Running GSEA C2 ... ")
  gsea_c2_all <- tryCatch({
    GSEA(geneList      = entrez_rank,
         TERM2GENE     = c2_entrez,
         minGSSize     = 10,
         maxGSSize     = 500,
         pvalueCutoff  = 0.05,
         pAdjustMethod = "BH",
         seed          = 42,
         verbose       = FALSE)
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", conditionMessage(e)))
    return(NULL)
  })
  
  if (!is.null(gsea_c2_all) && nrow(gsea_c2_all) > 0) {
    n_up   <- sum(gsea_c2_all@result$NES > 0 & gsea_c2_all@result$p.adjust < 0.05, na.rm = TRUE)
    n_down <- sum(gsea_c2_all@result$NES < 0 & gsea_c2_all@result$p.adjust < 0.05, na.rm = TRUE)
    cat(sprintf("%d terms (%d up, %d down)\n", nrow(gsea_c2_all), n_up, n_down))
    
    csv_out <- file.path(ENRICH_DIR, "GSEA_C2.csv")
    write.csv(as.data.frame(gsea_c2_all), csv_out, row.names = FALSE)
  } else {
    cat("0 terms\n")
  }
}
cat("\n")

# ============================================================
# 9. Generate GSEA plots
# ============================================================
cat("=== [9] Generating GSEA plots ===\n\n")

# --- Helper: ridgeplot for GSEA results ---
save_gsea_ridgeplot <- function(gse_obj, label) {
  if (is.null(gse_obj) || nrow(gse_obj) == 0) return(invisible())
  n_show <- min(nrow(gse_obj), 20)
  
  # S4-safe subset: use @result slot
  gse_plot <- gse_obj
  gse_plot@result <- gse_plot@result[1:n_show, , drop = FALSE]
  
  cat(sprintf("    %s ridgeplot ... ", label))
  p <- tryCatch({
    ridgeplot(gse_plot, showCategory = n_show) +
      labs(title = paste("GSEA", label)) +
      theme_bw(base_size = 11)
  }, error = function(e) { cat(sprintf("FAIL (%s)\n", conditionMessage(e))); NULL })
  
  if (!is.null(p)) {
    f <- file.path(FIG_DIR, sprintf("GSEA_%s_ridgeplot.pdf", label))
    pdf(f, width = 14, height = max(6, n_show * 0.35))
    print(p)
    dev.off()
    cat("OK\n")
  }
}

# --- Helper: dotplot for GSEA results ---
save_gsea_dotplot <- function(gse_obj, label) {
  if (is.null(gse_obj) || nrow(gse_obj) == 0) return(invisible())
  n_show <- min(nrow(gse_obj), 20)
  
  gse_plot <- gse_obj
  gse_plot@result <- gse_plot@result[1:n_show, , drop = FALSE]
  
  cat(sprintf("    %s dotplot ... ", label))
  p <- tryCatch({
    dotplot(gse_plot, showCategory = n_show) +
      labs(title = paste("GSEA", label)) +
      theme_bw(base_size = 11)
  }, error = function(e) { cat(sprintf("FAIL (%s)\n", conditionMessage(e))); NULL })
  
  if (!is.null(p)) {
    f <- file.path(FIG_DIR, sprintf("GSEA_%s_dotplot.pdf", label))
    pdf(f, width = 14, height = max(6, n_show * 0.35))
    print(p)
    dev.off()
    cat("OK\n")
  }
}

# --- Helper: gseaplot for top pathways ---
save_gsea_top_curves <- function(gse_obj, label, n_top = 5) {
  if (is.null(gse_obj) || nrow(gse_obj) == 0) return(invisible())
  
  n_top <- min(nrow(gse_obj), n_top)
  # Use @result to safely get gene set IDs
  top_ids <- head(gse_obj@result$ID, n_top)
  top_desc <- head(gse_obj@result$Description, n_top)
  
  cat(sprintf("    %s gseaplot (%d top): ", label, n_top))
  for (i in seq_along(top_ids)) {
    gs_id <- top_ids[i]
    # Use safe sub-description
    safe_desc <- gsub("[^A-Za-z0-9_-]", "", top_desc[i])
    safe_desc <- substr(safe_desc, 1, 40)
    
    p <- tryCatch({
      gseaplot(gse_obj, geneSetID = gs_id, title = top_desc[i])
    }, error = function(e) NULL)
    
    if (!is.null(p)) {
      f <- file.path(FIG_DIR, sprintf("GSEA_%s_curve_%s.pdf", label, safe_desc))
      pdf(f, width = 10, height = 7)
      print(p)
      dev.off()
      cat(".")
    }
  }
  cat(sprintf(" %d done\n", sum(sapply(top_ids, function(id) {
    !is.null(tryCatch(gseaplot(gse_obj, geneSetID = id, title = ""), error = function(e) NULL))
  }))))
}

# --- Generate all GSEA plots ---

# GO plots
for (ont in names(gsea_go_results)) {
  gse <- gsea_go_results[[ont]]
  label <- paste0("GO_", ont)
  save_gsea_ridgeplot(gse, label)
  save_gsea_dotplot(gse, label)
  save_gsea_top_curves(gse, label, n_top = 5)
}

# KEGG plots
if (!is.null(gsea_kegg)) {
  save_gsea_ridgeplot(gsea_kegg, "KEGG")
  save_gsea_dotplot(gsea_kegg, "KEGG")
  save_gsea_top_curves(gsea_kegg, "KEGG", n_top = 5)
}

# Hallmark plots
if (!is.null(gsea_h_all)) {
  save_gsea_ridgeplot(gsea_h_all, "Hallmark")
  save_gsea_dotplot(gsea_h_all, "Hallmark")
  save_gsea_top_curves(gsea_h_all, "Hallmark", n_top = 5)
}

# C2 plots
if (!is.null(gsea_c2_all)) {
  save_gsea_ridgeplot(gsea_c2_all, "C2")
  save_gsea_dotplot(gsea_c2_all, "C2")
  save_gsea_top_curves(gsea_c2_all, "C2", n_top = 5)
}

cat("\n")

# ============================================================
# 10. Summary Report
# ============================================================
cat("=== [10] Generating summary report ===\n\n")

report <- c(
  "==========================================================================",
  "GSEA on 253-Meta Gene Ranking",
  paste("Date:", Sys.Date()),
  "==========================================================================",
  "",
  sprintf("Total genes ranked: %d", length(entrez_rank)),
  sprintf("Positive (up in resistant): %d", sum(entrez_rank > 0)),
  sprintf("Negative (down in resistant): %d", sum(entrez_rank < 0)),
  sprintf("Ranking metric: mean_logfc"),
  sprintf("Range: [%.4f, %.4f]", min(entrez_rank), max(entrez_rank)),
  "",
  "--- GSEA GO ---"
)

for (ont in names(gsea_go_results)) {
  r <- gsea_go_results[[ont]]
  n <- if (!is.null(r)) nrow(r) else 0
  n_up <- if (!is.null(r)) sum(r@result$NES > 0 & r@result$p.adjust < 0.05, na.rm = TRUE) else 0
  n_down <- if (!is.null(r)) sum(r@result$NES < 0 & r@result$p.adjust < 0.05, na.rm = TRUE) else 0
  report <- c(report, sprintf("  GO-%-8s: %d terms (%d up, %d down)", ont, n, n_up, n_down))
}

report <- c(report, "",
            "--- GSEA KEGG ---")
if (!is.null(gsea_kegg)) {
  n_up <- sum(gsea_kegg@result$NES > 0 & gsea_kegg@result$p.adjust < 0.05, na.rm = TRUE)
  n_down <- sum(gsea_kegg@result$NES < 0 & gsea_kegg@result$p.adjust < 0.05, na.rm = TRUE)
  report <- c(report, sprintf("  KEGG: %d pathways (%d up, %d down)", nrow(gsea_kegg), n_up, n_down))
}

report <- c(report, "",
            "--- GSEA Hallmark ---")
if (!is.null(gsea_h_all)) {
  n_up <- sum(gsea_h_all@result$NES > 0 & gsea_h_all@result$p.adjust < 0.05, na.rm = TRUE)
  n_down <- sum(gsea_h_all@result$NES < 0 & gsea_h_all@result$p.adjust < 0.05, na.rm = TRUE)
  report <- c(report, sprintf("  Hallmark: %d terms (%d up, %d down)", nrow(gsea_h_all), n_up, n_down))
}

report <- c(report, "",
            "--- GSEA C2 (KEGG_legacy) ---")
if (!is.null(gsea_c2_all)) {
  n_up <- sum(gsea_c2_all@result$NES > 0 & gsea_c2_all@result$p.adjust < 0.05, na.rm = TRUE)
  n_down <- sum(gsea_c2_all@result$NES < 0 & gsea_c2_all@result$p.adjust < 0.05, na.rm = TRUE)
  report <- c(report, sprintf("  C2: %d terms (%d up, %d down)", nrow(gsea_c2_all), n_up, n_down))
}

# Add top GO-BP terms
report <- c(report, "",
            "--- Top 10 GSEA GO-BP ---")
if (!is.null(gsea_go_results[["BP"]]) && nrow(gsea_go_results[["BP"]]) > 0) {
  top10 <- head(as.data.frame(gsea_go_results[["BP"]]), 10)
  for (i in seq_len(nrow(top10))) {
    report <- c(report, sprintf("  %s | NES=%.3f | p.adj=%.2e | setSize=%d",
                                top10$Description[i], top10$NES[i],
                                top10$p.adjust[i], top10$setSize[i]))
  }
}

report <- c(report, "",
            "--- Top 10 GSEA KEGG ---")
if (!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
  top10 <- head(as.data.frame(gsea_kegg), 10)
  for (i in seq_len(nrow(top10))) {
    report <- c(report, sprintf("  %s | NES=%.3f | p.adj=%.2e | setSize=%d",
                                top10$Description[i], top10$NES[i],
                                top10$p.adjust[i], top10$setSize[i]))
  }
}

report <- c(report, "",
            "--- Top 10 GSEA Hallmark ---")
if (!is.null(gsea_h_all) && nrow(gsea_h_all) > 0) {
  top10 <- head(as.data.frame(gsea_h_all), 10)
  for (i in seq_len(nrow(top10))) {
    report <- c(report, sprintf("  %s | NES=%.3f | p.adj=%.2e | setSize=%d",
                                top10$Description[i], top10$NES[i],
                                top10$p.adjust[i], top10$setSize[i]))
  }
}

report <- c(report, "",
            "==========================================================================")

writeLines(report, file.path(ENRICH_DIR, "gsea_summary_253.txt"))
cat(paste(report, collapse = "\n"), "\n")

# ============================================================
# 11. Done
# ============================================================
cat("\n============================================================\n")
cat("GSEA ANALYSIS v18 COMPLETE\n")
cat("============================================================\n")
cat(sprintf("  Tables:  %s\n", ENRICH_DIR))
cat(sprintf("  Figures: %s\n", FIG_DIR))
cat(sprintf("  Summary: %s\n", file.path(ENRICH_DIR, "gsea_summary_253.txt")))
cat("\n")
