# ============================================================
# 17_PATHWAY_ENRICHMENT_253.R
# Pathway enrichment on 253 meta-analysis significant genes
# Input: deg_meta_results.csv (all 25568 genes, filter stouffer_fdr < 0.05 for 253 sig)
# Output: GO/KEGG/MSigDB enrichment tables + figures
# XELOX Resistance Study
# ============================================================

Sys.setenv(HOME = "/tmp", TMPDIR = "/tmp", TMP = "/tmp", TEMP = "/tmp",
           R_LIBS_USER = "/path/to/Rlibs",
           R_USER_CACHE_DIR = "/tmp/Rcache")  # avoid Chinese chars in default cache path
# Also override msigdbr cache dir via options
options(msigdbr_cache_dir = "/tmp/Rcache/msigdbr")
# Set locale to C to avoid issues with non-ASCII paths in read.csv
suppressWarnings(Sys.setlocale("LC_ALL", "C"))
.libPaths(c("/path/to/Rlibs", .libPaths()))

# Use short (ASCII-only) temp paths for I/O
# NOTE: deg_meta_results.csv contains ALL 25568 genes from the meta-analysis.
#       We filter stouffer_fdr < 0.05 to obtain the 253 significant genes.
TEMP_INPUT      <- "/tmp/xelox_input/deg_meta_results.csv"
ENRICH_DIR      <- "/tmp/xelox_output/tables"
FIG_DIR         <- "/tmp/xelox_output/figures"

dir.create(ENRICH_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("Pathway Enrichment on 253 Meta-Analysis Significant Genes\n")
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
  library(ggrepel)
  library(msigdbr)
})

cat(sprintf("  clusterProfiler v%s\n", packageVersion("clusterProfiler")))
cat(sprintf("  enrichplot     v%s\n", packageVersion("enrichplot")))
cat(sprintf("  msigdbr        v%s\n\n", packageVersion("msigdbr")))

# ============================================================
# 1. Load 253 significant genes
# ============================================================
cat("=== [1] Loading 253 significant meta-analysis genes ===\n\n")

dat <- read.csv(TEMP_INPUT, stringsAsFactors = FALSE, check.names = FALSE)

# Keep only Stouffer FDR < 0.05
# NOTE: use [[ for column access to avoid partial matching issues with $
#       CRITICAL: must also exclude NAs — dat[[col]] < 0.05 returns NA for
#       missing values, and dat[NA, ] silently includes those rows.
fdr_col <- "stouffer_fdr"
sig <- dat[!is.na(dat[[fdr_col]]) & dat[[fdr_col]] < 0.05, ]
sig <- sig[order(sig$stouffer_fdr), ]

cat(sprintf("  Significant genes (Stouffer FDR < 0.05): %d\n", nrow(sig)))
cat(sprintf("  Consensus up:   %d\n", sum(sig$direction_consensus == "up", na.rm = TRUE)))
cat(sprintf("  Consensus down: %d\n", sum(sig$direction_consensus == "down", na.rm = TRUE)))
cat(sprintf("  Mixed/NA direction: %d\n",
            sum(!sig$direction_consensus %in% c("up", "down"), na.rm = TRUE)))
cat("\n")

# Split by direction (handle NAs in direction_consensus)
genes_up   <- sig$gene_symbol[!is.na(sig$direction_consensus) & sig$direction_consensus == "up"]
genes_down <- sig$gene_symbol[!is.na(sig$direction_consensus) & sig$direction_consensus == "down"]
genes_all  <- sig$gene_symbol

# Save clean gene lists
writeLines(genes_all,  file.path(ENRICH_DIR, "gene_list_all.txt"))
writeLines(genes_up,   file.path(ENRICH_DIR, "gene_list_up.txt"))
writeLines(genes_down, file.path(ENRICH_DIR, "gene_list_down.txt"))

cat(sprintf("  All genes:   %d\n", length(genes_all)))
cat(sprintf("  Up genes:    %d\n", length(genes_up)))
cat(sprintf("  Down genes:  %d\n", length(genes_down)))

# ============================================================
# 2. Gene symbol -> Entrez ID mapping
# ============================================================
cat("\n=== [2] Mapping to Entrez IDs ===\n\n")

gene_tiers <- list(All = genes_all, Up = genes_up, Down = genes_down)

entrez_lists <- list()
unmapped_list <- list()

for (nm in names(gene_tiers)) {
  gl <- gene_tiers[[nm]]
  if (length(gl) < 3) {
    cat(sprintf("  %-8s: SKIP (< 3 genes)\n", nm))
    entrez_lists[[nm]] <- character(0)
    unmapped_list[[nm]] <- character(0)
    next
  }
  
  map <- tryCatch({
    bitr(gl, fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    cat(sprintf("    WARNING: bitr failed: %s\n", conditionMessage(e)))
    # Fallback: try with clean symbols
    clean <- grep("^[A-Za-z0-9.-]+$", gl, value = TRUE)
    bitr(clean, fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Hs.eg.db)
  })
  
  mapped_entrez <- unique(map$ENTREZID)
  unmapped <- setdiff(gl, unique(map$SYMBOL))
  
  entrez_lists[[nm]] <- mapped_entrez
  unmapped_list[[nm]] <- unmapped
  
  cat(sprintf("  %-8s: %d symbols -> %d Entrez IDs (unmapped: %d)\n",
              nm, length(gl), length(mapped_entrez), length(unmapped)))
  if (length(unmapped) > 0) {
    cat(sprintf("           e.g.: %s\n", paste(head(unmapped, 5), collapse = ", ")))
  }
}

# Save unmapped
all_unmapped <- unique(unlist(unmapped_list))
if (length(all_unmapped) > 0) {
  writeLines(all_unmapped, file.path(ENRICH_DIR, "unmapped_symbols.txt"))
}
cat("\n")

# ============================================================
# 3. GO Enrichment (BP / CC / MF)
# ============================================================
cat("=== [3] GO ORA: BP / CC / MF ===\n\n")

go_ontologies <- c("BP", "CC", "MF")
go_results <- list()

for (nm in names(entrez_lists)) {
  genes <- entrez_lists[[nm]]
  if (length(genes) < 5) {
    cat(sprintf("  %-8s: SKIP GO (< 5 Entrez IDs)\n", nm))
    for (ont in go_ontologies) {
      go_results[[paste(nm, ont, sep = "_")]] <- NULL
    }
    next
  }
  
  for (ont in go_ontologies) {
    cat(sprintf("  %-8s GO-%s ... ", nm, ont))
    
    ego <- tryCatch({
      enrichGO(gene          = genes,
               OrgDb         = org.Hs.eg.db,
               keyType       = "ENTREZID",
               ont           = ont,
               pvalueCutoff  = 0.05,
               pAdjustMethod = "BH",
               minGSSize     = 5,
               maxGSSize     = NULL,   # CP 4.20 default maxGSSize=500 filters large BP terms
               readable      = TRUE)
    }, error = function(e) {
      cat(sprintf("ERROR: %s\n", conditionMessage(e)))
      return(NULL)
    })
    
    key <- paste(nm, ont, sep = "_")
    if (!is.null(ego) && nrow(ego) > 0) {
      go_results[[key]] <- ego
      cat(sprintf("%d terms\n", nrow(ego)))
      
      csv_out <- file.path(ENRICH_DIR, sprintf("GO_%s_%s.csv", nm, ont))
      write.csv(as.data.frame(ego), csv_out, row.names = FALSE)
    } else {
      go_results[[key]] <- NULL
      cat("0 terms\n")
    }
  }
}

cat("\n")

# ============================================================
# 3b. KEGG Cache Patch (pre-downloaded files, no internet)
# ============================================================
cat("=== [3b] Loading cached KEGG data ===\n\n")

KEGG_CACHE_DIR <- "/tmp/tcga_review"
kegg_link_file <- file.path(KEGG_CACHE_DIR, "kegg_link_hsa_pathway.txt")
kegg_list_file <- file.path(KEGG_CACHE_DIR, "kegg_list_pathway_hsa.txt")
kegg_conv_file <- file.path(KEGG_CACHE_DIR, "kegg_conv_ncbi_hsa.txt")

# Ensure clusterProfiler KEGG cache directory exists (avoids garbled-chinese-path
# error when enrichKEGG tries to auto-create via dir.create)
cp_cache <- file.path(tools::R_user_dir("clusterProfiler", "cache"))
suppressWarnings(dir.create(cp_cache, recursive = TRUE))
cat(sprintf("  clusterProfiler KEGG cache dir: %s\n", cp_cache))
cat(sprintf("  kegg_category.rda present: %s\n",
            file.exists(file.path(cp_cache, "kegg_category.rda"))))

if (file.exists(kegg_link_file) && file.exists(kegg_list_file)) {
  tryCatch({
    # Read pre-downloaded KEGG data
    link_lines <- readLines(kegg_link_file, warn = FALSE)
    list_lines <- readLines(kegg_list_file, warn = FALSE)
    
    # Parse pathway→gene mapping: "path:hsa00010\thsa:10327"
    link_parts <- strsplit(link_lines, "\t")
    link_mat <- do.call(rbind, link_parts)
    kegg_extid <- data.frame(
      from = gsub("^[^:]+:", "", link_mat[, 1]),  # "hsa00010" -> "00010"
      to   = gsub("^[^:]+:", "", link_mat[, 2]),   # "10327"
      stringsAsFactors = FALSE
    )
    
    # Parse pathway→name mapping: "hsa01100\tMetabolic pathways - H. sapiens"
    list_parts <- strsplit(list_lines, "\t")
    list_mat <- do.call(rbind, list_parts)
    kegg_names <- data.frame(
      from = list_mat[, 1],
      to   = sub("\\s-\\s[^-]+$", "", list_mat[, 2]),
      stringsAsFactors = FALSE
    )
    
    # Populate clusterProfiler KEGG cache environment
    KEGG_Env <- clusterProfiler:::get_KEGG_Env()
    assign("organism",    "hsa",       envir = KEGG_Env)
    assign("_type_",      "KEGG",      envir = KEGG_Env)
    assign("KEGGPATHID2EXTID", kegg_extid, envir = KEGG_Env)
    assign("KEGGPATHID2NAME",  kegg_names, envir = KEGG_Env)
    
    # Also pre-cache ncbi-geneid conversion (enrichKEGG keyType="ncbi-geneid" needs it)
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
      cat(sprintf("  KEGG idconv (ncbi-geneid): %d rows\n", nrow(idconv)))
    }
    
    cat(sprintf("  KEGGPATHID2EXTID: %d rows\n", nrow(kegg_extid)))
    cat(sprintf("  KEGGPATHID2NAME:  %d rows\n", nrow(kegg_names)))
    cat("  [OK] KEGG cache populated from pre-downloaded files\n\n")
  }, error = function(e) {
    cat(sprintf("  [WARN] KEGG cache patch failed: %s\n", conditionMessage(e)))
    cat("  Will try online (may timeout)\n\n")
  })
} else {
  cat("  [WARN] KEGG cached files not found, will try online\n\n")
}

# ============================================================
# 4. KEGG Pathway Enrichment
# ============================================================
cat("=== [4] KEGG ORA ===\n\n")

kegg_results <- list()

for (nm in names(entrez_lists)) {
  genes <- entrez_lists[[nm]]
  if (length(genes) < 5) {
    cat(sprintf("  %-8s: SKIP KEGG (< 5 Entrez IDs)\n", nm))
    kegg_results[[nm]] <- NULL
    next
  }
  
  cat(sprintf("  %-8s KEGG ... ", nm))
  kk <- tryCatch({
    enrichKEGG(gene          = genes,
               organism       = "hsa",
               keyType        = "ncbi-geneid",
               pvalueCutoff   = 0.05,
               pAdjustMethod  = "BH",
               minGSSize      = 5,
               maxGSSize      = NULL)
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", conditionMessage(e)))
    return(NULL)
  })
  
  if (!is.null(kk) && nrow(kk) > 0) {
    kegg_results[[nm]] <- kk
    cat(sprintf("%d pathways\n", nrow(kk)))
    
    csv_out <- file.path(ENRICH_DIR, sprintf("KEGG_%s.csv", nm))
    write.csv(as.data.frame(kk), csv_out, row.names = FALSE)
  } else {
    kegg_results[[nm]] <- NULL
    cat("0 pathways\n")
  }
}

cat("\n")

# ============================================================
# Helper: GMT file parser (standalone, no internet required)
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
# 5. MSigDB Hallmark + C2 (Curated) enrichment
# ============================================================
cat("=== [5] MSigDB Hallmark + C2 Enrichment ===\n\n")

# --- 5a. Hallmark (from local GMT, no msigdbr dependency) ---
cat("--- [5a] Hallmark gene sets (GMT) ---\n\n")

GMT_HALLMARK <- "/tmp/h.all.v2024.1.symbols.gmt"
h_gmt <- tryCatch({
  read_gmt_simple(GMT_HALLMARK)
}, error = function(e) {
  cat(sprintf("    ERROR loading Hallmark GMT: %s\n", conditionMessage(e)))
  return(NULL)
})

if (is.null(h_gmt) || nrow(h_gmt) == 0) {
  cat("  SKIP Hallmark: failed to load GMT file\n\n")
  m_tier_results <- list()
} else {
  # Convert gene symbols to Entrez IDs for enricher()
  h_symbols <- unique(h_gmt$gene)
  h_map <- tryCatch({
    bitr(h_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    cat(sprintf("    WARNING bitr mapping failed: %s\n", conditionMessage(e)))
    data.frame(SYMBOL = character(0), ENTREZID = character(0), stringsAsFactors = FALSE)
  })
  # Merge and reformat for enricher(TERM2GENE = term + gene)
  h_gmt_merged <- merge(h_gmt, h_map, by.x = "gene", by.y = "SYMBOL", all.y = FALSE)
  h_gmt_entrez <- data.frame(
    term = h_gmt_merged$term,
    gene = h_gmt_merged$ENTREZID,
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("  GMT loaded: %d gene sets, %d symbols -> %d Entrez IDs\n",
              length(unique(h_gmt$term)), length(h_symbols),
              length(unique(h_gmt_entrez$gene))))
  cat(sprintf("  Unmapped symbols: %d\n", length(h_symbols) - nrow(h_map)))

  m_tier_results <- list()

  for (nm in names(entrez_lists)) {
    genes <- entrez_lists[[nm]]
    if (length(genes) < 5) {
      cat(sprintf("  %-8s: SKIP (< 5 Entrez IDs)\n", nm))
      m_tier_results[[nm]] <- NULL
      next
    }
    
    em <- tryCatch({
      enricher(gene          = genes,
               universe       = NULL,
               pvalueCutoff   = 0.05,
               pAdjustMethod  = "BH",
               minGSSize      = 5,
               TERM2GENE      = h_gmt_entrez)
    }, error = function(e) {
      cat(sprintf("    ERROR: %s\n", conditionMessage(e)))
      return(NULL)
    })
    
    if (!is.null(em) && nrow(em) > 0) {
      m_tier_results[[nm]] <- em
      cat(sprintf("  %-8s Hallmark: %d terms\n", nm, nrow(em)))
      
      csv_out <- file.path(ENRICH_DIR, sprintf("Hallmark_%s.csv", nm))
      write.csv(as.data.frame(em), csv_out, row.names = FALSE)
    } else {
      m_tier_results[[nm]] <- NULL
      cat(sprintf("  %-8s Hallmark: 0 terms\n", nm))
    }
  }
}  # end else (Hallmark)

# --- 5b. C2 Curated (from local GMT, CP:KEGG legacy) ---
cat("--- [5b] C2 Curated (GMT: c2.cp.kegg_legacy) ---\n\n")

GMT_C2 <- "/tmp/c2.cp.kegg_legacy.v2024.1.symbols.gmt"
c2_gmt <- tryCatch({
  read_gmt_simple(GMT_C2)
}, error = function(e) {
  cat(sprintf("  ERROR loading C2 GMT: %s\n", conditionMessage(e)))
  return(NULL)
})

if (is.null(c2_gmt) || nrow(c2_gmt) == 0) {
  cat("  SKIP C2: failed to load GMT file\n\n")
} else {
  # Convert symbols to Entrez
  c2_symbols <- unique(c2_gmt$gene)
  c2_map <- tryCatch({
    bitr(c2_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    cat(sprintf("    WARNING bitr mapping failed: %s\n", conditionMessage(e)))
    data.frame(SYMBOL = character(0), ENTREZID = character(0), stringsAsFactors = FALSE)
  })
  c2_merged <- merge(c2_gmt, c2_map, by.x = "gene", by.y = "SYMBOL", all.y = FALSE)
  c2_entrez <- data.frame(
    term = c2_merged$term,
    gene = c2_merged$ENTREZID,
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("  GMT loaded: %d gene sets, %d symbols -> %d Entrez IDs\n",
              length(unique(c2_gmt$term)), length(c2_symbols),
              length(unique(c2_entrez$gene))))
  cat(sprintf("  Unmapped symbols: %d\n", length(c2_symbols) - nrow(c2_map)))
  
  # Run C2 enrichment for All genes only
  cat("  Running C2 enricher (All genes)... ")
  em_c2 <- tryCatch({
    enricher(gene          = entrez_lists[["All"]],
             universe       = NULL,
             pvalueCutoff   = 0.05,
             pAdjustMethod  = "BH",
             minGSSize      = 5,
             TERM2GENE      = c2_entrez)
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", conditionMessage(e)))
    return(NULL)
  })
  
  if (!is.null(em_c2) && nrow(em_c2) > 0) {
    cat(sprintf("%d terms\n", nrow(em_c2)))
    m_tier_results[["C2_All"]] <- em_c2
    
    csv_out <- file.path(ENRICH_DIR, "C2_All.csv")
    write.csv(as.data.frame(em_c2), csv_out, row.names = FALSE)
  } else {
    cat("0 terms\n")
  }
}  # end else (C2)

cat("\n")

# ============================================================
# 6. Generate Plots
# ============================================================
cat("=== [6] Generating plots ===\n\n")

# Helper: plot GO/KEGG results
save_enrich_plots <- function(ego, prefix, suffix, n_show = 20) {
  if (is.null(ego) || nrow(ego) == 0) return(invisible())
  
  n_show <- min(nrow(ego), n_show)
  # CRITICAL: direct S4 slot assignment, NOT head() or [1:n_show,]
  #   DOSE/enrichit overwrites [.enrichResult and head.enrichResult S3
  #   methods, causing them to silently return a plain data.frame and
  #   breaking all enrichplot functions. Direct slot assignment preserves
  #   the S4 enrichResult class structure.
  ego_plot <- ego
  ego_plot@result <- ego_plot@result[1:n_show, , drop = FALSE]
  plot_title <- paste(prefix, suffix, sep = " ")
  
  # Barplot
  cat(sprintf("    %s %s: barplot ... ", prefix, suffix))
  p <- tryCatch(barplot(ego_plot, showCategory = n_show, title = plot_title),
                error = function(e) { cat(sprintf("FAIL (%s)\n", conditionMessage(e))); NULL })
  if (!is.null(p)) {
    f <- file.path(FIG_DIR, sprintf("%s_%s_barplot.pdf", prefix, suffix))
    pdf(f, width = 11, height = max(6, n_show * 0.35))
    print(p)
    dev.off()
    cat(sprintf("OK -> %s\n", basename(f)))
  }
  
  # Dotplot
  cat(sprintf("    %s %s: dotplot ... ", prefix, suffix))
  p <- tryCatch(dotplot(ego_plot, showCategory = n_show, title = plot_title),
                error = function(e) { cat(sprintf("FAIL (%s)\n", conditionMessage(e))); NULL })
  if (!is.null(p)) {
    f <- file.path(FIG_DIR, sprintf("%s_%s_dotplot.pdf", prefix, suffix))
    pdf(f, width = 11, height = max(6, n_show * 0.35))
    print(p)
    dev.off()
    cat(sprintf("OK -> %s\n", basename(f)))
  }
  
  # Emapplot (network) - only if >= 3 terms
  if (n_show >= 3) {
    cat(sprintf("    %s %s: emapplot ... ", prefix, suffix))
    p <- tryCatch({
      emapplot(pairwise_termsim(ego_plot), showCategory = min(n_show, 15))
    }, error = function(e) { cat(sprintf("FAIL (%s)\n", conditionMessage(e))); NULL })
    if (!is.null(p)) {
      f <- file.path(FIG_DIR, sprintf("%s_%s_emapplot.pdf", prefix, suffix))
      pdf(f, width = 12, height = 10)
      print(p)
      dev.off()
      cat(sprintf("OK -> %s\n", basename(f)))
    }
  }
  
  # Cnetplot (gene-concept) - only if <= 20 terms
  if (n_show <= 20) {
    cat(sprintf("    %s %s: cnetplot ... ", prefix, suffix))
    p <- tryCatch(cnetplot(ego_plot, showCategory = min(n_show, 10),
                           node_label = "category"),
                  error = function(e) { cat(sprintf("FAIL (%s)\n", conditionMessage(e))); NULL })
    if (!is.null(p)) {
      f <- file.path(FIG_DIR, sprintf("%s_%s_cnetplot.pdf", prefix, suffix))
      pdf(f, width = 16, height = 14)
      print(p)
      dev.off()
      cat(sprintf("OK -> %s\n", basename(f)))
    }
  }
  
  # Heatplot
  if (n_show >= 3) {
    cat(sprintf("    %s %s: heatplot ... ", prefix, suffix))
    p <- tryCatch(heatplot(ego_plot, showCategory = n_show),
                  error = function(e) { cat(sprintf("FAIL (%s)\n", conditionMessage(e))); NULL })
    if (!is.null(p)) {
      f <- file.path(FIG_DIR, sprintf("%s_%s_heatplot.pdf", prefix, suffix))
      pdf(f, width = 14, height = max(6, n_show * 0.4))
      print(p)
      dev.off()
      cat(sprintf("OK -> %s\n", basename(f)))
    }
  }
}

# GO plots
for (nm in names(go_results)) {
  ego <- go_results[[nm]]
  if (is.null(ego)) next
  
  parts <- strsplit(nm, "_")[[1]]
  ont <- parts[length(parts)]
  set_name <- paste(parts[-length(parts)], collapse = "_")
  save_enrich_plots(ego, set_name, ont)
}

# KEGG plots
for (nm in names(kegg_results)) {
  kk <- kegg_results[[nm]]
  if (is.null(kk)) next
  save_enrich_plots(kk, nm, "KEGG")
}

# MSigDB Hallmark plots
for (nm in names(m_tier_results)) {
  em <- m_tier_results[[nm]]
  if (is.null(em)) next
  
  # Only do dotplot + barplot for MSigDB (emapplot/cnetplot may not work
  # well with enricher output)
  
  if (grepl("^All_", nm) || nm %in% names(entrez_lists)) {
    n_show <- min(nrow(em), 20)
    # CRITICAL: direct S4 slot assignment (see save_enrich_plots for rationale)
    em_plot <- em
    em_plot@result <- em_plot@result[1:n_show, , drop = FALSE]
    
    p <- tryCatch(dotplot(em_plot, showCategory = n_show,
                          title = paste("MSigDB", nm)),
                  error = function(e) NULL)
    if (!is.null(p)) {
      f <- file.path(FIG_DIR, sprintf("MSigDB_%s_dotplot.pdf", nm))
      pdf(f, width = 14, height = max(8, n_show * 0.35))
      print(p)
      dev.off()
    }
    
    p <- tryCatch(barplot(em_plot, showCategory = n_show,
                          title = paste("MSigDB", nm)),
                  error = function(e) NULL)
    if (!is.null(p)) {
      f <- file.path(FIG_DIR, sprintf("MSigDB_%s_barplot.pdf", nm))
      pdf(f, width = 14, height = max(8, n_show * 0.35))
      print(p)
      dev.off()
    }
  }
}

# ---- Comparison dotplot: GO-BP across All / Up / Down ----
cat("  Creating comparison dotplot (GO-BP: All vs Up vs Down)...\n")

comp_go <- list()
for (nm in c("All", "Up", "Down")) {
  key <- paste(nm, "BP", sep = "_")
  if (!is.null(go_results[[key]]) && nrow(go_results[[key]]) > 0) {
    comp_go[[nm]] <- go_results[[key]]
  }
}

if (length(comp_go) >= 2) {
  tryCatch({
    p <- dotplot(comp_go, showCategory = 10,
                 title = "GO-BP: All vs Up vs Down (253 Meta Genes)",
                 label_format = 60) +
      theme(legend.position = "right")
    f <- file.path(FIG_DIR, "GO_BP_All_Up_Down_comparison.pdf")
    pdf(f, width = 14, height = 10)
    print(p)
    dev.off()
    cat("  Comparison dotplot saved\n")
  }, error = function(e) {
    cat(sprintf("  Comparison dotplot skipped: %s\n", conditionMessage(e)))
  })
}
cat("\n")

# ============================================================
# 7. Summary Report
# ============================================================
cat("=== [7] Generating summary report ===\n\n")

report <- c(
  "========================================================================",
  "Pathway Enrichment on 253 Meta-Analysis Significant Genes",
  paste("Date:", Sys.Date()),
  "========================================================================",
  "",
  sprintf("Total significant genes: %d", nrow(sig)),
  sprintf("Up-regulated (consensus): %d", length(genes_up)),
  sprintf("Down-regulated (consensus): %d", length(genes_down)),
  "",
  "--- GO ORA ---"
)

for (nm in names(go_results)) {
  n <- if (!is.null(go_results[[nm]])) nrow(go_results[[nm]]) else 0
  report <- c(report, sprintf("  %-16s: %d terms", nm, n))
}

report <- c(report, "",
            "--- KEGG ORA ---")
for (nm in names(kegg_results)) {
  n <- if (!is.null(kegg_results[[nm]])) nrow(kegg_results[[nm]]) else 0
  report <- c(report, sprintf("  %-16s: %d pathways", nm, n))
}

report <- c(report, "",
            "--- MSigDB Hallmark ---")
for (nm in names(m_tier_results)) {
  if (grepl("^C2_", nm)) next
  n <- if (!is.null(m_tier_results[[nm]])) nrow(m_tier_results[[nm]]) else 0
  report <- c(report, sprintf("  %-16s: %d terms", nm, n))
}

report <- c(report, "",
            "--- MSigDB C2 (All genes only) ---")
if (!is.null(m_tier_results[["C2_All"]])) {
  report <- c(report, sprintf("  C2_All: %d terms", nrow(m_tier_results[["C2_All"]])))
}

report <- c(report, "",
            "--- Top 10 GO-BP (All genes) ---")
if (!is.null(go_results[["All_BP"]]) && nrow(go_results[["All_BP"]]) > 0) {
  top10 <- head(go_results[["All_BP"]], 10)
  for (i in seq_len(nrow(top10))) {
    report <- c(report, sprintf("  %s | p.adj=%.2e | Count=%d",
                                top10$Description[i], top10$p.adjust[i], top10$Count[i]))
  }
}

report <- c(report, "",
            "--- Top 10 KEGG (All genes) ---")
if (!is.null(kegg_results[["All"]]) && nrow(kegg_results[["All"]]) > 0) {
  top10 <- head(kegg_results[["All"]], 10)
  for (i in seq_len(nrow(top10))) {
    report <- c(report, sprintf("  %s | p.adj=%.2e | Count=%d",
                                top10$Description[i], top10$p.adjust[i], top10$Count[i]))
  }
}

report <- c(report, "",
            "--- Top 10 Hallmark (All genes) ---")
if (!is.null(m_tier_results[["All"]])) {
  top10 <- head(m_tier_results[["All"]], 10)
  for (i in seq_len(nrow(top10))) {
    report <- c(report, sprintf("  %s | p.adj=%.2e | Count=%d",
                                top10$Description[i], top10$p.adjust[i], top10$Count[i]))
  }
}

report <- c(report, "",
            "========================================================================")

writeLines(report, file.path(ENRICH_DIR, "enrichment_summary_253.txt"))
cat(paste(report, collapse = "\n"), "\n")

# ============================================================
# 8. Done
# ============================================================
cat("\n============================================================\n")
cat("ENRICHMENT ANALYSIS v17 COMPLETE\n")
cat("============================================================\n")
cat(sprintf("  Tables:  %s\n", ENRICH_DIR))
cat(sprintf("  Figures: %s\n", FIG_DIR))
cat(sprintf("  Summary: %s\n",
            file.path(ENRICH_DIR, "enrichment_summary_253.txt")))
cat("\n")
