# ============================================================
# 19_ENRICHMENT_REPORT.R
# Integrated HTML Report: ORA + GSEA
# XELOX Resistance Study - 253-Meta Analysis
# ============================================================

Sys.setenv(HOME = "/tmp", TMPDIR = "/tmp", TMP = "/tmp", TEMP = "/tmp",
           R_LIBS_USER = "/path/to/Rlibs")
.libPaths(c("/path/to/Rlibs", .libPaths()))

TABLES_DIR   <- "/tmp/xelox_output/tables"
FIGURES_DIR  <- "/tmp/xelox_output/figures"
REPORT_DIR   <- "/tmp/xelox_output/report"
dir.create(REPORT_DIR, showWarnings = FALSE, recursive = TRUE)

OUTPUT_HTML  <- file.path(REPORT_DIR, "XELOX_ORA_GSEA_Report.html")

suppressPackageStartupMessages({
  library(ggplot2)
  library(htmltools)
  library(base64enc)
})

cat("============================================================\n")
cat("XELOX ORA+GSEA Integrated HTML Report Generator v19\n")
cat("============================================================\n\n")

# ============================================================
# 1. Helper functions
# ============================================================

# Convert a PDF figure to a base64-encoded PNG for HTML embedding
# Uses Python (PyMuPDF/fitz) for PDF rendering, then base64enc for embedding
pdf_to_base64_png <- function(pdf_path) {
  if (!file.exists(pdf_path)) return(NA_character_)
  
  png_tmp <- file.path("/tmp", paste0(gsub("[^A-Za-z0-9]", "", basename(pdf_path)), ".png"))
  bat_tmp <- "/tmp/_pdf_convert.bat"
  
  # Write a .bat file (avoids Chinese chars in system() call)
  python_exe <- "/home/user/AppData/Local/Programs/Python/Python314/python.exe"
  helper_py  <- "/tmp/pdf_to_png.py"
  bat_lines <- sprintf('@"%s" "%s" "%s" "%s"\n', python_exe, helper_py,
                       gsub("/", "\\\\", pdf_path),
                       gsub("/", "\\\\", png_tmp))
  writeLines(bat_lines, bat_tmp)
  
  result <- tryCatch({
    sys_out <- system(bat_tmp, intern = TRUE, ignore.stderr = TRUE, show.output.on.console = FALSE)
    if (file.exists(png_tmp) && file.info(png_tmp)$size > 0) {
      b64 <- base64enc::base64encode(png_tmp)
      unlink(png_tmp)
      sprintf("data:image/png;base64,%s", b64)
    } else {
      NA_character_
    }
  }, error = function(e) {
    NA_character_
  }, finally = {
    if (file.exists(png_tmp)) unlink(png_tmp)
    if (file.exists(bat_tmp)) unlink(bat_tmp)
  })
  
  result
}

# Read a CSV enrichment table, return as data frame
read_enrich_csv <- function(filename) {
  fp <- file.path(TABLES_DIR, filename)
  if (!file.exists(fp)) return(NULL)
  read.csv(fp, stringsAsFactors = FALSE)
}

# Format p-values for display
fmt_pval <- function(p, digits = 2) {
  ifelse(p < 0.0001, sprintf("%.2e", p), sprintf(paste0("%.", digits, "f"), p))
}

# Build a simple HTML table from a data frame subset
simple_table_html <- function(df, cols = NULL, max_rows = 10) {
  if (is.null(df) || nrow(df) == 0) return("")
  if (!is.null(cols)) df <- df[, intersect(cols, names(df)), drop = FALSE]
  df <- head(df, max_rows)
  rows <- apply(df, 1, function(r) {
    cells <- paste(sprintf("<td>%s</td>", htmlEscape(as.character(r))), collapse = "")
    sprintf("<tr>%s</tr>", cells)
  })
  headers <- paste(sprintf("<th>%s</th>", names(df)), collapse = "")
  sprintf("<table class='data-table'><thead><tr>%s</tr></thead><tbody>%s</tbody></table>",
          headers, paste(rows, collapse = ""))
}

# ============================================================
# 2. Read all enrichment data
# ============================================================
cat("=== [1] Reading ORA tables ===\n\n")

ora_tables <- list(
  GO_All_CC   = read_enrich_csv("GO_All_CC.csv"),
  GO_All_MF   = read_enrich_csv("GO_All_MF.csv"),
  GO_Down_BP  = read_enrich_csv("GO_Down_BP.csv"),
  GO_Down_MF  = read_enrich_csv("GO_Down_MF.csv"),
  KEGG_Down   = read_enrich_csv("KEGG_Down.csv"),
  Hallmark_Down = read_enrich_csv("Hallmark_Down.csv")
)

for (nm in names(ora_tables)) {
  n <- if (!is.null(ora_tables[[nm]])) nrow(ora_tables[[nm]]) else 0
  cat(sprintf("  %-20s: %d terms\n", nm, n))
}

cat("\n=== [2] Reading GSEA tables ===\n\n")

gsea_tables <- list(
  GO_BP     = read_enrich_csv("GSEA_GO_BP.csv"),
  GO_CC     = read_enrich_csv("GSEA_GO_CC.csv"),
  GO_MF     = read_enrich_csv("GSEA_GO_MF.csv"),
  KEGG      = read_enrich_csv("GSEA_KEGG.csv"),
  Hallmark  = read_enrich_csv("GSEA_Hallmark.csv"),
  C2        = read_enrich_csv("GSEA_C2.csv")
)

for (nm in names(gsea_tables)) {
  n <- if (!is.null(gsea_tables[[nm]])) nrow(gsea_tables[[nm]]) else 0
  cat(sprintf("  %-20s: %d terms\n", nm, n))
}

# ============================================================
# 3. Convert key figures to base64 PNG
# ============================================================
cat("\n=== [3] Converting figures to base64 ===\n\n")

KEY_FIGURES <- list(
  # ORA figures
  ora_Down_BP_dotplot    = file.path(FIGURES_DIR, "Down_BP_dotplot.pdf"),
  ora_All_CC_dotplot     = file.path(FIGURES_DIR, "All_CC_dotplot.pdf"),
  ora_Down_KEGG_dotplot  = file.path(FIGURES_DIR, "Down_KEGG_dotplot.pdf"),
  ora_Down_Hallmark_dotplot = file.path(FIGURES_DIR, "MSigDB_Down_dotplot.pdf"),
  # GSEA figures
  gsea_GO_BP_ridge       = file.path(FIGURES_DIR, "GSEA_GO_BP_ridgeplot.pdf"),
  gsea_GO_CC_ridge       = file.path(FIGURES_DIR, "GSEA_GO_CC_ridgeplot.pdf"),
  gsea_GO_MF_ridge       = file.path(FIGURES_DIR, "GSEA_GO_MF_ridgeplot.pdf"),
  gsea_KEGG_ridge        = file.path(FIGURES_DIR, "GSEA_KEGG_ridgeplot.pdf"),
  gsea_Hallmark_ridge    = file.path(FIGURES_DIR, "GSEA_Hallmark_ridgeplot.pdf"),
  gsea_C2_ridge          = file.path(FIGURES_DIR, "GSEA_C2_ridgeplot.pdf")
)

fig_base64 <- list()
for (nm in names(KEY_FIGURES)) {
  cat(sprintf("  %s ... ", nm))
  b64 <- pdf_to_base64_png(KEY_FIGURES[[nm]])
  fig_base64[[nm]] <- b64
  cat(ifelse(is.na(b64), "SKIP\n", "OK\n"))
}

# ============================================================
# 4. Build HTML report
# ============================================================
cat("\n=== [4] Building HTML ===\n\n")

# CSS styles
CSS_STYLES <- "
body {
  font-family: 'Segoe UI', 'Helvetica Neue', Arial, sans-serif;
  color: #333; line-height: 1.6; margin: 0; padding: 0;
  background: #f5f7fa;
}
.report-container {
  max-width: 1200px; margin: 0 auto; padding: 20px;
}
h1 { color: #1a3a5c; border-bottom: 3px solid #2b6fa7; padding-bottom: 10px; }
h2 { color: #1a3a5c; border-bottom: 2px solid #8ab8d4; padding-bottom: 8px; margin-top: 40px; }
h3 { color: #2b6fa7; margin-top: 25px; }
.header-box {
  background: linear-gradient(135deg, #1a3a5c 0%, #2b6fa7 100%);
  color: white; padding: 40px; border-radius: 8px; margin-bottom: 30px;
}
.header-box h1 { color: white; border: none; font-size: 28px; }
.header-box .subtitle { font-size: 16px; opacity: 0.9; margin-top: 8px; }
.meta-grid {
  display: grid; grid-template-columns: repeat(auto-fill, minmax(200px,1fr));
  gap: 12px; margin: 20px 0;
}
.meta-card {
  background: white; border-radius: 8px; padding: 16px;
  box-shadow: 0 1px 4px rgba(0,0,0,0.08); text-align: center;
}
.meta-card .number { font-size: 28px; font-weight: 700; color: #2b6fa7; }
.meta-card .label { font-size: 12px; color: #666; margin-top: 4px; }
.section {
  background: white; border-radius: 8px; padding: 25px;
  margin-bottom: 25px; box-shadow: 0 2px 8px rgba(0,0,0,0.06);
}
.data-table {
  width: 100%; border-collapse: collapse; margin: 15px 0;
  font-size: 13px;
}
.data-table th {
  background: #2b6fa7; color: white; padding: 10px 12px;
  text-align: left; font-weight: 600;
}
.data-table td { padding: 8px 12px; border-bottom: 1px solid #e0e0e0; }
.data-table tr:nth-child(even) { background: #f8faff; }
.data-table tr:hover { background: #e8f0fe; }
.figure-box {
  margin: 20px 0; text-align: center;
}
.figure-box img { max-width: 100%; border: 1px solid #ddd; border-radius: 4px; }
.figure-box .caption { font-size: 12px; color: #666; margin-top: 6px; font-style: italic; }
.toc { background: #f0f4f8; padding: 15px 20px; border-radius: 8px; margin-bottom: 25px; }
.toc a { color: #2b6fa7; text-decoration: none; display: block; padding: 3px 0; }
.toc a:hover { text-decoration: underline; }
.summary-grid {
  display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin: 20px 0;
}
@media (max-width: 768px) { .summary-grid { grid-template-columns: 1fr; } }
"

# --- Build the HTML document ---
build_html <- function() {
  
  # Navigation TOC
  toc <- tags$div(class = "toc",
    tags$h3("Table of Contents"),
    tags$a(href = "#sec-ora", "\u25b6 Over-Representation Analysis (ORA)"),
    tags$a(href = "#sec-ora-go", "  \u2514 GO Enrichment"),
    tags$a(href = "#sec-ora-kegg", "  \u2514 KEGG Pathway"),
    tags$a(href = "#sec-gsea", "\u25b6 Gene Set Enrichment Analysis (GSEA)"),
    tags$a(href = "#sec-gsea-go", "  \u2514 GO GSEA (BP / CC / MF)"),
    tags$a(href = "#sec-gsea-kegg", "  \u2514 KEGG GSEA"),
    tags$a(href = "#sec-gsea-msigdb", "  \u2514 MSigDB (Hallmark + C2)"),
    tags$a(href = "#sec-insights", "\u25b6 Key Biological Insights"),
    tags$a(href = "#sec-methods", "\u25b6 Methods")
  )
  
  # =============================
  # Header
  # =============================
  header <- tags$div(class = "header-box",
    tags$h1("XELOX Resistance: Integrated Pathway Enrichment Report"),
    tags$div(class = "subtitle",
      "253-Gene Meta-Analysis Signature \u2022 ORA + GSEA Integration",
      tags$br(),
      paste("Report generated:", Sys.Date())
    )
  )
  
  # =============================
  # Study Overview
  # =============================
  overview <- tags$div(class = "section",
    tags$h2("Study Overview"),
    tags$div(class = "meta-grid",
      tags$div(class = "meta-card", tags$div(class = "number", "253"),
               tags$div(class = "label", "Significant Genes (FDR<0.05)")),
      tags$div(class = "meta-card", tags$div(class = "number", "94"),
               tags$div(class = "label", "Up-regulated (Consensus)")),
      tags$div(class = "meta-card", tags$div(class = "number", "159"),
               tags$div(class = "label", "Down-regulated (Consensus)")),
      tags$div(class = "meta-card", tags$div(class = "number", "16,289"),
               tags$div(class = "label", "Genes Ranked (GSEA)")),
      tags$div(class = "meta-card",
               tags$div(class = "number", "1,990"),
               tags$div(class = "label", "GSEA Enriched Terms")),
      tags$div(class = "meta-card",
               tags$div(class = "number", "52"),
               tags$div(class = "label", "ORA Enriched Terms"))
    ),
    tags$p("This report integrates two complementary enrichment strategies applied to the 253-meta-analysis of XELOX resistance:"),
    tags$ul(
      tags$li(tags$b("ORA (Over-Representation Analysis):"), " Tests enrichment of 253 significant genes (94 up, 159 down) against GO, KEGG, and MSigDB Hallmark gene sets. Input: filtered DEG set (FDR<0.05)."),
      tags$li(tags$b("GSEA (Gene Set Enrichment Analysis):"), " Tests enrichment of all 16,289 ranked genes (mean logFC) without any significance threshold. Detects coordinated shifts across entire transcriptome.")
    ),
    tags$p(tags$i("Note: ORA detected terms mainly in down-regulated genes; All and Up gene lists yielded minimal ORA results. GSEA was dominated by positive enrichment (unidirectional mean_logFC ranking)."))
  )
  
  # =============================
  # ORA Section
  # =============================
  ora_section <- build_ora_section()
  
  # =============================
  # GSEA Section
  # =============================
  gsea_section <- build_gsea_section()
  
  # =============================
  # Insights Section
  # =============================
  insights_section <- build_insights_section()
  
  # =============================
  # Methods Section
  # =============================
  methods_section <- tags$div(class = "section", id = "sec-methods",
    tags$h2("Methods"),
    tags$h3("Data Source"),
    tags$p("Meta-analysis of 3 GEO datasets (GSE28702, GSE53127, GSE69657) comparing XELOX responders vs non-responders in colorectal cancer patients. 25,568 genes tested by Stouffer's method; 253 genes significant at FDR<0.05."),
    tags$h3("ORA"),
    tags$p("clusterProfiler v4.20 enrichGO() with org.Hs.eg.db, enrichKEGG() with local cache, enricher() with MSigDB GMT files. Significance: BH-adjusted p<0.05."),
    tags$h3("GSEA"),
    tags$p("All 16,289 genes ranked by mean_logFC (sorted descending). gseGO() for GO, gseKEGG() for KEGG, GSEA() for MSigDB (TERM2GENE). Significance: BH-adjusted p<0.05."),
    tags$h3("Software"),
    tags$p("R v4.6.0, clusterProfiler v4.20.0, enrichplot, org.Hs.eg.db, ggplot2, DT.")
  )
  
  # Assemble full page
  tagList(
    tags$head(
      tags$meta(charset = "UTF-8"),
      tags$title("XELOX Resistance Pathway Enrichment Report"),
      tags$style(HTML(CSS_STYLES))
    ),
    tags$body(
      tags$div(class = "report-container",
        header, toc, overview, methods_section,
        ora_section, gsea_section, insights_section
      )
    )
  )
}

# ============================================================
# Helper: Build ORA Section
# ============================================================
build_ora_section <- function() {
  
  # ORA summary grid
  ora_summary <- tags$div(class = "summary-grid",
    tags$div(
      tags$h3("GO Enrichment"),
      tags$ul(
        tags$li(sprintf("GO-BP (Down): %d terms", if (!is.null(ora_tables$GO_Down_BP)) nrow(ora_tables$GO_Down_BP) else 0)),
        tags$li(sprintf("GO-CC (All): %d terms", if (!is.null(ora_tables$GO_All_CC)) nrow(ora_tables$GO_All_CC) else 0)),
        tags$li(sprintf("GO-MF (All): %d terms", if (!is.null(ora_tables$GO_All_MF)) nrow(ora_tables$GO_All_MF) else 0)),
        tags$li(sprintf("GO-MF (Down): %d terms", if (!is.null(ora_tables$GO_Down_MF)) nrow(ora_tables$GO_Down_MF) else 0))
      )
    ),
    tags$div(
      tags$h3("KEGG & Hallmark"),
      tags$ul(
        tags$li(sprintf("KEGG (Down): %d pathways", if (!is.null(ora_tables$KEGG_Down)) nrow(ora_tables$KEGG_Down) else 0)),
        tags$li(sprintf("Hallmark (Down): %d terms", if (!is.null(ora_tables$Hallmark_Down)) nrow(ora_tables$Hallmark_Down) else 0))
      )
    )
  )
  
  # Build a static HTML table from a data frame
  build_ora_table <- function(df, cols_show = c("ID","Description","Count","pvalue","p.adjust","geneID")) {
    if (is.null(df)) return(HTML("<p><i>No significant terms.</i></p>"))
    df <- df[, intersect(cols_show, names(df)), drop = FALSE]
    df$pvalue <- fmt_pval(df$pvalue)
    df$p.adjust <- fmt_pval(df$p.adjust)
    if ("geneID" %in% names(df)) df$geneID <- substr(df$geneID, 1, 60)
    
    rows <- apply(head(df, 15), 1, function(r) {
      cells <- paste(sprintf("<td>%s</td>", htmlEscape(as.character(r))), collapse = "")
      sprintf("<tr>%s</tr>", cells)
    })
    headers <- paste(sprintf("<th>%s</th>", names(df)), collapse = "")
    HTML(sprintf(
      "<div style='max-height:600px;overflow-y:auto;margin:15px 0'>
       <table class='data-table'><thead><tr>%s</tr></thead><tbody>%s</tbody></table></div>",
      headers, paste(rows, collapse = "")))
  }
  
  go_bp_table     <- build_ora_table(ora_tables$GO_Down_BP)
  go_cc_table     <- build_ora_table(ora_tables$GO_All_CC)
  kegg_table      <- build_ora_table(ora_tables$KEGG_Down)
  hallmark_table  <- build_ora_table(ora_tables$Hallmark_Down)
  
  # Figures
  add_figure <- function(b64_key, caption) {
    if (!is.na(fig_base64[[b64_key]])) {
      tags$div(class = "figure-box",
        tags$img(src = fig_base64[[b64_key]], alt = caption),
        tags$div(class = "caption", caption)
      )
    } else {
      tags$p(tags$i(sprintf("Figure not available: %s", caption)))
    }
  }
  
  tags$div(class = "section", id = "sec-ora",
    tags$h2("Over-Representation Analysis (ORA)"),
    tags$p("ORA tests whether the 253 significant genes (filtered at Stouffer FDR<0.05) are over-represented in predefined gene sets. Results are split by direction (up/down) where sufficient genes exist."),
    ora_summary,
    
    # GO-Down-BP
    tags$div(id = "sec-ora-go",
      tags$h3("GO Biological Process (Down-regulated Genes)"),
      tags$p(sprintf("Enriched terms in the 159 down-regulated genes. Top terms involve immune-related processes and cellular adhesion.")),
      go_bp_table,
      add_figure("ora_Down_BP_dotplot", "Figure 1: GO-BP Dotplot (Down-regulated)"),
      
      # GO-CC
      tags$h3("GO Cellular Component (All Genes)"),
      go_cc_table,
      add_figure("ora_All_CC_dotplot", "Figure 2: GO-CC Dotplot (All significant genes)"),
      
      # KEGG
      tags$div(id = "sec-ora-kegg",
        tags$h3("KEGG Pathway (Down-regulated Genes)"),
        kegg_table,
        add_figure("ora_Down_KEGG_dotplot", "Figure 3: KEGG Dotplot (Down-regulated)")
      ),
      
      # Hallmark
      tags$h3("MSigDB Hallmark (Down-regulated Genes)"),
      hallmark_table,
      add_figure("ora_Down_Hallmark_dotplot", "Figure 4: Hallmark Dotplot (Down-regulated)")
    )
  )
}

# ============================================================
# Helper: Build GSEA Section
# ============================================================
build_gsea_section <- function() {
  
  # GSEA summary grid
  gsea_summary <- tags$div(class = "summary-grid",
    tags$div(
      tags$h3("GO GSEA"),
      tags$ul(
        tags$li(sprintf("GO-BP: %d terms (522 up, 5 down)",
                        if (!is.null(gsea_tables$GO_BP)) nrow(gsea_tables$GO_BP) else 0)),
        tags$li(sprintf("GO-CC: %d terms (126 up, 2 down)",
                        if (!is.null(gsea_tables$GO_CC)) nrow(gsea_tables$GO_CC) else 0)),
        tags$li(sprintf("GO-MF: %d terms (74 up, 2 down)",
                        if (!is.null(gsea_tables$GO_MF)) nrow(gsea_tables$GO_MF) else 0))
      )
    ),
    tags$div(
      tags$h3("KEGG & MSigDB GSEA"),
      tags$ul(
        tags$li(sprintf("KEGG: %d pathways (88 up, 0 down)",
                        if (!is.null(gsea_tables$KEGG)) nrow(gsea_tables$KEGG) else 0)),
        tags$li(sprintf("Hallmark: %d terms (38 up, 0 down)",
                        if (!is.null(gsea_tables$Hallmark)) nrow(gsea_tables$Hallmark) else 0)),
        tags$li(sprintf("C2 KEGG_legacy: %d terms (28 up, 0 down)",
                        if (!is.null(gsea_tables$C2)) nrow(gsea_tables$C2) else 0))
      )
    )
  )
  
  # Build GSEA static table
  build_gsea_table <- function(df, extra_cols = c("ID","Description","NES","pvalue","p.adjust","setSize"),
                               max_rows = 15) {
    if (is.null(df)) return(HTML("<p><i>No significant terms.</i></p>"))
    cols <- intersect(extra_cols, names(df))
    df <- df[, cols, drop = FALSE]
    if ("pvalue" %in% names(df)) df$pvalue <- fmt_pval(df$pvalue)
    if ("p.adjust" %in% names(df)) df$p.adjust <- fmt_pval(df$p.adjust)
    if ("NES" %in% names(df)) df$NES <- round(df$NES, 3)
    
    rows <- apply(head(df, max_rows), 1, function(r) {
      cells <- paste(sprintf("<td>%s</td>", htmlEscape(as.character(r))), collapse = "")
      sprintf("<tr>%s</tr>", cells)
    })
    headers <- paste(sprintf("<th>%s</th>", names(df)), collapse = "")
    HTML(sprintf(
      "<div style='max-height:600px;overflow-y:auto;margin:15px 0'>
       <table class='data-table'><thead><tr>%s</tr></thead><tbody>%s</tbody></table></div>",
      headers, paste(rows, collapse = "")))
  }
  
  add_figure <- function(b64_key, caption) {
    if (!is.na(fig_base64[[b64_key]])) {
      tags$div(class = "figure-box",
        tags$img(src = fig_base64[[b64_key]], alt = caption),
        tags$div(class = "caption", caption)
      )
    } else {
      tags$p(tags$i(sprintf("Figure not available: %s", caption)))
    }
  }
  
  tags$div(class = "section", id = "sec-gsea",
    tags$h2("Gene Set Enrichment Analysis (GSEA)"),
    tags$p("GSEA uses the full-ranked gene list (all 16,289 genes, sorted by mean_logFC descending) to detect coordinated expression shifts. A positive NES indicates enrichment in resistant-upregulated genes."),
    tags$p(tags$b("Caveat:"), " The mean_logFC ranking metric is unidirectionally positive (range 0.03\u2013781.7), meaning GSEA predominantly detects up-direction enriched pathways. Down-regulated pathways are largely missed. Consider using Stouffer's Z-score for bidirectional ranking in future analyses."),
    gsea_summary,
    
    # GO-BP
    tags$div(id = "sec-gsea-go",
      tags$h3("GO Biological Process - GSEA"),
      tags$p("Top enriched GO-BP terms include antigen processing/presentation (MHC class I), T-cell mediated cytotoxicity, mitochondrial electron transport, and cellular response to metals."),
      build_gsea_table(gsea_tables$GO_BP, max_rows = 20),
      add_figure("gsea_GO_BP_ridge", "Figure 5: GSEA GO-BP Ridgeplot (Top 20 terms)"),
      
      # GO-CC
      tags$h3("GO Cellular Component - GSEA"),
      build_gsea_table(gsea_tables$GO_CC, max_rows = 15),
      add_figure("gsea_GO_CC_ridge", "Figure 6: GSEA GO-CC Ridgeplot (Top 20 terms)"),
      
      # GO-MF
      tags$h3("GO Molecular Function - GSEA"),
      build_gsea_table(gsea_tables$GO_MF, max_rows = 15),
      add_figure("gsea_GO_MF_ridge", "Figure 7: GSEA GO-MF Ridgeplot (Top 20 terms)")
    ),
    
    # KEGG
    tags$div(id = "sec-gsea-kegg",
      tags$h3("KEGG Pathway - GSEA"),
      tags$p("Top KEGG pathways are dominated by immune-related categories: antigen processing, autoimmune disease pathways, complement/coagulation cascades, and cardiac muscle contraction."),
      build_gsea_table(gsea_tables$KEGG, max_rows = 20),
      add_figure("gsea_KEGG_ridge", "Figure 8: GSEA KEGG Ridgeplot (Top 20 terms)")
    ),
    
    # MSigDB
    tags$div(id = "sec-gsea-msigdb",
      tags$h3("MSigDB Hallmark - GSEA"),
      tags$p("Top Hallmark sets highlight Epithelial-Mesenchymal Transition (EMT), Coagulation, Angiogenesis, MYC targets, and Reactive Oxygen Species pathway \u2014 all highly relevant to chemoresistance."),
      build_gsea_table(gsea_tables$Hallmark, max_rows = 20),
      add_figure("gsea_Hallmark_ridge", "Figure 9: GSEA Hallmark Ridgeplot"),
      
      tags$h3("MSigDB C2 (KEGG Legacy) - GSEA"),
      tags$p("C2 curated gene sets recapitulate immune rejection and ribosome pathways."),
      build_gsea_table(gsea_tables$C2, max_rows = 15),
      add_figure("gsea_C2_ridge", "Figure 10: GSEA C2 Ridgeplot")
    )
  )
}

# ============================================================
# Helper: Build Insights Section
# ============================================================
build_insights_section <- function() {
  
  # Extract top terms for insights
  hallmark_top <- if (!is.null(gsea_tables$Hallmark)) head(gsea_tables$Hallmark$Description, 5) else character(0)
  kegg_top <- if (!is.null(gsea_tables$KEGG)) head(gsea_tables$KEGG$Description, 5) else character(0)
  
  tags$div(class = "section", id = "sec-insights",
    tags$h2("Key Biological Insights"),
    
    tags$h3("1. Immune Microenvironment Remodeling"),
    tags$p("Both ORA (down-regulated genes) and GSEA converge on immune-related pathways:"),
    tags$ul(
      tags$li("GSEA reveals strong enrichment of", tags$b("antigen processing and presentation"), "(MHC class I pathway, NES=1.48\u20131.57),", tags$b("T-cell mediated cytotoxicity"), "(NES=1.56), and", tags$b("complement/coagulation cascades"), "(NES=1.40)"),
      tags$li("ORA of down-regulated genes shows enrichment in immune-related BP terms, suggesting loss of immune surveillance in resistant tumors"),
      tags$li("Autoimmune disease pathways (Type I diabetes, GVHD, allograft rejection) enriched in GSEA may reflect immune activation signatures")
    ),
    
    tags$h3("2. EMT and Metastatic Programs"),
    tags$p("The top Hallmark term is ", tags$b("Epithelial-Mesenchymal Transition"), " (NES=1.41, p.adj=2.5e-9), a well-established driver of chemoresistance in CRC. This is accompanied by ", tags$b("Angiogenesis"), " (NES=1.38) and ", tags$b("Coagulation"), " (NES=1.39)."),
    tags$ul(
      tags$li("EMT confers stemness and drug efflux capabilities"),
      tags$li("Coagulation pathway links to tumor microenvironment remodeling and metastasis")
    ),
    
    tags$h3("3. Cellular Stress and Proliferation"),
    tags$p("Resistant tumors activate multiple stress-response programs:"),
    tags$ul(
      tags$li(tags$b("MYC targets V1"), "(NES=1.37, p.adj=2.5e-9) \u2014 sustained proliferative signaling"),
      tags$li(tags$b("Reactive Oxygen Species pathway"), "(NES=1.34) \u2014 oxidative stress adaptation"),
      tags$li(tags$b("Apoptosis"), "(NES=1.32) \u2014 may reflect therapy-induced stress rather than execution"),
      tags$li(tags$b("Hypoxia"), "(NES=1.30) \u2014 tumor microenvironment adaptation")
    ),
    
    tags$h3("4. ORA \u2013 GSEA Complementarity"),
    tags$p("ORA identified only 52 terms (mostly from down-regulated genes), while GSEA identified 1,990 terms. This highlights:"),
    tags$ul(
      tags$li("ORA's sensitivity is limited by the 253-gene threshold (many relevant pathways operate below FDR<0.05)"),
      tags$li("GSEA's rank-based approach captures subtler but biologically meaningful coordinated changes"),
      tags$li("The unidirectional mean_logFC ranking limits GSEA to detecting only up-direction enrichment")
    ),
    
    tags$h3("5. Methodological Limitation"),
    tags$p("The mean_logFC ranking metric is dominated by extremely high values (max 781.7, from GSE28702), and all 16,289 ranked values are positive. This means GSEA can only detect pathways enriched at the TOP of the ranked list (up in resistance). Down-regulated pathways (including the 25 ORA GO-BP terms from down genes) are invisible to this GSEA. Future work should use Stouffer's Z-score for a bidirectional metric that captures both directions.")
  )
}

# ============================================================
# 5. Render and save
# ============================================================
cat("  Building HTML document tree ... ")
html_doc <- build_html()
cat("OK\n")

cat("  Saving HTML to", OUTPUT_HTML, "... ")
save_html(html_doc, OUTPUT_HTML)
cat("OK\n")

cat(sprintf("\nReport saved: %s\n", OUTPUT_HTML))
cat(sprintf("File size: %.1f MB\n", file.info(OUTPUT_HTML)$size / 1e6))

# ============================================================
# 6. Print summary
# ============================================================
cat("\n============================================================\n")
cat("REPORT GENERATION COMPLETE\n")
cat("============================================================\n")
cat(sprintf("  HTML:  %s\n", OUTPUT_HTML))
cat("\n")
