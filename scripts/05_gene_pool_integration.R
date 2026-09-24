# ============================================================
# 05_GENE_POOL_INTEGRATION.R (v2 - Gene Symbol Aware)
# Capecitabine metabolism gene enrichment & final gene pool
# Uses mapped gene symbols (Step 6 output)
# XELOX Resistance Study
# ============================================================

Sys.setenv(TMPDIR = "/tmp", TMP = "/tmp", TEMP = "/tmp")
.libPaths(c("/path/to/Rlibs", .libPaths()))

PROJECT_ROOT <- "/path/to/xelox_project"
DATA_GEO_DIR  <- file.path(PROJECT_ROOT, "data", "geo")
DATA_PROC_DIR <- file.path(PROJECT_ROOT, "data", "processed")
RESULTS_TAB_DIR <- file.path(PROJECT_ROOT, "results", "tables")

cat("=== XELOX Resistance Study: Gene Pool Integration (v2 - Gene Symbol Aware) ===\n\n")

# ============================================================
# 1. Capecitabine metabolism pathway genes
# ============================================================
cat("[1] Loading capecitabine metabolism gene set\n")

# Comprehensive list of capecitabine/5-FU metabolism genes
cape_metabolism_genes <- c(
  # --- Activating enzymes ---
  "CES1", "CES2",       # Carboxylesterase - capecitabine -> 5'-DFCR
  "CDA",                # Cytidine deaminase - 5'-DFCR -> 5'-DFUR
  "TYMP",               # Thymidine phosphorylase (also TP) - 5'-DFUR -> 5-FU
  "UCK1", "UCK2",       # Uridine-cytidine kinase
  "TK1",                # Thymidine kinase 1
  
  # --- Catabolic/deactivating enzymes ---
  "DPYD",               # Dihydropyrimidine dehydrogenase - 5-FU -> DHFU (rate-limiting)
  "DPYS",               # Dihydropyrimidinase
  "UPB1",               # Beta-ureidopropionase
  
  # --- Target enzymes ---
  "TYMS",               # Thymidylate synthase - PRIMARY TARGET of 5-FU
  "DHFR",               # Dihydrofolate reductase
  
  # --- 5-FU nucleotide metabolism ---
  "RRM1", "RRM2",       # Ribonucleotide reductase
  "RRM2B",
  "NME1", "NME2",       # Nucleoside diphosphate kinase
  "UMPS",               # Uridine monophosphate synthetase
  
  # --- DNA repair / resistance ---
  "ERCC1", "ERCC2",     # NER pathway - oxaliplatin resistance
  "XRCC1",              # Base excision repair
  "GSTP1",              # Glutathione S-transferase - detoxification
  "ABCC2", "ABCC3",     # ATP-binding cassette transporters
  "ABCG2",
  
  # --- Oxaliplatin metabolism ---
  "ATP7A", "ATP7B",     # Copper transporters - oxaliplatin efflux
  "GSTT1", "GSTM1",     # Glutathione conjugation
  
  # --- Folate metabolism ---
  "MTHFR", "MTR", "MTRR",  # Folate cycle
  "SLC19A1",            # Reduced folate carrier
  
  # --- Cell cycle / apoptosis ---
  "TP53", "CDKN1A",     # p53 pathway
  "BAX", "BCL2",        # Apoptosis
  "BIRC5",              # Survivin
  
  # --- Wnt / stemness ---
  "CTNNB1", "MYC", "AXIN2",  # Wnt signaling
  
  # --- Additional oxaliplatin resistance ---
  "MGMT",               # DNA repair
  "MLH1", "MSH2", "MSH6",  # Mismatch repair
  
  # --- EMT markers ---
  "CDH1", "CDH2", "VIM", "SNAI1", "ZEB1",
  
  # --- Immune microenvironment ---
  "CD274", "PDCD1", "CTLA4"
)

cat(sprintf("  Capecitabine metabolism genes defined: %d\n", length(cape_metabolism_genes)))
cape_upper <- toupper(cape_metabolism_genes)

# ============================================================
# 2. Load MAPPED DEG results (gene symbols from Step 6)
# ============================================================
cat("\n[2] Loading mapped DEG results (gene symbols)\n")

deg_files <- list.files(RESULTS_TAB_DIR, pattern = "^DEG_.*_mapped\\.csv$", full.names = TRUE)

all_deg <- list()
for (f in deg_files) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  
  # Extract GSE ID from filename (e.g., "DEG_GSE39582_mapped" -> "GSE39582")
  base_name <- basename(f)
  gse_id <- gsub("^DEG_|_mapped\\.csv$", "", base_name)
  
  # Filter to rows with valid gene symbols
  df <- df[!is.na(df$gene_symbol) & df$gene_symbol != "", ]
  
  # Collapse duplicates: keep the row with minimum P.Value per gene symbol
  df <- df[order(df$P.Value), ]
  df <- df[!duplicated(df$gene_symbol), ]
  
  all_deg[[gse_id]] <- df
  cat(sprintf("  %s: %d unique gene symbols with P.Value\n", gse_id, nrow(df)))
}

# ============================================================
# 3. Load GPL570 annotation for WGCNA probe->gene mapping
# ============================================================
cat("\n[3] Loading GPL570 annotation for WGCNA gene mapping\n")

parse_gpl_annot <- function(filepath) {
  # Parse .annot.gz file: skip header lines before !platform_table_begin
  all_lines <- readLines(gzfile(filepath), warn = FALSE)
  data_start <- grep("^!platform_table_begin", all_lines)
  
  if (length(data_start) == 0) {
    stop("Could not find !platform_table_begin in ", filepath)
  }
  
  # Header is the line right after !platform_table_begin
  header_row <- data_start[1] + 1
  header <- strsplit(all_lines[header_row], "\t")[[1]]
  
  # Data rows
  data_rows <- all_lines[(header_row + 1):length(all_lines)]
  # Remove trailing !platform_table_end if present
  data_rows <- data_rows[!grepl("^!platform_table_end", data_rows)]
  
  # Parse
  con <- textConnection(paste(data_rows, collapse = "\n"))
  tbl <- read.delim(con, header = FALSE, stringsAsFactors = FALSE,
                    quote = "", check.names = FALSE)
  close(con)
  
  colnames(tbl) <- header
  return(tbl)
}

gpl570_file <- file.path(DATA_GEO_DIR, "GPL570.annot.gz")
if (file.exists(gpl570_file)) {
  gpl570_tbl <- parse_gpl_annot(gpl570_file)
  
  # Extract probe_id -> gene_symbol mapping
  probe_map <- data.frame(
    probe_id = toupper(trimws(gpl570_tbl$ID)),
    gene_symbol = trimws(gpl570_tbl$`Gene symbol`),
    stringsAsFactors = FALSE
  )
  
  # For multi-gene probes (e.g., "MIR4640///DDR1"), keep first symbol
  probe_map$gene_symbol <- gsub(" /// .*$", "", probe_map$gene_symbol)
  
  # Remove empty/unmapped
  probe_map <- probe_map[!is.na(probe_map$gene_symbol) & probe_map$gene_symbol != "" & probe_map$gene_symbol != "---", ]
  probe_map <- probe_map[!duplicated(probe_map$probe_id), ]
  
  cat(sprintf("  GPL570: %d probes mapped to gene symbols\n", nrow(probe_map)))
} else {
  cat("  WARNING: GPL570.annot.gz not found, WGCNA genes will remain as probe IDs\n")
  probe_map <- NULL
}

# ============================================================
# 4. Load WGCNA results and map probes to gene symbols
# ============================================================
cat("\n[4] Loading WGCNA results (with probe->gene mapping)\n")

# Hub genes
hub_file <- file.path(RESULTS_TAB_DIR, "WGCNA_hub_genes.csv")
hub_genes <- NULL
hub_symbols <- character(0)
if (file.exists(hub_file)) {
  hub_genes <- read.csv(hub_file, stringsAsFactors = TRUE)
  hub_probes <- toupper(trimws(as.character(hub_genes$Gene)))
  
  if (!is.null(probe_map)) {
    # Map hub gene probes to symbols
    hub_map <- probe_map[probe_map$probe_id %in% hub_probes, ]
    hub_symbols <- unique(hub_map$gene_symbol)
    
    # Report unmapped probes
    unmapped <- setdiff(hub_probes, probe_map$probe_id)
    cat(sprintf("  WGCNA hub genes: %d total probes, %d mapped to symbols, %d unmapped\n",
                length(hub_probes), length(hub_symbols), length(unmapped)))
    if (length(unmapped) > 0) {
      cat(sprintf("  Unmapped hub probes (first 10): %s\n",
                  paste(head(unmapped, 10), collapse = ", ")))
    }
  } else {
    hub_symbols <- hub_probes
    cat(sprintf("  WGCNA hub genes: %d (probe IDs, no mapping available)\n", length(hub_probes)))
  }
}

# All module membership (for WGCNA module gene set)
mm_file <- file.path(RESULTS_TAB_DIR, "WGCNA_all_gene_module_membership.csv")
wgcna_all <- NULL
wgcna_module_symbols <- character(0)
if (file.exists(mm_file)) {
  wgcna_all <- read.csv(mm_file, stringsAsFactors = TRUE)
  
  # Define which modules are biologically meaningful (exclude grey = unassigned)
  # From Step 4: modules found: turquoise(1), blue(2), brown(3), yellow(4)
  meaningful_modules <- c(1, 2, 3, 4)
  wgcna_in_modules <- wgcna_all[wgcna_all$ModuleLabel %in% meaningful_modules, ]
  
  module_probes <- toupper(trimws(as.character(wgcna_in_modules$Gene)))
  
  if (!is.null(probe_map)) {
    mod_map <- probe_map[probe_map$probe_id %in% module_probes, ]
    wgcna_module_symbols <- unique(mod_map$gene_symbol)
    cat(sprintf("  WGCNA module genes (excluding grey): %d probes, %d gene symbols\n",
                length(module_probes), length(wgcna_module_symbols)))
  } else {
    wgcna_module_symbols <- module_probes
    cat(sprintf("  WGCNA module genes (excluding grey): %d (probe IDs)\n", length(module_probes)))
  }
}

# ============================================================
# 5. Integrate gene pool using GENE SYMBOLS
# ============================================================
cat("\n[5] Integrating final gene pool (gene symbols)\n")

# Strategy:
# Tier 1 - DEG gene_symbol ∩ WGCNA module gene_symbol ∩ Capecitabine metabolism
# Tier 2 - DEG gene_symbol ∩ WGCNA module gene_symbol (excluding Tier 1)
# Tier 3 - DEG gene_symbol ∩ Capecitabine metabolism (excluding Tier 1/2)
# Tier 4 - Gene symbols appearing as DEG in 2+ datasets (excluding Tier 1-3)
# Tier 5 - WGCNA hub gene symbols (excluding Tier 1-4)

# Collect significant DEG gene symbols (P < 0.05 from any dataset)
deg_genes_list <- list()
for (name in names(all_deg)) {
  df <- all_deg[[name]]
  sig_genes <- toupper(trimws(df$gene_symbol[df$P.Value < 0.05]))
  sig_genes <- sig_genes[!is.na(sig_genes) & sig_genes != ""]
  deg_genes_list[[name]] <- sig_genes
  cat(sprintf("  %s: %d nominal DEGs (gene symbols)\n", name, length(sig_genes)))
}

# Union of all DEG gene symbols
all_deg_symbols <- unique(unlist(deg_genes_list))
cat(sprintf("\n  Union of DEG gene symbols (P<0.05): %d\n", length(all_deg_symbols)))

# --- Tier calculations using gene symbols ---

# Tier 1: DEG ∩ WGCNA module ∩ Capecitabine metabolism
tier1 <- intersect(intersect(all_deg_symbols, wgcna_module_symbols), cape_upper)
cat(sprintf("  Tier 1 (DEG ∩ WGCNA ∩ Metabolism): %d genes\n", length(tier1)))
if (length(tier1) > 0) {
  cat("    ")
  cat(tier1, sep = ", ")
  cat("\n")
}

# Tier 2: DEG ∩ WGCNA module (excluding Tier 1)
tier2 <- intersect(all_deg_symbols, wgcna_module_symbols)
tier2 <- setdiff(tier2, tier1)
cat(sprintf("  Tier 2 (DEG ∩ WGCNA, excl Tier 1): %d genes\n", length(tier2)))
if (length(tier2) > 0) {
  cat("    First 20: ")
  cat(head(tier2, 20), sep = ", ")
  cat("\n")
}

# Tier 3: DEG ∩ Metabolism (outside WGCNA)
tier3 <- intersect(all_deg_symbols, cape_upper)
tier3 <- setdiff(tier3, c(tier1, tier2))
cat(sprintf("  Tier 3 (DEG ∩ Metabolism, excl WGCNA): %d genes\n", length(tier3)))
if (length(tier3) > 0) {
  cat("    ")
  cat(tier3, sep = ", ")
  cat("\n")
}

# Tier 4: Gene symbols appearing in 2+ DEG datasets
deg_freq <- table(unlist(deg_genes_list))
multi_deg <- names(deg_freq[deg_freq >= 2])
tier4 <- setdiff(multi_deg, c(tier1, tier2, tier3))
cat(sprintf("  Tier 4 (Multi-dataset DEGs, ≥2): %d genes\n", length(tier4)))
if (length(tier4) > 0) {
  cat("    First 20: ")
  cat(head(tier4, 20), sep = ", ")
  cat("\n")
}

# Tier 5: WGCNA hub gene symbols (not in above)
tier5 <- setdiff(hub_symbols, c(tier1, tier2, tier3, tier4))
cat(sprintf("  Tier 5 (WGCNA hub only): %d genes\n", length(tier5)))
if (length(tier5) > 0) {
  cat("    First 20: ")
  cat(head(tier5, 20), sep = ", ")
  cat("\n")
}

# Final prioritized gene pool
final_pool <- data.frame(
  Gene = c(tier1, tier2, tier3, tier4, tier5),
  Priority = c(rep("Tier1_Capecitabine_DEG_WGCNA", length(tier1)),
               rep("Tier2_DEG_WGCNA", length(tier2)),
               rep("Tier3_DEG_Metabolism", length(tier3)),
               rep("Tier4_MultiDataset_DEG", length(tier4)),
               rep("Tier5_WGCNA_Hub", length(tier5))),
  stringsAsFactors = FALSE
)

cat(sprintf("\n  FINAL GENE POOL: %d genes\n", nrow(final_pool)))
cat("  Priority breakdown:\n")
print(table(final_pool$Priority))

write.csv(final_pool, file.path(RESULTS_TAB_DIR, "XELOX_final_gene_pool.csv"), row.names = FALSE)
cat("  Gene pool saved\n")

# ============================================================
# 6. Cross-reference with capecitabine metabolism
# ============================================================
cat("\n[6] Capecitabine metabolism cross-reference\n")

cape_in_deg <- cape_upper[cape_upper %in% all_deg_symbols]
cape_in_wgcna <- cape_upper[cape_upper %in% wgcna_module_symbols]
cape_in_both <- intersect(cape_in_deg, cape_in_wgcna)
cape_in_hub <- cape_upper[cape_upper %in% hub_symbols]

cat(sprintf("  Metabolism genes in DEG results: %d / %d\n",
            length(cape_in_deg), length(cape_metabolism_genes)))
cat(sprintf("  Metabolism genes in WGCNA modules: %d / %d\n",
            length(cape_in_wgcna), length(cape_metabolism_genes)))
cat(sprintf("  Metabolism genes in WGCNA hub: %d / %d\n",
            length(cape_in_hub), length(cape_metabolism_genes)))
cat(sprintf("  Metabolism genes in BOTH DEG+WGCNA: %d\n", length(cape_in_both)))

if (length(cape_in_deg) > 0) {
  cat("  Metabolism genes in DEG:\n    ")
  cat(cape_in_deg, sep = ", ")
  cat("\n")
}

if (length(cape_in_wgcna) > 0) {
  cat("  Metabolism genes in WGCNA modules:\n    ")
  cat(cape_in_wgcna, sep = ", ")
  cat("\n")
}

if (length(cape_in_hub) > 0) {
  cat("  Metabolism genes in WGCNA hub:\n    ")
  cat(cape_in_hub, sep = ", ")
  cat("\n")
}

# ============================================================
# 7. Generate report
# ============================================================
cat("\n[7] Generating integration report\n")

report_lines <- c(
  "=================================================================",
  "XELOX Resistance: Final Gene Pool Integration Report (v2 - Gene Symbols)",
  "=================================================================",
  "",
  sprintf("Capecitabine metabolism genes defined: %d", length(cape_metabolism_genes)),
  sprintf("DEG gene symbols (union, P<0.05, mapped): %d", length(all_deg_symbols)),
  sprintf("WGCNA module gene symbols (excl grey): %d", length(wgcna_module_symbols)),
  sprintf("WGCNA hub gene symbols: %d", length(hub_symbols)),
  sprintf("DEG datasets contributing: %d", length(deg_genes_list)),
  "",
  "Tiered integration (gene symbols):",
  sprintf("  Tier 1 (DEG ∩ WGCNA ∩ Metabolism): %d genes", length(tier1)),
  sprintf("  Tier 2 (DEG ∩ WGCNA): %d genes", length(tier2)),
  sprintf("  Tier 3 (DEG ∩ Metabolism): %d genes", length(tier3)),
  sprintf("  Tier 4 (Multi-dataset DEG, >=2): %d genes", length(tier4)),
  sprintf("  Tier 5 (WGCNA Hub only): %d genes", length(tier5)),
  "",
  sprintf("FINAL GENE POOL: %d genes", nrow(final_pool)),
  "",
  "--- Tier 1 Genes (highest confidence) ---",
  if (length(tier1) > 0) paste(tier1, collapse = ", ") else "(none)",
  "",
  "--- Tier 2 Genes (DEG ∩ WGCNA) ---",
  if (length(tier2) > 0) paste(tier2, collapse = ", ") else "(none)",
  "",
  "--- Tier 3 Genes (DEG ∩ Metabolism) ---",
  if (length(tier3) > 0) paste(tier3, collapse = ", ") else "(none)",
  "",
  "--- Capecitabine Metabolism in DEG ---",
  if (length(cape_in_deg) > 0) paste(cape_in_deg, collapse = ", ") else "(none)",
  "",
  "--- Capecitabine Metabolism in WGCNA ---",
  if (length(cape_in_wgcna) > 0) paste(cape_in_wgcna, collapse = ", ") else "(none)",
  "",
  "--- Capecitabine Metabolism in Hub ---",
  if (length(cape_in_hub) > 0) paste(cape_in_hub, collapse = ", ") else "(none)"
)

writeLines(report_lines, file.path(RESULTS_TAB_DIR, "gene_pool_integration_report.txt"))

cat("\n=== Gene Pool Integration Complete ===\n")
