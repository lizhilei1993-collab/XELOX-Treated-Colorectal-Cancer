# ============================================================
# 06_PROBE_TO_GENE_MAPPING.R
# Map Affymetrix probe IDs to gene symbols across all datasets
# XELOX Resistance Study
# ============================================================

Sys.setenv(TMPDIR = "C:/temp", TMP = "C:/temp", TEMP = "C:/temp")
.libPaths(c("C:/Rlibs", .libPaths()))

library(GEOquery)

# Quiet warnings
options(warn = -1)

PROJECT_ROOT <- "C:/xelox_research"
DATA_GEO_DIR  <- file.path(PROJECT_ROOT, "data", "geo")
DATA_PROC_DIR <- file.path(PROJECT_ROOT, "data", "processed")
RESULTS_TAB_DIR <- file.path(PROJECT_ROOT, "results", "tables")

cat("============================================================\n")
cat("Probe ID to Gene Symbol Mapping\n")
cat("============================================================\n\n")

# ============================================================
# 1. Identify which GPL platforms each dataset uses
# ============================================================
cat("[1] Checking dataset platforms\n\n")

# Series matrix files - extract platform info
sm_files <- list.files(DATA_GEO_DIR, pattern = "_series_matrix.txt.gz$", full.names = TRUE)

dataset_platforms <- list()

for (sm_file in sm_files) {
  gse_id <- gsub("_series_matrix.txt.gz", "", basename(sm_file))
  
  # Read just the first 200 lines to find platform
  lines <- readLines(sm_file, n = 200, warn = FALSE)
  platform_line <- grep("^!Series_platform_id", lines, value = TRUE)
  
  if (length(platform_line) > 0) {
    gpl_id <- trimws(gsub("^!Series_platform_id = ", "", platform_line[1]))
    dataset_platforms[[gse_id]] <- gpl_id
    cat(sprintf("  %s -> %s\n", gse_id, gpl_id))
  } else {
    cat(sprintf("  %s -> platform not found\n", gse_id))
  }
}

cat("\n")

# ============================================================
# 2. Parse GPL annotation files to create probe-gene maps
# ============================================================
cat("[2] Parsing GPL annotation files\n\n")

# Find all GPL files
gpl_files <- list.files(DATA_GEO_DIR, pattern = "\\.soft\\.gz$|\\.annot\\.gz$", full.names = TRUE)

# Known mappings: probe ID -> gene symbol
gpl_maps <- list()

for (gpl_file in gpl_files) {
  gpl_id <- gsub("\\.soft\\.gz$|\\.annot\\.gz$", "", basename(gpl_file))
  cat(sprintf("  Parsing %s...\n", gpl_id))
  
  # Read the GPL file
  # For .soft.gz, use GEOquery's parseGEO
  # For .annot.gz, read as table
  if (grepl("\\.soft\\.gz$", gpl_file)) {
    tryCatch({
      gpl <- parseGEO(gpl_file, GSElimits = NULL, getGPL = FALSE)
      
      # Extract probe-to-gene mapping
      if (is(gpl, "data.frame")) {
        tbl <- gpl
      } else if (is(gpl, "GPL")) {
        tbl <- Table(gpl)
      } else {
        cat(sprintf("    Unexpected type: %s\n", class(gpl)[1]))
        next
      }
      
      # Find ID and gene symbol columns
      id_col <- grep("^ID$|^Probe.Set.ID$|^Probe ID$|^ID_REF$", colnames(tbl), ignore.case = TRUE, value = TRUE)
      sym_cols <- grep("symbol|gene|gene.symbol|Gene.Symbol|Gene.symbol|GENE_SYMBOL|Symbol", 
                       colnames(tbl), ignore.case = FALSE, value = TRUE)
      
      if (length(id_col) > 0 && length(sym_cols) > 0) {
        # Create mapping: keep first gene symbol per probe
        mapping <- tbl[, c(id_col[1], sym_cols[1])]
        colnames(mapping) <- c("probe_id", "gene_symbol")
        
        # Clean: remove empty/NA symbols, split multiple symbols
        mapping$gene_symbol <- trimws(as.character(mapping$gene_symbol))
        mapping <- mapping[!is.na(mapping$gene_symbol) & mapping$gene_symbol != "" & mapping$gene_symbol != "---", ]
        
        # For probes with multiple gene symbols (e.g., "EGFR // EGF"), keep first
        mapping$gene_symbol <- gsub(" // .*$", "", mapping$gene_symbol)
        mapping$gene_symbol <- gsub(" / .*$", "", mapping$gene_symbol)
        
        # Remove duplicates (keep first)
        mapping <- mapping[!duplicated(mapping$probe_id), ]
        
        gpl_maps[[gpl_id]] <- mapping
        cat(sprintf("    Mapped %d probes to genes\n", nrow(mapping)))
        
        # Show a few examples
        cat(sprintf("    Examples: %s -> %s, %s -> %s, %s -> %s\n",
                    mapping$probe_id[1], mapping$gene_symbol[1],
                    mapping$probe_id[2], mapping$gene_symbol[2],
                    mapping$probe_id[3], mapping$gene_symbol[3]))
      } else {
        cat(sprintf("    No gene symbol column found. ID col: %s, Sym cols: %s\n",
                    paste(id_col, collapse=","), paste(sym_cols, collapse=",")))
        cat("    Available columns:", paste(colnames(tbl)[1:min(20, ncol(tbl))], collapse=", "), "\n")
      }
    }, error = function(e) {
      cat(sprintf("    ERROR: %s\n", conditionMessage(e)))
    })
  } else if (grepl("\\.annot\\.gz$", gpl_file)) {
    # .annot.gz files from GEO - need to skip header lines before !platform_table_begin
    tryCatch({
      # Read raw lines, find where data starts
      all_lines <- readLines(gzfile(gpl_file), warn = FALSE)
      data_start <- grep("^!platform_table_begin", all_lines)
      
      # Skip header lines (comment markers # and meta lines ^/!)
      # The actual data starts after !platform_table_begin + 1 (the header row that follows it)
      if (length(data_start) > 0) {
        # Find the header row (after !platform_table_begin)
        header_row <- data_start[1] + 1
        data_rows <- all_lines[(header_row + 1):length(all_lines)]
        
        # Parse header
        header <- strsplit(all_lines[header_row], "\t")[[1]]
        
        # Parse data
        con <- textConnection(paste(data_rows, collapse = "\n"))
        tbl <- read.delim(con, header = FALSE, stringsAsFactors = FALSE, 
                          quote = "", check.names = FALSE)
        close(con)
        
        if (ncol(tbl) == length(header)) {
          colnames(tbl) <- header
        } else {
          cat(sprintf("    Column mismatch: data=%d, header=%d\n", ncol(tbl), length(header)))
          # Try with first row as header
          tbl <- read.delim(con, stringsAsFactors = FALSE, check.names = FALSE)
        }
      } else {
        # Fallback: skip lines starting with ^, !, #
        clean_lines <- grep("^[^!#^]", all_lines, value = TRUE)
        con <- textConnection(paste(clean_lines, collapse = "\n"))
        tbl <- read.delim(con, stringsAsFactors = FALSE, check.names = FALSE)
        close(con)
      }
      
      cat(sprintf("    Read %d rows, columns: %s\n", nrow(tbl), paste(colnames(tbl)[1:min(8, ncol(tbl))], collapse=", ")))
      
      id_col <- grep("^ID$|Probe.Set.ID|Probe ID|ID_REF", colnames(tbl), ignore.case = TRUE, value = TRUE)
      sym_cols <- grep("Symbol|Gene.Symbol|GENE_SYMBOL|Gene.symbol|gene.symbol", 
                       colnames(tbl), value = TRUE)
      
      if (length(id_col) > 0 && length(sym_cols) > 0) {
        mapping <- tbl[, c(id_col[1], sym_cols[1])]
        colnames(mapping) <- c("probe_id", "gene_symbol")
        mapping$gene_symbol <- trimws(as.character(mapping$gene_symbol))
        mapping <- mapping[!is.na(mapping$gene_symbol) & mapping$gene_symbol != "" & mapping$gene_symbol != "---", ]
        mapping$gene_symbol <- gsub(" // .*$", "", mapping$gene_symbol)
        mapping <- mapping[!duplicated(mapping$probe_id), ]
        gpl_maps[[gpl_id]] <- mapping
        cat(sprintf("    Mapped %d probes to genes\n", nrow(mapping)))
      } else {
        cat(sprintf("    No gene symbol column. ID: %s, Sym: %s\n",
                    paste(id_col, collapse=","), paste(sym_cols, collapse=",")))
      }
    }, error = function(e) {
      cat(sprintf("    ERROR: %s\n", conditionMessage(e)))
    })
  }
}

cat("\n")

# ============================================================
# 3. Apply mappings to each DEG dataset
# ============================================================
cat("[3] Applying probe-to-gene mappings to DEG datasets\n\n")

# For each dataset, determine which GPL it uses, then map
deg_files <- list.files(RESULTS_TAB_DIR, pattern = "^DEG_.*_limma\\.csv$", full.names = TRUE)

all_mapped <- list()

for (deg_file in deg_files) {
  base_name <- gsub("\\.csv$", "", basename(deg_file))
  # Extract GSE ID
  gse_id <- gsub("^DEG_|_limma$", "", base_name)
  
  cat(sprintf("  Processing %s...\n", base_name))
  
  # Find which platform this dataset uses
  gpl_id <- dataset_platforms[[gse_id]]
  if (is.null(gpl_id)) {
    cat(sprintf("    WARNING: No platform found for %s, trying all GPL maps\n", gse_id))
  }
  
  # Load DEG results
  deg <- read.csv(deg_file, stringsAsFactors = FALSE)
  cat(sprintf("    Loaded %d genes (probe IDs)\n", nrow(deg)))
  
  # Try to map
  mapped_deg <- NULL
  
  if (!is.null(gpl_id) && gpl_id %in% names(gpl_maps)) {
    # Use specific platform map
    mapping <- gpl_maps[[gpl_id]]
    deg$probe_id <- toupper(trimws(deg$Gene))
    mapping$probe_id <- toupper(trimws(mapping$probe_id))
    
    mapped_deg <- merge(deg, mapping[, c("probe_id", "gene_symbol")], by = "probe_id", all.x = TRUE)
    cat(sprintf("    Platform %s: %d probes mapped to genes\n", gpl_id, 
                sum(!is.na(mapped_deg$gene_symbol))))
  } else {
    # Try ALL platform maps
    deg$probe_id <- toupper(trimws(deg$Gene))
    
    for (map_name in names(gpl_maps)) {
      mapping <- gpl_maps[[map_name]]
      mapping$probe_id <- toupper(trimws(mapping$probe_id))
      
      hits <- sum(deg$probe_id %in% mapping$probe_id)
      cat(sprintf("    Trying %s: %d probes match\n", map_name, hits))
      
      if (hits > 1000) {
        # This is likely the correct platform
        mapped_deg <- merge(deg, mapping[, c("probe_id", "gene_symbol")], by = "probe_id", all.x = TRUE)
        gpl_id <- map_name
        cat(sprintf("    -> Assigned platform %s: %d probes mapped\n", map_name,
                    sum(!is.na(mapped_deg$gene_symbol))))
        break
      }
    }
  }
  
  if (!is.null(mapped_deg)) {
    # Collapse multiple probes per gene: keep the one with lowest P.Value
    mapped_deg <- mapped_deg[order(mapped_deg$P.Value), ]
    mapped_deg <- mapped_deg[!duplicated(mapped_deg$gene_symbol) | is.na(mapped_deg$gene_symbol), ]
    
    # Separate mapped and unmapped
    mapped_genes <- mapped_deg[!is.na(mapped_deg$gene_symbol) & mapped_deg$gene_symbol != "", ]
    unmapped_probes <- mapped_deg[is.na(mapped_deg$gene_symbol) | mapped_deg$gene_symbol == "", ]
    
    cat(sprintf("    After collapsing: %d unique gene symbols, %d unmapped probes\n",
                nrow(mapped_genes), nrow(unmapped_probes)))
    
    # Save mapped results
    out_file <- file.path(RESULTS_TAB_DIR, sprintf("DEG_%s_mapped.csv", gse_id))
    write.csv(mapped_genes, out_file, row.names = FALSE)
    cat(sprintf("    Saved: %s\n", basename(out_file)))
    
    all_mapped[[gse_id]] <- mapped_genes
    
    # Save unmapped probes for reference
    if (nrow(unmapped_probes) > 0) {
      unmapped_file <- file.path(RESULTS_TAB_DIR, sprintf("DEG_%s_unmapped.csv", gse_id))
      write.csv(unmapped_probes, unmapped_file, row.names = FALSE)
      cat(sprintf("    Unmapped probes saved: %s\n", basename(unmapped_file)))
    }
  } else {
    cat(sprintf("    WARNING: No mapping available for %s\n", gse_id))
  }
  
  cat("\n")
}

# ============================================================
# 4. Summary
# ============================================================
cat("============================================================\n")
cat("Mapping Summary\n")
cat("============================================================\n\n")

total_mapped <- 0
total_unmapped <- 0

for (gse_id in names(all_mapped)) {
  mapped_df <- all_mapped[[gse_id]]
  mapped_count <- nrow(mapped_df)
  # Count how many had gene symbols before collapsing
  total_mapped <- total_mapped + mapped_count
  
  cat(sprintf("  %s: %d gene symbols\n", gse_id, mapped_count))
}

cat(sprintf("\n  Total gene symbols across all datasets: %d\n", total_mapped))

save.image(file.path(DATA_PROC_DIR, "06_probe_mapping.RData"))
cat("\nWorkspace saved\n")
