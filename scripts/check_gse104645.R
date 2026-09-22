library_path <- "C:/Rlibs"
.libPaths(c(library_path, .libPaths()))
Sys.setenv(TMPDIR = "C:/temp", TMP = "C:/temp", TEMP = "C:/temp")

m <- readRDS("C:/xelox_research/data/geo/GSE104645_expression.rds")
rns <- rownames(m)
cat("Total rows:", length(rns), "\n")

# Get all non-control probe IDs
non_ctrl <- rns[!grepl("^\\(\\+\\)", rns)]
cat("Non-control probes:", length(non_ctrl), "\n")
cat("First 15 non-control:", paste(non_ctrl[1:15], collapse=", "), "\n\n")

# Check if they match GPL6480 annotation
cat("=== Loading GPL6480 annotation ===\n")
annot <- read.table("C:/xelox_research/data/geo/GPL6480.annot.gz", 
  header=TRUE, sep="\t", comment.char="!", quote="\"", stringsAsFactors=FALSE)
cat("GPL6480 annot rows:", nrow(annot), "\n")
cat("GPL6480 annot first probe IDs:", paste(annot$ID[1:10], collapse=", "), "\n\n")

# How many GSE104645 probes match GPL6480?
match_count <- sum(non_ctrl %in% annot$ID)
cat("GSE104645 probes matching GPL6480:", match_count, "/", length(non_ctrl), "\n")

# Get gene symbols for matched probes
idx <- match(non_ctrl, annot$ID)
matched_genes <- annot$Gene.symbol[idx[!is.na(idx)]]
cat("Unique gene symbols from matched probes:", length(unique(matched_genes[matched_genes != ""])), "\n")

# Also check GSE39582
cat("\n=== GSE39582 ===\n")
e <- readRDS("C:/xelox_research/data/geo/GSE39582_eset.rds")
cat("Class:", class(e), "\n")
require(Biobase, quietly=TRUE, lib.loc=library_path)
mat <- exprs(e)
cat("Expr dim:", nrow(mat), "x", ncol(mat), "\n")
cat("Rownames[1:5]:", paste(rownames(mat)[1:5], collapse=", "), "\n")
