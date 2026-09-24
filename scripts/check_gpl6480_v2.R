library_path <- "/path/to/Rlibs"
.libPaths(c(library_path, .libPaths()))
Sys.setenv(TMPDIR = "/tmp", TMP = "/tmp", TEMP = "/tmp")

# Use scanning to skip header lines
cat("=== Reading GPL6480 annotation with header skip ===\n")
con <- gzfile("/path/to/xelox_project/data/geo/GPL6480.annot.gz", "rt")
# Read until we hit the table header line (after !platform_table_begin)
while (TRUE) {
  line <- readLines(con, n=1)
  if (length(line) == 0 || grepl("!platform_table_begin", line)) break
}
# Now read the column names
col_names <- readLines(con, n=1)
cat("Raw col_names line:", col_names, "\n")
col_names <- strsplit(col_names, "\t")[[1]]
cat("Parsed col_names:", paste(col_names, collapse=" | "), "\n")

# Read first 5 data rows
data_lines <- readLines(con, n=5)
close(con)

# Parse
data <- read.table(text = data_lines, sep = "\t", header = FALSE, 
  col.names = col_names, quote = "\"", stringsAsFactors = FALSE)
cat("\nFirst data row:\n")
print(data[1, c("ID", "Gene.symbol", "Gene.title")])

# Now check matching
m <- readRDS("/path/to/xelox_project/data/geo/GSE104645_expression.rds")
rns <- rownames(m)
non_ctrl <- rns[!grepl("^\\(\\+\\)", rns)]

cat("\n=== Probe matching ===\n")
idx <- match(non_ctrl, data$ID)
cat("Matched in first 5 rows:", sum(!is.na(idx)), "\n")

# Full annot read
cat("\nReading full annotation...\n")
annot_full <- read.table("/path/to/xelox_project/data/geo/GPL6480.annot.gz", 
  sep="\t", skip=29, header=TRUE, quote="\"", stringsAsFactors=FALSE)
cat("Full annot rows:", nrow(annot_full), "\n")
cat("Colnames:", paste(colnames(annot_full), collapse=" | "), "\n")

idx_full <- match(non_ctrl, annot_full$ID)
cat("GSE104645 probes matched in GPL6480:", sum(!is.na(idx_full)), "/", length(non_ctrl), "\n")

# Get gene symbols
matched_symbols <- annot_full$Gene.symbol[idx_full[!is.na(idx_full)]]
cat("Unique gene symbols:", length(unique(matched_symbols[matched_symbols != ""])), "\n")
