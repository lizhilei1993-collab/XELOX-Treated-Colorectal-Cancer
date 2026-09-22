library_path <- "C:/Rlibs"
.libPaths(c(library_path, .libPaths()))
Sys.setenv(TMPDIR = "C:/temp", TMP = "C:/temp", TEMP = "C:/temp")

# Check GPL6480 annotation format
cat("=== GPL6480 annotation raw header ===\n")
con <- gzfile("C:/xelox_research/data/geo/GPL6480.annot.gz", "rt")
lines <- readLines(con, n=30)
close(con)
for (l in lines) cat(l, "\n")

cat("\n=== Trying different read.table ===\n")
annot <- read.table("C:/xelox_research/data/geo/GPL6480.annot.gz", 
  header=TRUE, sep="\t", comment.char="", quote="\"", stringsAsFactors=FALSE, nrows=10)
cat("Colnames:", paste(colnames(annot), collapse=" | "), "\n")
cat("First ID col:", annot[1,1], "\n")
cat("Gene symbol col:", annot$Gene.symbol[1], "\n")

# Now check GSE104645 probe ID matching
m <- readRDS("C:/xelox_research/data/geo/GSE104645_expression.rds")
rns <- rownames(m)
non_ctrl <- rns[!grepl("^\\(\\+\\)", rns)]
cat("\nGSE104645 first probe ID:", non_ctrl[2], "\n")
cat("Matches first annot ID:", non_ctrl[2] == annot$ID[1], "\n")

# Check matching properly
idx <- match(non_ctrl, annot$ID)
cat("Matched:", sum(!is.na(idx)), "/", length(non_ctrl), "\n")
