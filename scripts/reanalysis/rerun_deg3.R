# Corrected five-cohort DEG re-run, on a verified log2 scale, with an
# explicit control so the repaired cohorts can be trusted.
#
# All paths are ASCII. R on this machine cannot open files under the Chinese
# account/project path (LC_CTYPE=C.UTF-8 fails at startup), so every input was
# staged into /path/to/xelox_reanalysis/in beforehand.
#
# The sample->resistance labelling is NOT re-derived here: it is taken from the
# project's own annotate_* functions (verbatim lift of
# 15_pathway_differential_analysis.R lines 55-191) applied to the same
# pathway-score objects the published pathway analysis used. So the cohort
# subsets are identical by construction, and their sizes are printed to be
# checked against the recorded 113 / 83 / 124 / 30 / 164.

Sys.setenv(TMPDIR = "/tmp", TEMP = "/tmp", TMP = "/tmp")
.libPaths(c("/path/to/Rlibs", .libPaths()))
suppressPackageStartupMessages({ library(limma); library(Biobase) })

IN  <- "/path/to/xelox_reanalysis/in/"
OUTD <- Sys.getenv("RERUN_OUT", "/path/to/xelox_reanalysis/rerun/")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
lg <- file(paste0(OUTD, "rerun_log.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) writeLines(sprintf(...), con = lg)

src <- readLines(paste0(IN, "annot_source.R"), encoding = "UTF-8", warn = FALSE)
i <- grep("Defining sample annotation functions", src)[1]
j <- grep("Defining differential analysis function", src)[1]
eval(parse(text = src[(i + 2):(j - 2)]), envir = globalenv())
say("lifted annotators: %s", paste(grep("^annotate_", ls(), value = TRUE), collapse = ", "))

to_matrix <- function(x) if (is(x, "ExpressionSet")) as.matrix(exprs(x)) else as.matrix(x)

analyze_deg_limma <- function(expr_mat, group_vec) {
  group_f <- factor(group_vec, levels = c("sensitive", "resistant"))
  design  <- model.matrix(~ group_f)
  if (Sys.getenv("SKIP_FILTER") == "1") {
    keep <- rep(TRUE, nrow(expr_mat))          # attribution arm: filter pinned off
  } else {
    keep <- apply(expr_mat, 1, median, na.rm = TRUE) >= median(expr_mat, na.rm = TRUE) &
            apply(expr_mat, 1, IQR, na.rm = TRUE) >= 0.5
  }
  e <- normalizeBetweenArrays(as.matrix(expr_mat[keep, , drop = FALSE]), method = "quantile")
  res <- topTable(eBayes(lmFit(e, design)), coef = 2, number = Inf, sort.by = "P")
  res$probe <- rownames(res)                      # the column 03 dropped
  res$Direction <- ifelse(res$logFC > 0, "up", "down")
  list(results = res, tested = nrow(e))
}

as_log2 <- function(m, cohort) {
  v <- as.numeric(m); v <- v[is.finite(v)]
  if (max(v) > 40) {
    if (Sys.getenv("SKIP_REPAIR") == "1") {
      say("%s   : max %-10.1f  LINEAR scale LEFT UNREPAIRED (isolation arm)", cohort, max(v))
      return(m)
    }
    if (min(v) < 0) { say("%s: linear scale with negative values, refusing to log2", cohort); return(NULL) }
    say("%s   : max %-10.1f  LINEAR -> log2(x+1) applied", cohort, max(v))
    return(log2(m + 1))
  }
  say("%s   : max %-10.2f  already log2, untouched", cohort, max(v))
  m
}

cohorts <- c("GSE39582", "GSE104645", "GSE28702", "GSE72970", "GSE69657")
res_store <- new.env()

for (cohort in cohorts) {
  m  <- to_matrix(readRDS(paste0(IN, cohort, "_m.rds")))
  ps <- readRDS(paste0(IN, cohort, "_pathway_scores.rds"))
  cf <- if (cohort == "GSE39582") paste0(IN, "GSE39582_xelox_groups.csv")
        else paste0(IN, cohort, "_clinical_data.csv")
  ann <- tryCatch(get(paste0("annotate_", cohort))(ps, cf),
                  error = function(e) { say("%s: annotator failed: %s", cohort, conditionMessage(e)); NULL })
  if (is.null(ann)) next

  g   <- tolower(ann$group)     # annotators return "Sensitive"/"Resistant"; the
  ids <- colnames(ann$pathway_matrix)   # limma design in 03 uses lowercase levels
  if (is.null(ids) || is.null(g) || length(ids) != length(g)) {
    say("%s: annotator output unusable (ids=%s group=%s)", cohort, length(ids), length(g)); next
  }
  if (is.null(names(g))) names(g) <- ids        # index by sample id from here on
  keep <- !is.na(g) & g %in% c("resistant", "sensitive")
  ids <- intersect(ids[keep], colnames(m))
  g <- g[ids]
  say("%s: labelled samples %d -> matched to expression columns %d (%d resistant / %d sensitive)",
      cohort, sum(keep), length(ids), sum(g == "resistant"), sum(g == "sensitive"))
  if (length(ids) < 10) { say("%s: too few matched samples, skipping", cohort); next }

  mm <- as_log2(m[, ids, drop = FALSE], cohort)
  if (is.null(mm)) next
  r <- analyze_deg_limma(mm, g[match(ids, names(g))])
  assign(cohort, r$results, envir = res_store)
  write.csv(r$results, paste0(OUTD, "DEG_", cohort, "_log2corrected.csv"), row.names = FALSE)
  say("%s   : %d probes tested | adj.P<0.05 = %d (%d up / %d down) | max |logFC| = %.3f",
      cohort, r$tested, sum(r$results$adj.P.Val < .05, na.rm = TRUE),
      sum(r$results$adj.P.Val < .05 & r$results$logFC > 0, na.rm = TRUE),
      sum(r$results$adj.P.Val < .05 & r$results$logFC < 0, na.rm = TRUE),
      max(abs(r$results$logFC), na.rm = TRUE))
  rm(m, mm, ps); invisible(gc())
}

saveRDS(as.list(res_store), paste0(OUTD, "corrected_deg_results.rds"))
say("\nDONE. corrected tables in %s", OUTD)
close(lg)
cat("done -", paste0(OUTD, "rerun_log.txt"), "\n")
