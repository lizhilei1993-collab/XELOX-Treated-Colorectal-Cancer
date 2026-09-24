out <- file("/path/to/xelox_reanalysis/scale_check.txt", open="wt", encoding="UTF-8")
ow <- function(...) writeLines(paste0(...), con=out)
for (g in c("GSE39582","GSE104645","GSE28702","GSE72970","GSE69657")) {
  f <- paste0("/path/to/xelox_reanalysis/in/", g, "_m.rds")
  x <- tryCatch(readRDS(f), error=function(e){ow(g," readRDS failed: ",conditionMessage(e)); NULL})
  if (is.null(x)) next
  m <- NULL
  for (fn in c("as.matrix","exprs")) {
    cand <- tryCatch({
      z <- if (fn=="exprs") x else x
      if (fn=="exprs") as.matrix(exprs(z)) else as.matrix(z)
    }, error=function(e) NULL)
    if (!is.null(cand)) { m <- cand; break }
  }
  if (is.null(m)) { ow(sprintf("%-11s class=%s  not coercible to matrix", g, class(x)[1])); next }
  v <- as.numeric(m); v <- v[is.finite(v)]
  ow(sprintf("%-11s class=%-9s dims=%dx%-5d min=%-8.2f med=%-7.2f p99=%-8.2f max=%-9.2f  %s",
     g, class(x)[1], nrow(m), ncol(m), min(v), median(v),
     quantile(v,.99), max(v),
     ifelse(max(v) > 40, "<== NOT log2 (linear scale)", "log2-consistent")))
  rm(m, v); invisible(gc())
}
close(out); cat("done\n")
