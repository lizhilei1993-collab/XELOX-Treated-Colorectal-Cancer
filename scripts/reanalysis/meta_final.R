# Five-cohort Stouffer meta-analysis, recomputed on the corrected inputs.
#
# Inputs are the rerun_final/ tables, i.e. the Methods' IQR filter active and the
# two linear-scale cohorts log2-repaired. Everything downstream of limma is the
# published pipeline's own rule, verified rather than assumed:
#   probe -> symbol  : the project's own GPL570 map (data/geo) and the GPL6480
#                      annotation, both staged as ASCII CSVs
#   gene collapse    : keep the lowest-P probe per symbol, as 06 does
#   statistic        : Stouffer, weights sqrt(n) with the n each cohort actually
#                      contributed (113/83/124/30/164) - this weighting reproduces
#                      all 17,517 finite stouffer_z values in the published table
#                      exactly, so the corrected number is comparable to theirs
#   multiplicity     : Benjamini-Hochberg over all tested genes

IN  <- "/path/to/xelox_reanalysis/in/"
# which arm to meta-analyse: filter+repair (rerun_final) or repair-only (rerun_nofilt)
CR  <- Sys.getenv("META_TABLES", "/path/to/xelox_reanalysis/rerun_final/")
OUT <- Sys.getenv("META_OUT",   "/path/to/xelox_reanalysis/rerun_final/")
lg <- file(paste0(OUT, "meta_final_log.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { l <- sprintf(...); writeLines(l, con = lg); cat(l, "\n") }

COH <- c("GSE104645", "GSE28702", "GSE72970", "GSE69657", "GSE39582")
N   <- c(GSE104645 = 113, GSE28702 = 83, GSE72970 = 124, GSE69657 = 30, GSE39582 = 164)
PLAT <- c(GSE104645 = "GPL6480", GSE28702 = "GPL570", GSE72970 = "GPL570",
          GSE69657 = "GPL570", GSE39582 = "GPL570")
W <- unname(sqrt(N[COH]))

zs <- function(p, lfc) sign(lfc) * qnorm(pmax(p / 2, .Machine$double.eps), lower.tail = FALSE)
sto  <- function(P, L) {
  Zm <- do.call(cbind, lapply(seq_along(P), function(i) {
    z <- zs(P[[i]], L[[i]]); z[is.na(z)] <- 0; z }))
  as.numeric(Zm %*% W) / sqrt(sum(W^2))
}
bh <- function(p) {
  # NA is not "highly significant": an absent gene gets p = 1 and stays last in the
  # ordering. An earlier version left its slot at 0, which made every missing gene
  # look like q = 0 and reported the whole table as FDR < 0.10.
  p[is.na(p)] <- 1
  o <- order(p); n <- length(p)
  q <- p[o] * n / seq_len(n)
  out <- numeric(n); out[o] <- pmin(1, rev(cummin(rev(q))))
  out
}

gene_of <- function(cohort) {
  d <- read.csv(paste0(CR, "DEG_", cohort, "_log2corrected.csv"), stringsAsFactors = FALSE)
  m <- read.csv(paste0(IN, "probe2sym_", PLAT[cohort], ".csv"), stringsAsFactors = FALSE)
  d$probe_id <- toupper(trimws(d$probe))
  m$probe_id <- toupper(trimws(m$probe_id))
  d <- merge(d, m[, c("probe_id", "gene_symbol")], by = "probe_id", all.x = TRUE)
  d <- d[!is.na(d$gene_symbol) & d$gene_symbol != "", ]
  d <- d[order(d$P.Value), ]
  d[!duplicated(d$gene_symbol), c("gene_symbol", "logFC", "P.Value")]   # 06's rule
}

per <- lapply(COH, gene_of)
for (i in seq_along(COH))
  say("%-10s genes after collapse %6d   nominal P<0.05 %5d   adj.P<0.05 %4d",
      COH[i], nrow(per[[i]]), sum(per[[i]]$P.Value < .05),
      sum(read.csv(paste0(CR, "DEG_", COH[i], "_log2corrected.csv"))$adj.P.Val < .05))

allg <- sort(unique(unlist(lapply(per, function(x) x$gene_symbol))))
idx  <- lapply(per, function(x) match(allg, x$gene_symbol))
P <- lapply(seq_along(per), function(i) { x <- per[[i]]; x$P.Value[idx[[i]]] })
L <- lapply(seq_along(per), function(i) { x <- per[[i]]; x$logFC[idx[[i]]] })
names(P) <- names(L) <- COH
z <- sto(P, L); p <- 2 * pnorm(abs(z), lower.tail = FALSE); fdr <- bh(p)
n_coh <- rowSums(!do.call(cbind, lapply(P, is.na)))   # cohorts contributing each gene

res <- data.frame(gene_symbol = allg, stouffer_z = z, stouffer_p = p, stouffer_fdr = fdr)
res <- res[order(res$stouffer_p), ]
write.csv(res, paste0(OUT, "deg_meta_corrected.csv"), row.names = FALSE)

say("\ncorrected meta: %d genes tested, %d at BH FDR<0.05, %d at FDR<0.10, %d at nominal p<0.05",
    nrow(res), sum(res$stouffer_fdr < .05), sum(res$stouffer_fdr < .10),
    sum(res$stouffer_p < .05))
say("top 10 by p: %s", paste(sprintf("%s(%.1e)", head(res$gene_symbol, 10), head(res$stouffer_p, 10)), collapse = " "))

# ---- compare with what the manuscript actually reports ------------------------
pub <- read.csv(paste0(IN, "deg_meta_results.csv"), stringsAsFactors = FALSE)
psig <- pub[which(pub$stouffer_fdr < .05), ]
say("\npublished table recomputed on its own stored p-values: %d genes at FDR<0.05 (table's own count of significant: see deg_meta_filtered)",
    nrow(psig))
inter <- intersect(res$gene_symbol[res$stouffer_fdr < .05], psig$gene_symbol)
say("overlap corrected-significant with published-significant: %d of %d corrected / %d published  (Jaccard %.3f)",
    length(inter), sum(res$stouffer_fdr < .05), nrow(psig),
    length(inter) / length(union(res$gene_symbol[res$stouffer_fdr < .05], psig$gene_symbol)))
for (g in c("LDLRAD4", "BAG5", "PCGF5", "PRDM2")) {
  a <- res[res$gene_symbol == g, ]; b <- psig[psig$gene_symbol == g, ]
  say("  %-8s corrected z=%6.2f p=%.2e fdr=%.3f | published z=%6.2f significant=%s",
      g, a$stouffer_z, a$stouffer_p, a$stouffer_fdr,
      if (nrow(b)) b$stouffer_z else NA, if (nrow(b)) "yes" else "no")
}
close(lg)
