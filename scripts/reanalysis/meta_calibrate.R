# Which weighting reproduces the published stouffer_z column?
#
# Nothing here is assumed. The published meta-analysis table stores every input
# the statistic needs (per-cohort p-values and signed logFC) alongside its own
# stouffer_z, so the arithmetic can be recovered by trying candidate weightings
# and keeping the one that reproduces the recorded numbers. Only then can the
# corrected five-cohort tables be put through the same statistic and compared
# with the published 253.

# staged to an ASCII path: R on this machine cannot open the account path
IN_PUB <- "/path/to/xelox_reanalysis/in/"
OUT    <- "/path/to/xelox_reanalysis/rerun_final/"
lg <- file(paste0(OUT, "meta_calibrate_log.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { l <- sprintf(...); writeLines(l, con = lg); cat(l, "\n") }

COH <- c("GSE104645", "GSE28702", "GSE72970", "GSE69657", "GSE39582")
# sample counts per cohort: as the corrected run actually used them, and as the
# series total them (the published run is not documented either way)
N_USED  <- c(GSE104645 = 113, GSE28702 = 83, GSE72970 = 124, GSE69657 = 30, GSE39582 = 164)
N_PUB   <- c(GSE104645 = 193, GSE28702 = 83, GSE72970 = 124, GSE69657 = 30, GSE39582 = 585)

# signed one-sided z from a two-sided p: direction from logFC, magnitude qnorm(1-p/2)
zs <- function(p, lfc) sign(lfc) * qnorm(pmax(p / 2, .Machine$double.eps), lower.tail = FALSE)

stouffer_z <- function(P, L, w) {
  # vectorised: one signed-z matrix, then a weighted column sum
  Zm <- do.call(cbind, lapply(seq_along(P), function(i) {
    z <- sign(L[[i]]) * qnorm(pmax(P[[i]] / 2, .Machine$double.eps), lower.tail = FALSE)
    z[is.na(z)] <- 0; z
  }))
  as.numeric(Zm %*% w) / sqrt(sum(w^2))
}

pub <- read.csv(paste0(IN_PUB, "deg_meta_results.csv"), stringsAsFactors = FALSE)
P <- setNames(lapply(COH, function(c) pub[[paste0("pval_", c)]]), COH)
L <- setNames(lapply(COH, function(c) pub[[paste0("logfc_", c)]]), COH)
obs <- pub$stouffer_z
say("published table: %d rows, %d with a finite stouffer_z", nrow(pub), sum(is.finite(obs)))

best <- NULL
for (nm in c("equal", "sqrt_n_used", "sqrt_n_published", "n_used", "n_published")) {
  w <- switch(nm, equal = rep(1, 5), sqrt_n_used = unname(sqrt(N_USED[COH])),
              sqrt_n_published = unname(sqrt(N_PUB[COH])), n_used = unname(N_USED[COH]),
              n_published = unname(N_PUB[COH]))
  s <- stouffer_z(P, L, w)
  keep <- is.finite(obs) & is.finite(s)
  mad <- median(abs(obs[keep] - s[keep])); rr <- cor(obs[keep], s[keep])
  ex  <- sum(abs(obs[keep] - s[keep]) < 0.01)
  say("%-17s median|diff|=%8.4f  r=%.6f  reproduced %d/%d", nm, mad, rr, ex, sum(keep))
  if (is.null(best) || ex > best$ex) best <- list(nm = nm, ex = ex, n = sum(keep), w = w)
}
say("CHOSEN weighting: %s (reproduces %d/%d within 0.01)", best$nm, best$ex, best$n)
saveRDS(list(name = best$nm, w = best$w), paste0(OUT, "stouffer_weighting.rds"))
close(lg)
