# Exports the three sheets that make up the rebuilt Supplementary Table S1.
#
# Everything here is the corrected arm: the two linear-scale cohorts (GSE28702,
# GSE69657) log2(x+1)-repaired, the Methods' IQR probe filter NOT applied (it
# postdates the published tables and is now disclaimed in the Methods), and the
# Stouffer weighting that reproduces the published stouffer_z column exactly
# (sqrt(n), n = 113/83/124/30/164; verified 17,517/17,517 in meta_calibrate.R).
#
# Sheet layout deliberately mirrors the old S1 column names so that any text or
# script referring to logfc_GSE*/pval_GSE*/dir_GSE*/stouffer_* still resolves.
# The old fisher_p/fisher_fdr columns are dropped: the Methods describe a
# Stouffer combination and the Fisher column was never quoted anywhere.

IN  <- "/path/to/xelox_reanalysis/in/"
CR  <- Sys.getenv("META_TABLES", "/path/to/xelox_reanalysis/rerun_nofilt/")
OUT <- Sys.getenv("META_OUT",    "/path/to/xelox_reanalysis/meta_nofilt/")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
lg <- file(paste0(OUT, "meta_export_log.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { l <- sprintf(...); writeLines(l, con = lg); cat(l, "\n") }

COH  <- c("GSE104645", "GSE28702", "GSE72970", "GSE69657", "GSE39582")
N    <- c(GSE104645 = 113, GSE28702 = 83, GSE72970 = 124, GSE69657 = 30, GSE39582 = 164)
PLAT <- c(GSE104645 = "GPL6480", GSE28702 = "GPL570", GSE72970 = "GPL570",
          GSE69657 = "GPL570", GSE39582 = "GPL570")
W <- unname(sqrt(N[COH]))

zs <- function(p, lfc) sign(lfc) * qnorm(pmax(p / 2, .Machine$double.eps), lower.tail = FALSE)
sto <- function(P, L) {
  Zm <- do.call(cbind, lapply(seq_along(P), function(i) {
    z <- zs(P[[i]], L[[i]]); z[is.na(z)] <- 0; z }))
  as.numeric(Zm %*% W) / sqrt(sum(W^2))
}
bh <- function(p) {
  p[is.na(p)] <- 1
  o <- order(p); n <- length(p)
  q <- p[o] * n / seq_len(n)
  out <- numeric(n); out[o] <- pmin(1, rev(cummin(rev(q))))
  out
}

gene_of <- function(cohort) {
  d <- read.csv(paste0(CR, "DEG_", cohort, "_log2corrected.csv"), stringsAsFactors = FALSE)
  m <- read.csv(paste0(IN, "probe2sym_", PLAT[cohort], ".csv"), stringsAsFactors = FALSE)
  d$probe_id <- toupper(trimws(d$probe)); m$probe_id <- toupper(trimws(m$probe_id))
  d <- merge(d, m[, c("probe_id", "gene_symbol")], by = "probe_id", all.x = TRUE)
  d <- d[!is.na(d$gene_symbol) & d$gene_symbol != "", ]
  d <- d[order(d$P.Value), ]
  d[!duplicated(d$gene_symbol), c("gene_symbol", "logFC", "P.Value")]
}

per <- lapply(COH, gene_of)
allg <- sort(unique(unlist(lapply(per, function(x) x$gene_symbol))))
idx  <- lapply(per, function(x) match(allg, x$gene_symbol))
P <- lapply(seq_along(per), function(i) per[[i]]$P.Value[idx[[i]]])
L <- lapply(seq_along(per), function(i) per[[i]]$logFC[idx[[i]]])

z <- sto(P, L); p <- 2 * pnorm(abs(z), lower.tail = FALSE); fdr <- bh(p)
Lm <- do.call(cbind, L); Pm <- do.call(cbind, P)
n_coh <- rowSums(!is.na(Pm))
n_up   <- rowSums(Lm > 0, na.rm = TRUE)
n_down <- rowSums(Lm < 0, na.rm = TRUE)
consensus <- ifelse(n_up > n_down, "up", ifelse(n_down > n_up, "down", "tie"))
cons_frac <- pmax(n_up, n_down) / n_coh
mean_lfc  <- rowMeans(Lm, na.rm = TRUE)

full <- data.frame(gene_symbol = allg, stringsAsFactors = FALSE)
for (i in seq_along(COH)) {
  full[[paste0("logfc_", COH[i])]] <- Lm[, i]
  full[[paste0("pval_", COH[i])]]  <- Pm[, i]
  full[[paste0("dir_",  COH[i])]]  <- ifelse(is.na(Lm[, i]), NA, ifelse(Lm[, i] > 0, "up", "down"))
}
full$n_cohorts           <- n_coh
full$n_up                <- n_up
full$n_down              <- n_down
full$direction_consensus <- consensus
full$consensus_frac      <- cons_frac
full$stouffer_z          <- z
full$stouffer_p          <- p
full$stouffer_fdr        <- fdr
full$mean_logfc          <- mean_lfc
full <- full[order(full$stouffer_p), ]

write.csv(full, paste0(OUT, "s1_meta_all_genes.csv"), row.names = FALSE)
sig <- full[full$stouffer_fdr < 0.05, ]
write.csv(sig, paste0(OUT, "s1_meta_signature.csv"), row.names = FALSE)

dm <- data.frame(gene_symbol = full$gene_symbol, stringsAsFactors = FALSE)
for (i in seq_along(COH)) dm[[paste0("dir_", COH[i])]] <- full[[paste0("dir_", COH[i])]]
dm$n_cohorts <- full$n_cohorts
dm$n_concordant <- pmax(full$n_up, full$n_down)
dm$direction_consensus <- full$direction_consensus
dm$concordance <- sprintf("%d/%d", dm$n_concordant, dm$n_cohorts)
dm <- dm[order(dm$gene_symbol), ]
write.csv(dm, paste0(OUT, "s1_direction_matrix.csv"), row.names = FALSE)

say("genes exported          : %d", nrow(full))
say("signature (FDR<0.05)    : %d  -> %s", nrow(sig), paste(sig$gene_symbol, collapse = ", "))
say("  up %d / down %d ; |z| range %.2f-%.2f",
    sum(sig$direction_consensus == "up"), sum(sig$direction_consensus == "down"),
    min(abs(sig$stouffer_z)), max(abs(sig$stouffer_z)))
say("  concordant in 5/5     : %d ; in 4/5 : %d",
    sum(dm$concordance == "5/5" & dm$gene_symbol %in% sig$gene_symbol),
    sum(dm$concordance == "4/5" & dm$gene_symbol %in% sig$gene_symbol))
say("max |logfc| anywhere    : %.3f  (was >400 in the unrepaired linear-scale cohorts)",
    max(abs(Lm), na.rm = TRUE))
say("rows with any |logfc|>5 : %d  (was 3,895 in the published S1)",
    sum(rowSums(!is.na(Lm) & abs(Lm) > 5) > 0))
say("direction matrix rows   : %d", nrow(dm))
say("wrote s1_meta_all_genes.csv, s1_meta_signature.csv, s1_direction_matrix.csv to %s", OUT)
close(lg)
