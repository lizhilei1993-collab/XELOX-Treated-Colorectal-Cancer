IN <- "/path/to/xelox_reanalysis/in/"; CR <- "/path/to/xelox_reanalysis/rerun_nofilt/"
COH <- c("GSE104645","GSE28702","GSE72970","GSE69657","GSE39582")
PLAT <- c(GSE104645="GPL6480", GSE28702="GPL570", GSE72970="GPL570", GSE69657="GPL570", GSE39582="GPL570")
sig <- read.csv("/path/to/xelox_reanalysis/meta_nofilt/deg_meta_corrected.csv")
top <- sig[sig$stouffer_fdr < 0.05, ]
out <- data.frame(gene_symbol = top$gene_symbol, z = top$stouffer_z, p = top$stouffer_p, fdr = top$stouffer_fdr)
for (c in COH) {
  d <- read.csv(paste0(CR, "DEG_", c, "_log2corrected.csv"), stringsAsFactors = FALSE)
  m <- read.csv(paste0(IN, "probe2sym_", PLAT[c], ".csv"), stringsAsFactors = FALSE)
  d$probe_id <- toupper(trimws(d$probe)); m$probe_id <- toupper(trimws(m$probe_id))
  d <- merge(d, m[, c("probe_id","gene_symbol")], by = "probe_id")
  d <- d[order(d$P.Value), ]; d <- d[!duplicated(d$gene_symbol), ]
  i <- match(out$gene_symbol, d$gene_symbol)
  out[[paste0("logfc_", c)]] <- d$logFC[i]; out[[paste0("p_", c)]] <- d$P.Value[i]
  out[[paste0("dir_", c)]]   <- ifelse(is.na(i), NA, ifelse(d$logFC[i] > 0, "up", "down"))
}
pc <- paste0("p_", COH); out$n_observed <- apply(out[, pc], 1, function(r) sum(!is.na(r)))
dc <- paste0("dir_", COH)
out$dir_up   <- apply(out[, dc], 1, function(x) sum(x == "up", na.rm = TRUE))
out$dir_down <- apply(out[, dc], 1, function(x) sum(x == "down", na.rm = TRUE))
out$dir_consensus <- paste0(out$pmax_dir <- pmax(out$dir_up, out$dir_down), "/", out$n_observed)
out$meets_4of5 <- out$pmax_dir >= 4
print(out[, c("gene_symbol","z","p","fdr","dir_up","dir_down","n_observed","dir_consensus","meets_4of5")], row.names = FALSE)
write.csv(out, "/path/to/xelox_reanalysis/meta_nofilt/signature_9_percohort.csv", row.names = FALSE)
cat("\nup:", sum(out$dir_up > out$dir_down), "  down:", sum(out$dir_down > out$dir_up),
    "  meeting >=4/5 directional consensus:", sum(out$meets_4of5), "\n")
