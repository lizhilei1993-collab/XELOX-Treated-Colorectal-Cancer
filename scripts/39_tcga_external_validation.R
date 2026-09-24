
## TEMP must be ASCII: the Windows user profile path contains CJK characters,
## which breaks HDF5Array::.onLoad (registered in GSVA's namespace).
Sys.setenv(TMPDIR = "/tmp", TMP = "/tmp", TEMP = "/tmp")
.libPaths(c("/path/to/Rlibs", .libPaths()))
suppressPackageStartupMessages({
  library(survival); library(GSVA); library(GSEABase); library(timeROC)
})
options(warn = -1)
IN  <- "/path/to/xelox_reanalysis/input"
OUT <- "/path/to/xelox_reanalysis/tcga_val"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ts <- function() format(Sys.time(), "%H:%M:%S")
say <- function(...) cat("[", ts(), "] ", paste0(...), "\n", sep = "")
t0 <- Sys.time()

## ============ frozen model (from cox_final_model.txt) ============
COEF <- c(
  HALLMARK_TGF_BETA_SIGNALING          =  0.6269,
  HALLMARK_WNT_BETA_CATENIN_SIGNALING  =  0.3765,
  KEGG_ECM_RECEPTOR_INTERACTION        =  0.5019,
  KEGG_TGF_BETA_SIGNALING_PATHWAY      = -0.7095,
  KEGG_PATHWAYS_IN_CANCER              = -0.4886,
  HALLMARK_MYC_TARGETS_V2              = -0.3843,
  KEGG_COLORECTAL_CANCER               =  0.4188
)
LOC_COEF <- 0.7624
say("frozen model: ", length(COEF), " pathways + location_distal")

## ============ 44 gene sets (18 HALLMARK + 26 KEGG) ============
read_gmt <- function(f) {
  ln <- readLines(f, warn = FALSE)
  out <- list()
  for (l in ln) { p <- strsplit(l, "\t")[[1]]; g <- p[-(1:2)]; g <- g[nchar(g) > 0]
    if (length(g)) out[[p[1]]] <- g }
  out
}
H <- read_gmt("/tmp/h.all.v2024.1.symbols.gmt")
K <- read_gmt("/tmp/c2.cp.kegg_legacy.v2024.1.symbols.gmt")
kegg_targets <- c("KEGG_NUCLEOTIDE_EXCISION_REPAIR","KEGG_BASE_EXCISION_REPAIR","KEGG_MISMATCH_REPAIR",
 "KEGG_HOMOLOGOUS_RECOMBINATION","KEGG_P53_SIGNALING_PATHWAY","KEGG_APOPTOSIS","KEGG_CELL_CYCLE",
 "KEGG_ABC_TRANSPORTERS","KEGG_GLUTATHIONE_METABOLISM","KEGG_DRUG_METABOLISM_CYTOCHROME_P450",
 "KEGG_DRUG_METABOLISM_OTHER_ENZYMES","KEGG_METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450",
 "KEGG_WNT_SIGNALING_PATHWAY","KEGG_TGF_BETA_SIGNALING_PATHWAY","KEGG_NOTCH_SIGNALING_PATHWAY",
 "KEGG_FOCAL_ADHESION","KEGG_ECM_RECEPTOR_INTERACTION","KEGG_MAPK_SIGNALING_PATHWAY",
 "KEGG_MTOR_SIGNALING_PATHWAY","KEGG_JAK_STAT_SIGNALING_PATHWAY","KEGG_CHEMOKINE_SIGNALING_PATHWAY",
 "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY","KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
 "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY","KEGG_COLORECTAL_CANCER","KEGG_PATHWAYS_IN_CANCER")
hall_targets <- c("HALLMARK_DNA_REPAIR","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_APOPTOSIS",
 "HALLMARK_ANGIOGENESIS","HALLMARK_HYPOXIA","HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_GLYCOLYSIS",
 "HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_WNT_BETA_CATENIN_SIGNALING",
 "HALLMARK_P53_PATHWAY","HALLMARK_PI3K_AKT_MTOR_SIGNALING","HALLMARK_MYC_TARGETS_V1",
 "HALLMARK_MYC_TARGETS_V2","HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_MTORC1_SIGNALING",
 "HALLMARK_TGF_BETA_SIGNALING")
gs <- list()
for (n in kegg_targets) if (n %in% names(K)) gs[[n]] <- K[[n]]
for (n in hall_targets) if (n %in% names(H)) gs[[n]] <- H[[n]]
say("gene sets loaded: ", length(gs))
missing <- setdiff(names(COEF), names(gs))
if (length(missing)) say("!! MISSING from GMT: ", paste(missing, collapse = ", "))

## ============ expression -> symbols ============
cat("loading expression ...\n")
ex <- readRDS(file.path(IN, "tcga_logcpm.rds"))
gi <- read.csv(file.path(IN, "tcga_gene_info.csv"), stringsAsFactors = FALSE)
cat("expr:", nrow(ex), "x", ncol(ex), "  gene_info:", nrow(gi), "\n")
map <- setNames(gi$gene_name, gi$gene_id)
sym <- map[rownames(ex)]
keep <- !is.na(sym) & sym != ""
ex <- ex[keep, , drop = FALSE]; sym <- sym[keep]
## collapse duplicates by max variance row
ord <- order(-apply(ex, 1, var))
ex <- ex[ord, , drop = FALSE]; sym <- sym[ord]
ex <- ex[!duplicated(sym), , drop = FALSE]; rownames(ex) <- sym[!duplicated(sym)]
cat("after symbol mapping + dedup:", nrow(ex), "genes x", ncol(ex), "samples\n")

need <- unique(unlist(gs))
cat("genes in 44 sets present in matrix:", sum(need %in% rownames(ex)), "/", length(need), "\n")

## ============ ssGSEA (same params as training: minSize=10, maxSize=500) ============
cat("running ssGSEA ...\n")
gsc <- lapply(names(gs), function(n) GeneSet(setName = n, geneIds = gs[[n]]))
gsc <- GeneSetCollection(gsc)
pw <- gsva(ssgseaParam(ex, gsc, minSize = 10, maxSize = 500, verbose = FALSE))
say("ssGSEA done: ", nrow(pw), " x ", ncol(pw))
write.csv(pw, file.path(OUT, "tcga_pathway_scores.csv"))

## ============ clinical ============
cl <- read.csv(file.path(IN, "tcga_clinical_clean.csv"), stringsAsFactors = FALSE)
cl$patient <- substr(cl$barcode, 1, 12)
samples <- colnames(pw); pw_pat <- substr(samples, 1, 12)
m <- match(pw_pat, cl$patient)
say("samples matched to clinical: ", sum(!is.na(m)), "/", length(samples))

PROX <- c("Cecum", "Ascending Colon", "Hepatic Flexure", "Transverse Colon")
DIST <- c("Splenic Flexure", "Descending Colon", "Sigmoid Colon", "Rectosigmoid Junction")
locd <- ifelse(cl$anatomic_subdivision[m] %in% DIST, 1,
        ifelse(cl$anatomic_subdivision[m] %in% PROX, 0, NA))
say("location_distal: distal=", sum(locd == 1, na.rm = TRUE),
    " proximal=", sum(locd == 0, na.rm = TRUE), " unknown=", sum(is.na(locd)))

df <- data.frame(
  sample = samples,
  os_time = as.numeric(cl$os_time[m]),
  os_event = as.numeric(cl$os_event[m]),
  location_distal = locd,
  stringsAsFactors = FALSE
)
for (n in names(COEF)) { if (n %in% rownames(pw)) df[[n]] <- as.numeric(pw[n, ]) else df[[n]] <- 0 }
df$score <- as.numeric(as.matrix(df[, names(COEF)]) %*% COEF + LOC_COEF * ifelse(is.na(df$location_distal), 0, df$location_distal))
df <- df[!is.na(df$os_time) & !is.na(df$os_event) & df$os_time > 0, ]
say("analysis set n = ", nrow(df), "  deaths = ", sum(df$os_event))
write.csv(df, file.path(OUT, "tcga_scored.csv"), row.names = FALSE)

## ============ analysis ============
cidx <- function(s, t, e) {
  num <- 0; den <- 0
  for (i in which(e == 1)) { j <- which(t > t[i]); if (!length(j)) next
    num <- num + sum(s[i] > s[j]) + 0.5 * sum(s[i] == s[j]); den <- den + length(j) }
  list(C = if (den) num / den else NA_real_, npair = den)
}
res <- list()
x <- cidx(df$score, df$os_time, df$os_event)
z <- as.numeric(scale(df$score))
f <- coxph(Surv(os_time, os_event) ~ z, data = df)
ci <- exp(confint(f))
res[[1]] <- data.frame(Analysis = "Frozen score", n = nrow(df), Events = sum(df$os_event),
  C_index = round(x$C, 4), Comparable_pairs = x$npair,
  HR_per_SD = round(exp(coef(f)), 3), CI_low = round(ci[1], 3), CI_high = round(ci[2], 3),
  p = signif(summary(f)$coefficients[1, "Pr(>|z|)"], 3), Calib_slope = NA)

g <- coxph(Surv(os_time, os_event) ~ I(df$score), data = df)
say("")
say("===== EXTERNAL VALIDATION (TCGA-COADREAD, OS) =====")
say("n = ", nrow(df), "   deaths = ", sum(df$os_event))
say(sprintf("  C-index = %.4f  (comparable pairs = %d)", x$C, x$npair))
say(sprintf("  HR per SD = %.3f (%.3f-%.3f)  p = %.3g",
    exp(coef(f)), ci[1], ci[2], summary(f)$coefficients[1, "Pr(>|z|)"]))
say(sprintf("  calibration slope = %.3f (ideal 1.0)", coef(g)[1]))
res[[1]]$Calib_slope <- round(coef(g)[1], 3)

## risk groups (median split) + KM log-rank
grp <- ifelse(df$score > median(df$score), "High", "Low")
sd_ <- survdiff(Surv(os_time, os_event) ~ grp, data = df)
p_lr <- 1 - pchisq(sd_$chisq, length(sd_$n) - 1)
say(sprintf("  median-split log-rank p = %.4g", p_lr))
res[[1]]$Logrank_p <- signif(p_lr, 3)

## time-dependent AUC
tr <- timeROC(T = df$os_time, delta = df$os_event, marker = df$score,
              cause = 1, times = c(365, 1095, 1825), iid = FALSE)
say("  time-dependent AUC: 1y = ", round(tr$AUC[1], 3),
    "  3y = ", round(tr$AUC[2], 3), "  5y = ", round(tr$AUC[3], 3))
res[[1]]$AUC_1y <- round(tr$AUC[1], 3); res[[1]]$AUC_3y <- round(tr$AUC[2], 3); res[[1]]$AUC_5y <- round(tr$AUC[3], 3)

## location as a component: score without location
df$score_noloc <- as.numeric(as.matrix(df[, names(COEF)]) %*% COEF)
x2 <- cidx(df$score_noloc, df$os_time, df$os_event)
say(sprintf("  [control] score WITHOUT tumor location: C = %.4f", x2$C))
res[[2]] <- data.frame(Analysis = "Score without location", n = nrow(df), Events = sum(df$os_event),
  C_index = round(x2$C, 4), Comparable_pairs = x2$npair, HR_per_SD = NA, CI_low = NA,
  CI_high = NA, p = NA, Calib_slope = NA, Logrank_p = NA, AUC_1y = NA, AUC_3y = NA, AUC_5y = NA)
## location alone
x3 <- cidx(ifelse(is.na(df$location_distal), 0, df$location_distal), df$os_time, df$os_event)
say(sprintf("  [control] tumor location alone:          C = %.4f  (n with location = %d)",
    x3$C, sum(!is.na(df$location_distal))))
res[[3]] <- data.frame(Analysis = "Tumor location alone", n = sum(!is.na(df$location_distal)),
  Events = sum(df$os_event[!is.na(df$location_distal)]), C_index = round(x3$C, 4),
  Comparable_pairs = x3$npair, HR_per_SD = NA, CI_low = NA, CI_high = NA, p = NA,
  Calib_slope = NA, Logrank_p = NA, AUC_1y = NA, AUC_3y = NA, AUC_5y = NA)

out <- do.call(rbind, res)
write.csv(out, file.path(OUT, "external_validation_TCGA.csv"), row.names = FALSE)
say("")
print(out, row.names = FALSE)

## KM figure + calibration
png(file.path(OUT, "fig_km_riskgroups.png"), width = 1800, height = 1500, res = 200)
plot(survfit(Surv(os_time, os_event) ~ grp, data = df), col = c("#C0392B", "#1F6FB4"), lwd = 2.5,
     xlab = "Days", ylab = "Overall survival", main = "TCGA-COADREAD: frozen score, median split")
legend("topright", legend = c("High score", "Low score"), col = c("#C0392B", "#1F6FB4"), lwd = 2.5, bty = "n")
dev.off()

## deciles-based calibration (Breslow baseline from the validation set = post-hoc recalibration, declared)
bh <- basehaz(g, centered = FALSE)
df$lp <- coef(g)[1] * df$score
S0 <- function(tt) exp(-approx(bh$time, bh$hazard, xout = tt, rule = 2)$y)
for (tt in c(1095, 1825)) {
  pred <- 1 - S0(tt)^exp(df$lp)
  q <- cut(pred, breaks = quantile(pred, seq(0, 1, 0.2)), include.lowest = TRUE)
  obs <- sapply(split(df, q), function(d) {
    f2 <- survfit(Surv(os_time, os_event) ~ 1, data = d); s <- summary(f2, times = tt)
    1 - s$surv[1] })
  pp <- sapply(split(pred, q), mean)
  say(sprintf("  calibration @ %dy: predicted %s", tt / 365,
      paste(sprintf("%.2f", pp), collapse = " ")))
  say(sprintf("                       observed  %s", paste(sprintf("%.2f", obs), collapse = " ")))
}
say("TOTAL ", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s")
say("ALL DONE")
