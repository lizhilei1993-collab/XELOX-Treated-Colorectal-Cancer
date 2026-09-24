# ============================================================
# 27_FIXED -- Gene vs Pathway LASSO-Cox FAIR COMPARISON (bug-fixed)
#
# WHY THIS FILE EXISTS
#   The original 27_gene_lasso_cox_fair_comparison.R hard-coded
#   s = "lambda.1se" when computing (a) the apparent C-index and
#   (b) the bootstrap predictions, while a silent fallback to
#   lambda.min was applied ONLY to coefficient extraction and the
#   feature count. At lambda.1se both feature levels select ZERO
#   variables, so the reported apparent C-index of 0.500 and the
#   optimism values (0.0945 vs 0.011) came from NULL models and
#   are not interpretable.
#
# WHAT IS FIXED
#   1. A single lambda rule is used CONSISTENTLY for selection,
#      apparent C-index and bootstrap predictions. No fallback.
#   2. Both lambda rules (min and 1se) are run side by side and
#      reported, so the reader can see that lambda.1se is a null
#      model at this sample size.
#   3. Number of selected features and the proportion of NULL
#      bootstrap replicates are reported (this is the diagnostic
#      that exposes the problem).
#   4. EPV = events / n_selected, computed consistently.
#   5. A complexity-matched block (plain Cox on the in-data top-k
#      univariate features, k on a fixed grid) is added, because
#      optimism is a monotone function of model complexity and a
#      cross-level comparison WITHOUT complexity matching is not
#      interpretable.
#   6. Fast vectorised Harrell C-index (risk convention). NOTE:
#      Hmisc::rcorr.cens is ~100x slower and dominated runtime.
#
# FOR THE MANUSCRIPT, PREFER scripts/27c_nested_cv_fair_comparison.R
#   which additionally removes feature-selection leakage by doing
#   all selection INSIDE training folds (nested CV).
#
# HOW TO RUN
#   Must be run with the working directory set to the project root
#   (the folder that contains data/ and results/). This avoids any
#   non-ASCII literal in the source, which R cannot resolve when
#   its locale is "C".
#
#   cd /d <project root>
#   Rscript scripts/27_gene_lasso_cox_fair_comparison_FIXED.R
# ============================================================

# Windows TEMP can sit under a CJK user profile; HDF5Array/GSVA fail there.
Sys.setenv(TMPDIR = "/tmp", TMP = "/tmp", TEMP = "/tmp")
.libPaths(c("/path/to/Rlibs", .libPaths()))
suppressPackageStartupMessages({ library(survival); library(glmnet) })
options(warn = -1)
set.seed(42)

PROJECT_ROOT <- getwd()
if (!dir.exists(file.path(PROJECT_ROOT, "results", "tables"))) {
  stop("Run this script with the working directory set to the project root.")
}
PAC <- file.path(PROJECT_ROOT, "results", "tables", "pathway_activity")
RES <- file.path(PROJECT_ROOT, "results", "tables")
OUT <- file.path(RES, "fair_comparison_CORRECTED")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

B_BOOT <- 300        # bootstrap replicates for optimism
K_INNER <- 10        # folds for cv.glmnet
GENE_POOL <- 500     # candidate genes (top by univariate Cox) -- as in original
K_GRID <- c(2, 3, 5, 7, 10, 15, 20, 30, 44, 50, 74)
LAM_RULES <- c("lambda.min", "lambda.1se")

say <- function(...) cat(paste0(...), "\n", sep = "")

## ---------------- load ----------------
ps <- as.matrix(read.csv(file.path(PAC, "GSE39582_pathway_scores.csv"),
                         row.names = 1, check.names = FALSE))
storage.mode(ps) <- "double"
cl <- read.csv(file.path(RES, "GSE39582_xelox_groups.csv"),
               stringsAsFactors = FALSE, check.names = FALSE)
rownames(cl) <- cl$sample_id
ge <- as.matrix(readRDS(file.path(PAC, "GSE39582_gene_expression.rds")))
storage.mode(ge) <- "double"

gu <- read.csv(file.path(RES, "gene_level", "gene_cox_univariate.csv"),
               stringsAsFactors = FALSE, check.names = FALSE)

## ---------------- cohort ----------------
common <- Reduce(intersect, list(colnames(ps), colnames(ge), cl$sample_id))
ok <- !is.na(cl[common, "rfs_event"]) & !is.na(cl[common, "rfs_delay"]) &
      cl[common, "rfs_delay"] > 0
cohort <- common[ok]
TIME  <- as.numeric(cl[cohort, "rfs_delay"])
EVENT <- as.numeric(cl[cohort, "rfs_event"])
X_pw  <- t(ps[, cohort, drop = FALSE])
X_pw  <- X_pw[, colSums(is.na(X_pw)) == 0, drop = FALSE]

gu_ok <- gu[!is.na(gu$p_value) & gu$Gene %in% rownames(ge), ]
gu_ok <- gu_ok[order(gu_ok$p_value), ]
topk  <- head(unique(gu_ok$Gene), GENE_POOL)
X_gn  <- t(ge[topk, cohort, drop = FALSE])
X_gn  <- X_gn[, colSums(is.na(X_gn)) == 0, drop = FALSE]

say("cohort n = ", length(cohort), "  events = ", sum(EVENT))
say("design: gene ", nrow(X_gn), "x", ncol(X_gn), " | pathway ", nrow(X_pw), "x", ncol(X_pw))

## ---------------- fast Harrell C (risk convention) ----------------
cidx <- function(s, t, e) {
  num <- 0; den <- 0
  for (i in which(e == 1)) {
    j <- which(t > t[i]); if (!length(j)) next
    num <- num + sum(s[i] > s[j]) + 0.5 * sum(s[i] == s[j]); den <- den + length(j)
  }
  if (den == 0) NA_real_ else num / den
}
## validate against an independent implementation
set.seed(1); sv <- rnorm(60); tt <- rexp(60); ee <- rbinom(60, 1, 0.7)
say("C-index check: mine = ", round(cidx(sv, tt, ee), 4),
    "   1 - survival::concordance = ",
    round(1 - as.numeric(survival::concordance(Surv(tt, ee) ~ sv)$concordance), 4))

fit_cv  <- function(X, t, e) tryCatch(cv.glmnet(X, Surv(t, e), family = "cox", alpha = 1,
                        nfolds = K_INNER, standardize = TRUE), error = function(z) NULL)
sel_n   <- function(cv, s) { cf <- as.matrix(coef(cv, s = s)); sum(cf[, 1] != 0) }
lp_at   <- function(cv, X, s) as.numeric(predict(cv, newx = X, s = s, type = "link"))

## ---------------- consistent-lambda fair comparison ----------------
run_level <- function(X, tag, rule) {
  cvf <- fit_cv(X, TIME, EVENT)
  if (is.null(cvf)) { say(tag, " ", rule, ": FULL FIT FAILED"); return(NULL) }
  n_full <- sel_n(cvf, rule)
  c_app  <- cidx(lp_at(cvf, X, rule), TIME, EVENT)      # SAME lambda as the model
  n <- nrow(X); ops <- numeric(0); nsel <- integer(0)
  for (b in seq_len(B_BOOT)) {
    idx <- sample(1:n, n, replace = TRUE)
    Xb <- X[idx, , drop = FALSE]; tb <- TIME[idx]; eb <- EVENT[idx]
    vb <- tryCatch(fit_cv(Xb, tb, eb), error = function(z) NULL)
    if (is.null(vb)) next
    nsel <- c(nsel, sel_n(vb, rule))                    # SAME lambda in bootstrap
    cb <- cidx(lp_at(vb, Xb, rule), tb, eb)
    co <- cidx(lp_at(vb, X,  rule), TIME, EVENT)
    if (!is.na(cb) && !is.na(co)) ops <- c(ops, cb - co)
  }
  o <- mean(ops, na.rm = TRUE)
  data.frame(Level = tag, Lambda = rule, N_features_full = n_full,
    EPV = round(sum(EVENT) / max(n_full, 1), 2),
    C_apparent = round(c_app, 4), Optimism = round(o, 4),
    C_corrected = round(c_app - o, 4),
    Mean_n_selected_boot = round(mean(nsel, na.rm = TRUE), 2),
    Pct_null_boot = round(mean(nsel == 0, na.rm = TRUE), 3),
    B_ok = length(ops), stringsAsFactors = FALSE)
}

acc <- list()
for (rule in LAM_RULES) {
  for (lv in c("Gene", "Pathway")) {
    r <- run_level(if (lv == "Gene") X_gn else X_pw, lv, rule)
    if (!is.null(r)) {
      acc[[length(acc) + 1]] <- r
      say(sprintf("%-8s | %-12s n_sel=%3d  EPV=%6.2f  C_app=%.4f  optimism=%.4f  C_corr=%.4f  boot-null=%4.1f%%",
        r$Level, r$Lambda, r$N_features_full, r$EPV, r$C_apparent,
        r$Optimism, r$C_corrected, 100 * r$Pct_null_boot))
    }
  }
}
res <- do.call(rbind, acc)
write.csv(res, file.path(OUT, "fair_comparison_CORRECTED.csv"), row.names = FALSE)

## ---------------- complexity-matched ----------------
say("")
say("--- complexity-matched (plain Cox, top-k univariate, boot B=", B_BOOT, ") ---")
univ_score <- function(X, t, e) {
  n <- nrow(X); ord <- order(t, -e)
  Xo <- X[ord, , drop = FALSE]; eo <- e[ord]
  rc <- function(v) rev(cumsum(rev(v)))
  cs1 <- apply(Xo, 2, rc); cs2 <- apply(Xo^2, 2, rc)
  nr <- rc(rep(1, n)); ev <- which(eo == 1)
  if (!length(ev)) return(rep(NA_real_, ncol(X)))
  s0 <- nr[ev]; s1 <- cs1[ev, , drop = FALSE]; s2 <- cs2[ev, , drop = FALSE]
  U <- colSums(Xo[ev, , drop = FALSE] - s1 / s0)
  I <- colSums(s2 / s0 - (s1 / s0)^2); I[I <= 0] <- NA_real_
  abs(U) / sqrt(I)
}
ord_pw <- names(sort(univ_score(X_pw, TIME, EVENT), decreasing = TRUE))
ord_gn <- names(sort(univ_score(X_gn, TIME, EVENT), decreasing = TRUE))
opt_plain <- function(X, k) {
  if (k > ncol(X)) return(NULL)
  d0 <- data.frame(X[, colnames(X)[1:k], drop = FALSE]); d0$t <- TIME; d0$e <- EVENT
  f <- tryCatch(coxph(Surv(t, e) ~ ., data = d0), error = function(z) NULL)
  if (is.null(f)) return(NULL)
  c_app <- cidx(as.numeric(predict(f, newdata = d0)), TIME, EVENT)
  n <- nrow(d0); ops <- numeric(0)
  for (b in seq_len(B_BOOT)) {
    idx <- sample(1:n, n, replace = TRUE); db <- d0[idx, , drop = FALSE]
    fb <- tryCatch(coxph(Surv(t, e) ~ ., data = db), error = function(z) NULL)
    if (is.null(fb)) next
    pb <- tryCatch(as.numeric(predict(fb, newdata = db)), error = function(z) NULL)
    po <- tryCatch(as.numeric(predict(fb, newdata = d0)), error = function(z) NULL)
    if (is.null(pb) || is.null(po)) next
    cb <- cidx(pb, db$t, db$e); co <- cidx(po, TIME, EVENT)
    if (!is.na(cb) && !is.na(co)) ops <- c(ops, cb - co)
  }
  data.frame(k = k, C_apparent = round(c_app, 4), Optimism = round(mean(ops, na.rm = TRUE), 4),
             C_corrected = round(c_app - mean(ops, na.rm = TRUE), 4),
             EPV = round(sum(EVENT) / k, 2), B_ok = length(ops))
}
cc <- list()
for (k in K_GRID) {
  r <- opt_plain(X_pw[, ord_pw, drop = FALSE], k)
  if (!is.null(r)) { r$Level <- "Pathway"; cc[[length(cc) + 1]] <- r }
  if (k <= ncol(X_gn)) {
    r <- opt_plain(X_gn[, ord_gn, drop = FALSE], k)
    if (!is.null(r)) { r$Level <- "Gene"; cc[[length(cc) + 1]] <- r }
  }
}
cm <- do.call(rbind, cc)[, c("Level", "k", "C_apparent", "Optimism", "C_corrected", "EPV", "B_ok")]
write.csv(cm, file.path(OUT, "complexity_matched_CORRECTED.csv"), row.names = FALSE)
say(""); print(cm, row.names = FALSE)
say(""); say("outputs written to: ", OUT)
