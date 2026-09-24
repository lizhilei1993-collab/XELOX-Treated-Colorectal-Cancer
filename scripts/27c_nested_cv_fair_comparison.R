
.libPaths(c("/path/to/Rlibs", .libPaths()))
suppressPackageStartupMessages({ library(survival); library(glmnet) })
options(warn = -1)
STAGE <- "/path/to/xelox_reanalysis"; IN <- file.path(STAGE, "input"); OUT <- file.path(STAGE, "nested")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ts <- function() format(Sys.time(), "%H:%M:%S")
say <- function(...) cat("[", ts(), "] ", paste0(...), "\n", sep = "")
t0 <- Sys.time()

## ===== PRE-SPECIFIED DESIGN =====
K_OUTER <- 5; N_REPEAT <- 5; K_INNER <- 5
LAM_RULES <- c("lambda.1se", "lambda.min")
GENE_POOL <- 500
K_GRID <- c(2, 3, 5, 7, 10, 15, 20, 30, 50)

ps <- as.matrix(read.csv(file.path(IN, "pathway_scores.csv"), row.names = 1, check.names = FALSE))
storage.mode(ps) <- "double"
cl <- read.csv(file.path(IN, "clinical.csv"), stringsAsFactors = FALSE, check.names = FALSE)
rownames(cl) <- cl$sample_id
ge <- as.matrix(readRDS(file.path(IN, "gene_expr.rds"))); storage.mode(ge) <- "double"

common <- Reduce(intersect, list(colnames(ps), colnames(ge), cl$sample_id))
ok <- !is.na(cl[common, "rfs_event"]) & !is.na(cl[common, "rfs_delay"]) & cl[common, "rfs_delay"] > 0
cohort <- common[ok]
TIME <- as.numeric(cl[cohort, "rfs_delay"]); EVENT <- as.numeric(cl[cohort, "rfs_event"])
X_pw <- t(ps[, cohort, drop = FALSE])
X_gn <- t(ge[, cohort, drop = FALSE])
loc_raw <- cl[cohort, "tumor_location"]
loc <- ifelse(loc_raw == "distal", 1, ifelse(loc_raw == "proximal", 0, NA))
say("cohort n=", length(cohort), " events=", sum(EVENT),
    " | gene ", ncol(X_gn), " | pathway ", ncol(X_pw),
    " | with location ", sum(!is.na(loc)))

cidx <- function(s, t, e) {
  num <- 0; den <- 0
  for (i in which(e == 1)) { j <- which(t > t[i]); if (!length(j)) next
    num <- num + sum(s[i] > s[j]) + 0.5 * sum(s[i] == s[j]); den <- den + length(j) }
  if (den == 0) NA_real_ else num / den
}
score_test_all <- function(X, t, e) {
  n <- nrow(X); ord <- order(t, -e)
  Xo <- X[ord, , drop = FALSE]; eo <- e[ord]
  revcum <- function(v) rev(cumsum(rev(v)))
  cs1 <- apply(Xo, 2, revcum); cs2 <- apply(Xo^2, 2, revcum)
  nrisk <- revcum(rep(1, n)); ev <- which(eo == 1)
  if (!length(ev)) return(rep(NA_real_, ncol(X)))
  s0 <- nrisk[ev]; s1 <- cs1[ev, , drop = FALSE]; s2 <- cs2[ev, , drop = FALSE]
  U <- colSums(Xo[ev, , drop = FALSE] - s1 / s0)
  I <- colSums(s2 / s0 - (s1 / s0)^2); I[I <= 0] <- NA_real_
  abs(U) / sqrt(I)
}
strat_folds <- function(e, K, seed) {
  set.seed(seed); f <- integer(length(e))
  for (v in sort(unique(e))) { idx <- which(e == v)
    f[idx] <- rep(1:K, length.out = length(idx))[sample(length(idx))] }
  f
}
fit_cv <- function(X, t, e, nf) tryCatch(cv.glmnet(X, Surv(t, e), family = "cox", alpha = 1,
  nfolds = nf, standardize = TRUE), error = function(x) NULL)
cal_slope <- function(lp, t, e) {
  f <- tryCatch(coxph(Surv(t, e) ~ lp), error = function(x) NULL)
  if (is.null(f)) NA_real_ else as.numeric(coef(f)[1])
}

## ---- validate fast score test vs coxph ----
set.seed(7); chk <- sample(ncol(X_gn), 10)
st <- score_test_all(X_gn[, chk, drop = FALSE], TIME, EVENT)
say("score-test check: mine vs coxph p-values")
for (i in seq_along(chk)) {
  f <- coxph(Surv(TIME, EVENT) ~ X_gn[, chk[i]])
  say(sprintf("   gene %-10s  mine=%.4f  coxph=%.4f",
      colnames(X_gn)[chk[i]], st[i], summary(f)$coefficients[1, "Pr(>|z|)"]))
}

## ===== SINGLE PASS over (repeat, fold) =====
rec <- list(); cmrec <- list()
for (rp in seq_len(N_REPEAT)) {
  foldid <- strat_folds(EVENT, K_OUTER, 2026 + rp)
  for (k in seq_len(K_OUTER)) {
    tr <- which(foldid != k); te <- which(foldid == k)
    ttr <- TIME[tr]; etr <- EVENT[tr]; tte <- TIME[te]; ete <- EVENT[te]
    ## ---- in-fold rankings (computed ONCE per fold) ----
    sc_g <- score_test_all(X_gn[tr, , drop = FALSE], ttr, etr)
    ord_g <- order(sc_g, decreasing = TRUE)
    sc_p <- score_test_all(X_pw[tr, , drop = FALSE], ttr, etr)
    ord_p <- order(sc_p, decreasing = TRUE)

    ## ---- arm1: gene, in-fold top500 + LASSO-Cox ----
    gsel <- colnames(X_gn)[ord_g[seq_len(min(GENE_POOL, length(ord_g)))]]
    Xtr <- X_gn[tr, gsel, drop = FALSE]; Xte <- X_gn[te, gsel, drop = FALSE]
    cvf <- fit_cv(Xtr, ttr, etr, K_INNER)
    if (!is.null(cvf)) for (rule in LAM_RULES) {
      cf <- as.matrix(coef(cvf, s = rule)); lp_tr <- as.numeric(predict(cvf, newx = Xtr, s = rule, type = "link"))
      lp_te <- as.numeric(predict(cvf, newx = Xte, s = rule, type = "link"))
      rec[[length(rec) + 1]] <- data.frame(Arm = "Gene_500", Repeat = rp, Fold = k, Lambda = rule,
        n_sel = sum(cf[, 1] != 0), C_train = cidx(lp_tr, ttr, etr), C_test = cidx(lp_te, tte, ete),
        cal_slope_test = cal_slope(lp_te, tte, ete), stringsAsFactors = FALSE)
    }

    ## ---- arm2: pathway, all 44 + LASSO-Cox ----
    cvf2 <- fit_cv(X_pw[tr, , drop = FALSE], ttr, etr, K_INNER)
    if (!is.null(cvf2)) for (rule in LAM_RULES) {
      cf <- as.matrix(coef(cvf2, s = rule))
      lp_tr <- as.numeric(predict(cvf2, newx = X_pw[tr, , drop = FALSE], s = rule, type = "link"))
      lp_te <- as.numeric(predict(cvf2, newx = X_pw[te, , drop = FALSE], s = rule, type = "link"))
      rec[[length(rec) + 1]] <- data.frame(Arm = "Pathway_44", Repeat = rp, Fold = k, Lambda = rule,
        n_sel = sum(cf[, 1] != 0), C_train = cidx(lp_tr, ttr, etr), C_test = cidx(lp_te, tte, ete),
        cal_slope_test = cal_slope(lp_te, tte, ete), stringsAsFactors = FALSE)
    }

    ## ---- complexity-matched: plain Cox on in-fold top-k ----
    for (kk in K_GRID) {
      if (kk > length(ord_g)) next
      sf <- colnames(X_gn)[ord_g[seq_len(kk)]]
      dtr <- data.frame(X_gn[tr, sf, drop = FALSE]); dtr$t <- ttr; dtr$e <- etr
      f <- tryCatch(coxph(Surv(t, e) ~ ., data = dtr), error = function(x) NULL)
      if (is.null(f)) next
      dte <- data.frame(X_gn[te, sf, drop = FALSE]); dte$t <- tte; dte$e <- ete
      ltr <- tryCatch(as.numeric(predict(f, newdata = dtr)), error = function(x) NULL)
      lte <- tryCatch(as.numeric(predict(f, newdata = dte)), error = function(x) NULL)
      if (is.null(ltr) || is.null(lte)) next
      cmrec[[length(cmrec) + 1]] <- data.frame(Arm = "Gene", k = kk,
        C_train = cidx(ltr, ttr, etr), C_test = cidx(lte, tte, ete), stringsAsFactors = FALSE)

      sff <- colnames(X_pw)[ord_p[seq_len(min(kk, length(ord_p)))]]
      dtr <- data.frame(X_pw[tr, sff, drop = FALSE]); dtr$t <- ttr; dtr$e <- etr
      f <- tryCatch(coxph(Surv(t, e) ~ ., data = dtr), error = function(x) NULL)
      if (is.null(f)) next
      dte <- data.frame(X_pw[te, sff, drop = FALSE]); dte$t <- tte; dte$e <- ete
      ltr <- tryCatch(as.numeric(predict(f, newdata = dtr)), error = function(x) NULL)
      lte <- tryCatch(as.numeric(predict(f, newdata = dte)), error = function(x) NULL)
      if (is.null(ltr) || is.null(lte)) next
      cmrec[[length(cmrec) + 1]] <- data.frame(Arm = "Pathway", k = kk,
        C_train = cidx(ltr, ttr, etr), C_test = cidx(lte, tte, ete), stringsAsFactors = FALSE)
    }
  }
  say("repeat ", rp, " / ", N_REPEAT, " done")
}

rec <- do.call(rbind, rec)
write.csv(rec, file.path(OUT, "perfold_lasso.csv"), row.names = FALSE)
agg <- do.call(rbind, lapply(split(rec, list(rec$Arm, rec$Lambda)), function(d) data.frame(
  Arm = d$Arm[1], Lambda = d$Lambda[1],
  C_train_mean = round(mean(d$C_train, na.rm = TRUE), 4),
  C_test_mean  = round(mean(d$C_test, na.rm = TRUE), 4),
  C_test_SE    = round(sd(d$C_test, na.rm = TRUE) / sqrt(nrow(d)), 4),
  Nested_optimism = round(mean(d$C_train - d$C_test, na.rm = TRUE), 4),
  C_test_median = round(median(d$C_test, na.rm = TRUE), 4),
  mean_n_sel = round(mean(d$n_sel), 2), pct_null_folds = round(mean(d$n_sel == 0), 3),
  mean_cal_slope = round(mean(d$cal_slope_test, na.rm = TRUE), 4), n_folds = nrow(d))))
agg <- agg[order(agg$Lambda, agg$Arm), ]
write.csv(agg, file.path(OUT, "agg_lasso_nested.csv"), row.names = FALSE)
say("--- nested LASSO results ---"); print(agg, row.names = FALSE)

cmdf <- do.call(rbind, cmrec)
write.csv(cmdf, file.path(OUT, "perfold_complexity.csv"), row.names = FALSE)
cmagg <- do.call(rbind, lapply(split(cmdf, list(cmdf$Arm, cmdf$k)), function(d) data.frame(
  Arm = d$Arm[1], k = d$k[1],
  C_train_mean = round(mean(d$C_train, na.rm = TRUE), 4),
  C_test_mean  = round(mean(d$C_test, na.rm = TRUE), 4),
  C_test_SE    = round(sd(d$C_test, na.rm = TRUE) / sqrt(nrow(d)), 4),
  Nested_optimism = round(mean(d$C_train - d$C_test, na.rm = TRUE), 4), n_folds = nrow(d))))
cmagg <- cmagg[order(cmagg$Arm, cmagg$k), ]
write.csv(cmagg, file.path(OUT, "complexity_matched_nested.csv"), row.names = FALSE)
say("--- complexity-matched nested ---"); print(cmagg, row.names = FALSE)

png(file.path(OUT, "fig_nested_C_vs_k.png"), width = 2400, height = 1600, res = 220)
par(mar = c(4.5, 4.5, 3, 1))
plot(NA, xlim = range(K_GRID), ylim = range(c(cmagg$C_test_mean, 0.5)) + c(-0.03, 0.04), log = "x",
     xlab = "Number of predictors (k, log scale)",
     ylab = "Nested CV C-index (held-out folds)", main = "Complexity-matched nested CV")
abline(h = 0.5, lty = 3, col = "grey40")
for (arm in c("Gene", "Pathway")) {
  d <- cmagg[cmagg$Arm == arm, ]; col <- if (arm == "Gene") "#C0392B" else "#1F6FB4"
  lines(d$k, d$C_test_mean, type = "b", pch = if (arm == "Gene") 19 else 17, lwd = 2.5, col = col)
  arrows(d$k, d$C_test_mean - d$C_test_SE, d$k, d$C_test_mean + d$C_test_SE,
         angle = 90, code = 3, length = 0.03, col = col)
}
legend("topleft", legend = c("Gene level", "Pathway level"), col = c("#C0392B", "#1F6FB4"),
       pch = c(19, 17), lwd = 2.5, bty = "n")
dev.off()
say("TOTAL ", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s")
say("ALL DONE")
