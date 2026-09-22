
.libPaths(c("C:/Rlibs", .libPaths()))
suppressPackageStartupMessages(library(survival))
options(warn = -1)
IN <- "C:/xelox_reanalysis/input"
say <- function(...) cat(paste0(...), "\n", sep = "")

cidx <- function(s, t, e) {
  num <- 0; den <- 0
  for (i in which(e == 1)) {
    j <- which(t > t[i]); if (!length(j)) next
    num <- num + sum(s[i] > s[j]) + 0.5 * sum(s[i] == s[j]); den <- den + length(j)
  }
  list(C = if (den) num / den else NA_real_, npair = den)
}
hr_per_sd <- function(s, t, e) {
  if (length(unique(e)) < 2) return(c(NA, NA, NA, NA))
  z <- as.numeric(scale(s))
  f <- tryCatch(coxph(Surv(t, e) ~ z), error = function(x) NULL)
  if (is.null(f)) return(c(NA, NA, NA, NA))
  ci <- tryCatch(exp(confint(f)), error = function(x) c(NA, NA))
  c(exp(coef(f)[1]), ci[1], ci[2], summary(f)$coefficients[1, "Pr(>|z|)"])
}

risk <- read.csv(file.path(IN, "gse72970_risk.csv"), stringsAsFactors = FALSE)
ffx  <- read.csv(file.path(IN, "gse72970_folfox.csv"), stringsAsFactors = FALSE)
cat("risk n =", nrow(risk), "  folfox n =", nrow(ffx), "\n\n")

say("========== A. WHOLE COHORT n=124 ==========")
say(sprintf("PFS events = %d", sum(risk$pfs_event)))
for (nm in c("risk_score", "risk_score_z")) {
  x <- cidx(risk[[nm]], risk$pfs_time, risk$pfs_event)
  h <- hr_per_sd(risk[[nm]], risk$pfs_time, risk$pfs_event)
  say(sprintf("  %-14s  C = %.4f (comparable pairs = %d)   HR/SD = %.3f (%.3f-%.3f)  p = %.3f",
      nm, x$C, x$npair, h[1], h[2], h[3], h[4]))
}
say("")
say("  events by regimen:")
print(table(risk$regimen, risk$pfs_event))
say("")
say("  risk score vs RESPONSE (R = responder 63, NR = non-responder 61):")
r_R  <- risk$risk_score[risk$response_status == "R"]
r_NR <- risk$risk_score[risk$response_status == "NR"]
say(sprintf("    R : n=%d  median %.4f", length(r_R),  median(r_R)))
say(sprintf("    NR: n=%d  median %.4f", length(r_NR), median(r_NR)))
w <- suppressWarnings(wilcox.test(r_R, r_NR))
say(sprintf("    Wilcoxon p = %.4f  -> higher score in: %s",
    w$p.value, if (median(r_R) > median(r_NR)) "R (responder)" else "NR (non-responder)"))
say("")
say("  risk group vs response:")
print(table(risk$risk_group, risk$response_status))
say("")
say("  median risk score by PFS event status:")
print(tapply(risk$risk_score, risk$pfs_event, median))

say("")
say("========== B. FOLFOX SUBGROUP n=32 ==========")
say(sprintf("PFS events = %d", sum(ffx$pfs_event)))
for (nm in c("risk_score_A", "risk_score_B")) {
  x <- cidx(ffx[[nm]], ffx$pfs_time, ffx$pfs_event)
  h <- hr_per_sd(ffx[[nm]], ffx$pfs_time, ffx$pfs_event)
  say(sprintf("  %-14s  C = %.4f (comparable pairs = %d)   HR/SD = %.3f (%.3f-%.3f)  p = %.3f",
      nm, x$C, x$npair, h[1], h[2], h[3], h[4]))
}
m <- match(ffx$sample, risk$sample_id)
ffx$nomogram <- risk$risk_score[m]
say(sprintf("  nomogram score aligned: %d / %d", sum(!is.na(ffx$nomogram)), nrow(ffx)))
x <- cidx(ffx$nomogram, ffx$pfs_time, ffx$pfs_event)
h <- hr_per_sd(ffx$nomogram, ffx$pfs_time, ffx$pfs_event)
say(sprintf("  %-14s  C = %.4f (comparable pairs = %d)   HR/SD = %.3f (%.3f-%.3f)  p = %.3f",
    "nomogram", x$C, x$npair, h[1], h[2], h[3], h[4]))
say("")
say("  pairwise correlations inside FOLFOX subgroup:")
say(sprintf("    cor(nomogram, A) = %.3f    cor(nomogram, B) = %.3f    cor(A, B) = %.3f",
    cor(ffx$nomogram, ffx$risk_score_A), cor(ffx$nomogram, ffx$risk_score_B),
    cor(ffx$risk_score_A, ffx$risk_score_B)))
say("")
say("  median nomogram score by response (0 = responder, 1 = non-responder):")
print(tapply(ffx$nomogram, ffx$response, median))
print(table(ffx$response))
