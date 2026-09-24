# ============================================================
# 29_BOOTSTRAP_COEF_DISTRIBUTION.R
# Bootstrap Coefficient Distribution for 8-Variable Nomogram
#
# Purpose:
#   The reviewer correctly noted that "96-100% directional stability"
#   is insufficient — coefficient magnitude stability must also be
#   reported. This script bootstraps the final 8-variable Cox
#   nomogram (B = 1000) and extracts:
#   1. Per-variable coefficient distribution (median, 2.5%, 97.5%)
#   2. Directional stability (proportion of positive/negative)
#   3. Selection frequency under backward AIC
#   4. Forest plot with full bootstrap distributions
#   5. Coefficient density/violin plots
#
# Model specification (from Nomogram, 79 events, 227 samples):
#   HALLMARK_TGF_BETA_SIGNALING
#   HALLMARK_WNT_BETA_CATENIN_SIGNALING
#   KEGG_ECM_RECEPTOR_INTERACTION
#   KEGG_TGF_BETA_SIGNALING_PATHWAY
#   KEGG_PATHWAYS_IN_CANCER
#   HALLMARK_MYC_TARGETS_V2
#   KEGG_COLORECTAL_CANCER
#   tumor_location (distal vs. proximal)
#
# Input:
#   results/tables/GSE39582_xelox_groups.csv
#   results/tables/pathway_activity/GSE39582_pathway_scores.rds
#   results/tables/nomogram/cox_final_results.csv
#
# Output:
#   results/tables/nomogram/bootstrap_coef_distribution.csv
#   results/figures/nomogram/bootstrap_coef_forest.pdf
#   results/figures/nomogram/bootstrap_coef_violin.pdf
# ============================================================

Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.936")
.libPaths(c("/path/to/Rlibs", .libPaths()))

PROJECT <- "X:/"
RESULTS_TAB_DIR <- file.path(PROJECT, "results", "tables")
RESULTS_FIG_DIR <- file.path(PROJECT, "results", "figures")
PA_DIR          <- file.path(RESULTS_TAB_DIR, "pathway_activity")
NOMO_TAB_DIR    <- "/tmp/output/tables/nomogram"
NOMO_FIG_DIR    <- "/tmp/output/figures/nomogram"

set.seed(42)

cat("============================================================\n")
cat("Script 29: Bootstrap Coefficient Distribution (Nomogram)\n")
cat("============================================================\n\n")

# ============================================================
# 0. Load packages
# ============================================================
suppressPackageStartupMessages({
  library(survival)
  library(rms)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

cat("=== [0] Packages loaded ===\n\n")


# ============================================================
# 1. Load data
# ============================================================
cat("=== [1] Loading data ===\n\n")

# Clinical data
clin <- read.csv(file.path(RESULTS_TAB_DIR, "GSE39582_xelox_groups.csv"),
                 stringsAsFactors = FALSE)
rownames(clin) <- clin$sample_id

# Pathway scores
pw_scores <- readRDS(file.path(PA_DIR, "GSE39582_pathway_scores.rds"))

cat(sprintf("  Clinical: %d samples\n", nrow(clin)))
cat(sprintf("  Pathway scores: %d pathways x %d samples\n", nrow(pw_scores), ncol(pw_scores)))

# Prepare survival data
common_samples <- intersect(clin$sample_id, colnames(pw_scores))
cox_df <- clin[common_samples, c("sample_id", "rfs_event", "rfs_delay", "tumor_location")]
cox_df$rfs_event <- as.numeric(cox_df$rfs_event)
cox_df$rfs_delay <- as.numeric(cox_df$rfs_delay)
cox_df <- cox_df[!is.na(cox_df$rfs_event) & !is.na(cox_df$rfs_delay) & cox_df$rfs_delay > 0, ]

# Pathway score matrix
pw_mat <- t(pw_scores[, cox_df$sample_id])
rownames(pw_mat) <- cox_df$sample_id

# Recode tumor_location to binary if needed
if (is.character(cox_df$tumor_location)) {
  cox_df$tumor_distal <- ifelse(grepl("distal|left|descending|sigmoid|rect", 
                                       cox_df$tumor_location, ignore.case = TRUE), 1, 0)
} else {
  cox_df$tumor_distal <- as.numeric(cox_df$tumor_location)
}

surv_obj <- Surv(cox_df$rfs_delay, cox_df$rfs_event)

n_events <- sum(cox_df$rfs_event)
n_total  <- nrow(cox_df)
cat(sprintf("  Valid: n=%d, events=%d (%.1f%%)\n", n_total, n_events, 100*n_events/n_total))


# ============================================================
# 2. Define the final 8-variable model
# ============================================================
cat("\n=== [2] Defining 8-variable nomogram model ===\n\n")

final_vars <- c(
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "KEGG_ECM_RECEPTOR_INTERACTION",
  "KEGG_TGF_BETA_SIGNALING_PATHWAY",
  "KEGG_PATHWAYS_IN_CANCER",
  "HALLMARK_MYC_TARGETS_V2",
  "KEGG_COLORECTAL_CANCER"
)

# Build design matrix
X <- cbind(pw_mat[, final_vars], tumor_distal = cox_df$tumor_distal)
colnames(X) <- c(final_vars, "Tumor location (distal vs. proximal)")

cat(sprintf("  Variables: %d pathways + 1 clinical\n", length(final_vars)))

# Fit the original model
fit_orig <- coxph(surv_obj ~ ., data = as.data.frame(X))
cat(sprintf("  Original C-index = %.4f\n", summary(fit_orig)$concordance["C"]))


# ============================================================
# 3. Bootstrap coefficient distribution
# ============================================================
cat("\n=== [3] Bootstrap coefficient distribution (B = 1000) ===\n\n")

B <- 1000
var_names <- colnames(X)
n_vars <- length(var_names)

boot_coefs <- matrix(NA, nrow = B, ncol = n_vars)
colnames(boot_coefs) <- var_names

for (b in seq_len(B)) {
  # Resample with replacement
  boot_idx <- sample(n_total, replace = TRUE)
  X_boot <- X[boot_idx, ]
  surv_boot <- surv_obj[boot_idx]

  # Fit Cox model
  fit_boot <- tryCatch({
    coxph(surv_boot ~ ., data = as.data.frame(X_boot))
  }, error = function(e) NULL)

  if (!is.null(fit_boot)) {
    coefs <- coef(fit_boot)
    for (j in seq_len(n_vars)) {
      if (var_names[j] %in% names(coefs)) {
        boot_coefs[b, j] <- coefs[var_names[j]]
      }
    }
  }

  if (b %% 200 == 0) cat(sprintf("  Bootstrap %d/%d...\n", b, B))
}

# Summarize
coef_summary <- data.frame(
  Variable = var_names,
  Original_Coef = coef(fit_orig),
  Bootstrap_Median = apply(boot_coefs, 2, median, na.rm = TRUE),
  Bootstrap_Q025 = apply(boot_coefs, 2, quantile, 0.025, na.rm = TRUE),
  Bootstrap_Q975 = apply(boot_coefs, 2, quantile, 0.975, na.rm = TRUE),
  Bootstrap_SD = apply(boot_coefs, 2, sd, na.rm = TRUE),
  Prop_Positive = colMeans(boot_coefs > 0, na.rm = TRUE),
  Prop_Negative = colMeans(boot_coefs < 0, na.rm = TRUE),
  Prop_Zero = colMeans(is.na(boot_coefs) | boot_coefs == 0),
  stringsAsFactors = FALSE
)

# Directional stability
coef_summary$Directional_Stability <- pmax(coef_summary$Prop_Positive, 
                                             coef_summary$Prop_Negative) * 100

# Express HR in exponentiated form
coef_summary$Original_HR <- exp(coef_summary$Original_Coef)
coef_summary$HR_Median <- exp(coef_summary$Bootstrap_Median)
coef_summary$HR_Q025 <- exp(coef_summary$Bootstrap_Q025)
coef_summary$HR_Q975 <- exp(coef_summary$Bootstrap_Q975)

cat("\n  Bootstrap Coefficient Distribution Summary:\n")
cat(sprintf("  %-45s %8s %8s %8s %8s %7s\n",
            "Variable", "Orig", "Median", "Q025", "Q975", "Dir%"))
cat(paste(rep("-", 85), collapse = ""), "\n")
for (i in seq_len(n_vars)) {
  cat(sprintf("  %-45s %8.3f %8.3f %8.3f %8.3f %6.1f%%\n",
              coef_summary$Variable[i],
              coef_summary$Original_Coef[i],
              coef_summary$Bootstrap_Median[i],
              coef_summary$Bootstrap_Q025[i],
              coef_summary$Bootstrap_Q975[i],
              coef_summary$Directional_Stability[i]))
}

write.csv(coef_summary, file.path(NOMO_TAB_DIR, "bootstrap_coef_distribution.csv"), row.names = FALSE)

# Also save as Supplementary Table S14 format
s14 <- coef_summary[, c("Variable", "Original_HR", "HR_Median", "HR_Q025",
                         "HR_Q975", "Directional_Stability")]
names(s14) <- c("Variable", "Original_HR", "Bootstrap_Median_HR",
                "Bootstrap_Q025_HR", "Bootstrap_Q975_HR", "Directional_Stability_Pct")
write.csv(s14, file.path(NOMO_TAB_DIR, "Supplementary_Table_S14.csv"), row.names = FALSE)


# ============================================================
# 4. Forest plot with bootstrap distributions
# ============================================================
cat("\n=== [4] Generating publication-quality forest plot ===\n\n")

# Prepare plot data (use HR scale for clinical interpretability)
plot_df <- coef_summary
plot_df$Variable <- factor(plot_df$Variable, 
                            levels = plot_df$Variable[order(plot_df$HR_Median, decreasing = TRUE)])

pdf(file.path(NOMO_FIG_DIR, "bootstrap_coef_forest.pdf"), width = 9, height = 5)

p1 <- ggplot(plot_df, aes(x = HR_Median, y = Variable)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  # 95% CI from bootstrap
  geom_errorbarh(aes(xmin = HR_Q025, xmax = HR_Q975), height = 0.2, 
                 color = "grey40", linewidth = 0.8) +
  # IQR box
  geom_errorbarh(aes(xmin = exp(Bootstrap_Median - 0.674 * Bootstrap_SD),
                      xmax = exp(Bootstrap_Median + 0.674 * Bootstrap_SD)),
                  height = 0.15, color = "steelblue", linewidth = 2) +
  # Bootstrap median point
  geom_point(aes(fill = ifelse(HR_Median > 1, "#E64B35", "#4DBBD5")),
             size = 4, shape = 21, color = "black", stroke = 0.5) +
  # Original estimate (diamond)
  geom_point(aes(x = Original_HR), shape = 18, size = 3.5, color = "black") +
  scale_fill_identity() +
  scale_x_log10(
    breaks = c(0.4, 0.5, 0.67, 0.8, 1.0, 1.25, 1.5, 2.0),
    labels = c("0.40", "0.50", "0.67", "0.80", "1.00", "1.25", "1.50", "2.00")
  ) +
  labs(
    title = "Bootstrap Coefficient Stability: 8-Variable Cox Nomogram",
    subtitle = sprintf(
      "GSE39582 XELOX cohort (n=%d, %d RFS events). B=1000 bootstrap resamples.\nPoints: bootstrap median | Box: IQR | Whiskers: 2.5th–97.5th %%ile | Diamond: original estimate",
      n_total, n_events
    ),
    x = "Hazard Ratio (log scale, 95% CI from bootstrap)",
    y = ""
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey30"),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10)
  )

# Add directional stability labels
p1 <- p1 + annotate("text",
  x = max(plot_df$HR_Q975, na.rm = TRUE) * 1.1,
  y = seq_len(n_vars),
  label = sprintf("%.1f%%", plot_df$Directional_Stability),
  size = 3, hjust = 0, color = "grey40"
)

print(p1)
dev.off()
cat(sprintf("  -> %s\n", file.path(NOMO_FIG_DIR, "bootstrap_coef_forest.pdf")))


# ============================================================
# 5. Coefficient density/violin plots
# ============================================================
cat("=== [5] Generating coefficient density plots ===\n\n")

# Convert bootstrap coefficients to long format
boot_long <- as.data.frame(boot_coefs)
boot_long$Bootstrap <- seq_len(B)
boot_long <- tidyr::pivot_longer(boot_long, -Bootstrap,
                                  names_to = "Variable", values_to = "Coefficient")
boot_long <- boot_long[!is.na(boot_long$Coefficient), ]

# Add original coefficient for reference
boot_long$Orig_Coef <- coef_summary$Original_Coef[match(boot_long$Variable, coef_summary$Variable)]

pdf(file.path(NOMO_FIG_DIR, "bootstrap_coef_violin.pdf"), width = 10, height = 6)

p2 <- ggplot(boot_long, aes(x = Variable, y = Coefficient)) +
  geom_violin(aes(fill = ifelse(median_cluster > 0, "#E64B3540", "#4DBBD540")),
              alpha = 0.5, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(y = Orig_Coef), shape = 18, size = 3, color = "black") +
  labs(
    title = "Bootstrap Coefficient Distribution: 8-Variable Nomogram (B=1000)",
    subtitle = sprintf("GSE39582 (n=%d, %d events). Triangle: original estimate. Lines: bootstrap median + IQR.",
                       n_total, n_events),
    x = "", y = "Log Hazard Ratio (coefficient)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "none"
  )

print(p2)
dev.off()
cat(sprintf("  -> %s\n", file.path(NOMO_FIG_DIR, "bootstrap_coef_violin.pdf")))


# ============================================================
# 6. Selection frequency under backward AIC
# ============================================================
cat("\n=== [6] Backward AIC selection frequency ===\n\n")

# In each bootstrap, run univariate Cox → backward AIC
# on ALL 44 pathways + tumor_location (replicating the original pipeline)
all_pathways <- colnames(pw_mat)
n_pathways <- length(all_pathways)

sel_freq <- setNames(rep(0, n_pathways + 1), c(all_pathways, "tumor_distal"))

for (b in seq_len(B)) {
  boot_idx <- sample(n_total, replace = TRUE)
  X_all <- cbind(pw_mat[boot_idx, all_pathways], 
                 tumor_distal = cox_df$tumor_distal[boot_idx])
  surv_boot <- surv_obj[boot_idx]

  # Step 1: Univariate Cox screening (p < 0.05)
  sig_vars <- character(0)
  for (v in c(all_pathways, "tumor_distal")) {
    if (sd(X_all[, v], na.rm = TRUE) < 1e-10) next
    fit_uni <- tryCatch({
      coxph(surv_boot ~ X_all[, v])
    }, error = function(e) NULL)
    if (!is.null(fit_uni)) {
      p_val <- summary(fit_uni)$coefficients[1, "Pr(>|z|)"]
      if (!is.na(p_val) && p_val < 0.05) sig_vars <- c(sig_vars, v)
    }
  }

  # Step 2: Backward AIC on significant variables
  if (length(sig_vars) >= 2) {
    X_sig <- X_all[, sig_vars, drop = FALSE]
    full_fit <- tryCatch({
      coxph(surv_boot ~ ., data = as.data.frame(X_sig))
    }, error = function(e) NULL)

    if (!is.null(full_fit)) {
      step_fit <- tryCatch({
        step(full_fit, direction = "backward", trace = 0)
      }, error = function(e) NULL)

      if (!is.null(step_fit)) {
        final_vars_sel <- names(coef(step_fit))
        for (v in final_vars_sel) {
          if (v %in% names(sel_freq)) sel_freq[v] <- sel_freq[v] + 1
        }
      }
    }
  } else if (length(sig_vars) == 1) {
    sel_freq[sig_vars] <- sel_freq[sig_vars] + 1
  }

  if (b %% 200 == 0) cat(sprintf("  AIC selection bootstrap %d/%d...\n", b, B))
}

# Selection frequency table
sel_freq_df <- data.frame(
  Variable = names(sel_freq),
  Selection_Frequency = sel_freq,
  Selection_Rate = 100 * sel_freq / B,
  stringsAsFactors = FALSE
)
sel_freq_df <- sel_freq_df[order(sel_freq_df$Selection_Frequency, decreasing = TRUE), ]
sel_freq_df <- sel_freq_df[sel_freq_df$Selection_Frequency > 0, ]

cat(sprintf("\n  %d / %d variables ever selected\n", nrow(sel_freq_df), length(sel_freq)))
cat("  Top 12 by selection frequency:\n")
head(sel_freq_df, 12) %>% 
  mutate(Display = sprintf("  %-45s %4d/%d (%.1f%%)",
                           Variable, Selection_Frequency, B, Selection_Rate)) %>%
  pull(Display) %>%
  cat(sep = "\n")

write.csv(sel_freq_df, file.path(NOMO_TAB_DIR, "bootstrap_selection_freq.csv"), row.names = FALSE)


# ============================================================
# 7. Selection frequency bar plot (top 15 + final 8 highlighted)
# ============================================================
cat("\n=== [7] Selection frequency visualization ===\n\n")

plot_sel <- sel_freq_df[1:min(15, nrow(sel_freq_df)), ]
plot_sel$Variable <- factor(plot_sel$Variable, 
                             levels = rev(plot_sel$Variable))
plot_sel$InFinal <- plot_sel$Variable %in% 
  c("HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
    "KEGG_ECM_RECEPTOR_INTERACTION", "KEGG_TGF_BETA_SIGNALING_PATHWAY",
    "KEGG_PATHWAYS_IN_CANCER", "HALLMARK_MYC_TARGETS_V2",
    "KEGG_COLORECTAL_CANCER", "tumor_distal")

pdf(file.path(NOMO_FIG_DIR, "bootstrap_selection_freq.pdf"), width = 9, height = 5)

p3 <- ggplot(plot_sel, aes(x = Selection_Rate, y = Variable, fill = InFinal)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "grey70"),
                    labels = c("TRUE" = "In Final Model", "FALSE" = "Not Selected"),
                    name = "") +
  labs(
    title = "Bootstrap Selection Frequency: Backward AIC Cox (B=1000)",
    subtitle = "Top 15 variables by selection rate. Red = in final 8-variable nomogram.",
    x = "Selection Rate (%)", y = ""
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

print(p3)
dev.off()
cat(sprintf("  -> %s\n", file.path(NOMO_FIG_DIR, "bootstrap_selection_freq.pdf")))


# ============================================================
# 8. Summary
# ============================================================
cat("\n============================================================\n")
cat("BOOTSTRAP COEFFICIENT DISTRIBUTION --- SUMMARY\n")
cat("============================================================\n\n")

cat(sprintf("  B = %d resamples, n = %d, events = %d\n", B, n_total, n_events))
cat(sprintf("  Original C-index = %.4f\n", summary(fit_orig)$concordance["C"]))
cat("\n  Directional stability summary:\n")
for (i in seq_len(n_vars)) {
  cat(sprintf("    %-45s %.1f%% (%d/%d positive)\n",
              var_names[i],
              coef_summary$Directional_Stability[i],
              sum(boot_coefs[, i] > 0, na.rm = TRUE), B))
}

cat(sprintf("\n  Mean C-index across bootstraps: %.4f\n",
            mean(apply(boot_coefs, 1, function(bc) {
              if (all(is.na(bc))) return(NA)
              lp <- as.matrix(X) %*% bc
              tryCatch(as.numeric(1 - rcorr.cens(lp, surv_obj)["C Index"]),
                       error = function(e) NA)
            }), na.rm = TRUE)))

cat("\n  FILES GENERATED:\n")
cat(sprintf("    %s/bootstrap_coef_distribution.csv\n", NOMO_TAB_DIR))
cat(sprintf("    %s/Supplementary_Table_S14.csv\n", NOMO_TAB_DIR))
cat(sprintf("    %s/bootstrap_selection_freq.csv\n", NOMO_TAB_DIR))
cat(sprintf("    %s/bootstrap_coef_forest.pdf\n", NOMO_FIG_DIR))
cat(sprintf("    %s/bootstrap_coef_violin.pdf\n", NOMO_FIG_DIR))
cat(sprintf("    %s/bootstrap_selection_freq.pdf\n", NOMO_FIG_DIR))

cat("\n=== Script 29 Complete ===\n")
