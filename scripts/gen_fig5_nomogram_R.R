################################################################################
#  Figure 5 — Nomogram + Calibration Curves (R / rms)
#  Panel A: Standard clinical nomogram (Points → Variables → Total → RFS)
#  Panel B: 12/36/60-month calibration curves
#
#  Dependencies: rms, survival, ggplot2, gridExtra, Hmisc
#  Install if needed:
#    install.packages(c("rms", "survival", "ggplot2", "gridExtra", "Hmisc"))
#
#  ★ HOW TO USE WITH YOUR REAL DATA ★
#    1. Replace the "DUMMY DATA" section with:
#         df <- read.csv("your_data.csv")
#       Ensure columns: time, status, and all 8 predictor columns
#    2. Adjust variable names in the rms::datadist() and cph() calls
#    3. Run the entire script
################################################################################

.libPaths(c("e:/基于可解释性机器学习的XELOX耐药分子指纹研究/r_libs", .libPaths()))
library(rms)
library(survival)
library(ggplot2)
library(gridExtra)

set.seed(42)  # reproducibility

# ============================================================================
#  DUMMY DATA (replace this entire section with your real data)
# ============================================================================
n <- 227  # match your cohort size

# Pathway scores (continuous, range ~0–1, matching ssGSEA-like distributions)
df <- data.frame(
  TGFb_HALLMARK       = rnorm(n, mean = 0.52, sd = 0.18),
  Wnt_HALLMARK        = rnorm(n, mean = 0.48, sd = 0.15),
  ECM_KEGG            = rnorm(n, mean = 0.50, sd = 0.17),
  Colorectal_KEGG     = rnorm(n, mean = 0.45, sd = 0.16),
  MYC_HALLMARK        = rnorm(n, mean = 0.55, sd = 0.14),
  Pathways_KEGG       = rnorm(n, mean = 0.50, sd = 0.13),
  TGFb_KEGG           = rbinom(n, 1, 0.35),          # binary pathway flag
  Location_Distal     = rbinom(n, 1, 0.42)            # 0=proximal, 1=distal
)

# Generate survival outcome matching real coefficients
# log(HR) from your Cox model:
lp_true <- with(df,
  0.627  * TGFb_HALLMARK     +   # HR = 1.872
  0.376  * Wnt_HALLMARK      +   # HR = 1.457
  0.502  * ECM_KEGG          +   # HR = 1.652
  (-0.419) * Colorectal_KEGG +   # HR = 1.520 (positive direction in Cox)
  (-0.384) * MYC_HALLMARK    +   # HR = 0.681
  (-0.488) * Pathways_KEGG   +   # HR = 0.614
  (-0.709) * TGFb_KEGG       +   # HR = 0.492
  0.763  * Location_Distal        # HR = 2.143
)

# Simulate Weibull survival times
shape <- 1.2
scale <- 60 * exp(-0.3 * lp_true)
df$time   <- rweibull(n, shape = shape, scale = scale)
df$status <- rbinom(n, 1, prob = plogis(0.5 * lp_true))  # ~35% event rate
df$time   <- pmin(df$time, 120)  # censor at 120 months
df$time   <- pmax(df$time, 1)    # minimum 1 month

cat("Dummy data summary:\n")
cat(sprintf("  n = %d, events = %d (%.1f%%)\n",
            nrow(df), sum(df$status), 100 * mean(df$status)))
cat(sprintf("  Median follow-up: %.1f months\n", median(df$time)))

# ============================================================================
#  ★★★ REPLACE ABOVE WITH YOUR REAL DATA ★★★
#  Example:
#    df <- read.csv("path/to/your/patient_data.csv")
#    # Ensure columns: time, status, TGFb_HALLMARK, Wnt_HALLMARK, ECM_KEGG,
#    #                 Colorectal_KEGG, MYC_HALLMARK, Pathways_KEGG,
#    #                 TGFb_KEGG, Location_Distal
# ============================================================================

# ============================================================================
#  FIT Cox MODEL (rms)
# ============================================================================
dd <- datadist(df)
options(datadist = "dd")

# Fit the Cox PH model
cox_fit <- cph(
  Surv(time, status) ~ TGFb_HALLMARK + Wnt_HALLMARK + ECM_KEGG +
                        Colorectal_KEGG + MYC_HALLMARK + Pathways_KEGG +
                        TGFb_KEGG + Location_Distal,
  data   = df,
  x = TRUE, y = TRUE,          # required for calibrate()
  surv = TRUE,                  # required for survival predictions
  time.inc = 60                 # primary time horizon
)

cat("\n=== Cox Model Summary ===\n")
print(cox_fit)

# ============================================================================
#  PANEL A: NOMOGRAM
# ============================================================================
# Generate the nomogram object
surv_f <- Survival(cox_fit)
nom <- nomogram(
  cox_fit,
  fun     = list(
    function(x) surv_f(12, x),
    function(x) surv_f(36, x),
    function(x) surv_f(60, x)
  ),
  fun.at  = seq(0.1, 0.9, by = 0.1),
  funlabel = c("12-Month RFS", "36-Month RFS", "60-Month RFS"),
  lp      = TRUE,                # show linear predictor axis
  conf.int = FALSE,
  abbrev  = FALSE                # show full variable names
)

# ============================================================================
#  PANEL B: CALIBRATION CURVES (ggplot2 version for publication quality)
# ============================================================================
# Use rms::calibrate for bootstrap-corrected calibration
# B = 200 bootstrap resamples (increase to 500 for final publication)
cal_12 <- calibrate(cox_fit, u = 12, cmethod = "KM", m = 40, B = 200)
cal_36 <- calibrate(cox_fit, u = 36, cmethod = "KM", m = 40, B = 200)
cal_60 <- calibrate(cox_fit, u = 60, cmethod = "KM", m = 40, B = 200)

# Convert to data frames for ggplot2
cal_to_df <- function(cal_obj, label) {
  data.frame(
    predicted = cal_obj[, "mean.predicted"],
    observed  = cal_obj[, "KM.corrected"],
    se        = cal_obj[, "std.err"],
    timepoint = label
  )
}

df_cal <- rbind(
  cal_to_df(cal_12, "12 months"),
  cal_to_df(cal_36, "36 months"),
  cal_to_df(cal_60, "60 months")
)
df_cal <- na.omit(df_cal)

# Color scheme
cal_colors <- c("12 months" = "#3498DB", "36 months" = "#E74C3C", "60 months" = "#27AE60")

p_cal <- ggplot(df_cal, aes(x = predicted, y = observed, color = timepoint)) +
  # Ideal line
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#AAAAAA", linewidth = 0.6) +
  # Calibration curves with error bars
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5, shape = 16) +
  geom_errorbar(aes(ymin = pmax(observed - se, 0), ymax = pmin(observed + se, 1)),
                width = 0.02, linewidth = 0.5, alpha = 0.6) +
  # Scales
  scale_color_manual(values = cal_colors, name = "Time Point") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  # Labels
  labs(
    x = "Predicted RFS Probability",
    y = "Observed RFS Probability",
    title = "Bootstrap-Corrected Calibration Curves"
  ) +
  # Theme
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(
    plot.title       = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(b = 10)),
    axis.title       = element_text(size = 10),
    axis.text        = element_text(size = 9),
    legend.position  = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = "#CCCCCC", linewidth = 0.3),
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(10, 15, 10, 10)
  ) +
  coord_fixed(ratio = 1)

# ============================================================================
#  COMBINE AND EXPORT
# ============================================================================
# Output path — adjust as needed
out_dir  <- file.path("results", "figures_png", "v5.11")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Export Panel A: Nomogram (large standalone PNG) ---
nom_file <- file.path(out_dir, "Fig5A_nomogram_R.png")
png(nom_file, width = 10, height = 14, units = "in", res = 300, family = "Arial")
par(mar = c(2, 1, 3, 1))  # bottom, left, top, right
plot(nom, xfrac = 0.25,     # fraction of plot for variable names
     cex.var   = 0.9,       # variable label size
     cex.axis  = 0.7,       # axis label size
     lmgp      = 0.3)       # label margin
title("Nomogram for Predicting 12/36/60-Month RFS", cex.main = 1.3, line = 1.5)
mtext("LASSO-Cox (lambda.1se)  |  7-Pathway Signature + Tumor Location  |  C-index = 0.676",
      side = 3, line = 0.3, cex = 0.75, col = "#888888")
dev.off()
cat(sprintf("\nPanel A saved: %s\n", nom_file))

# --- Export Panel B: Calibration (ggplot) ---
cal_file <- file.path(out_dir, "Fig5B_calibration_R.png")
ggsave(cal_file, p_cal, width = 6, height = 5.5, dpi = 300)
cat(sprintf("Panel B saved: %s\n", cal_file))

# --- Export combined figure (A left, B right) ---
combo_file <- file.path(out_dir, "Fig5_nomogram_calibration_R.png")
png(combo_file, width = 16, height = 10, units = "in", res = 300, family = "Arial")

# Set up side-by-side layout
layout_matrix <- matrix(c(1, 2), nrow = 1)
layout(layout_matrix, widths = c(1.4, 1))  # nomogram wider than calibration

# Panel A
par(mar = c(2, 1, 4, 1))
plot(nom, xfrac = 0.28, cex.var = 0.85, cex.axis = 0.65, lmgp = 0.3)
title("A  Nomogram for Predicting RFS", cex.main = 1.2, adj = 0, line = 2.5,
      font.main = 2)
mtext("LASSO-Cox  |  7 Pathways + Location  |  C-index = 0.676",
      side = 3, line = 0.8, adj = 0, cex = 0.65, col = "#888888")

# Panel B
par(mar = c(5, 5, 4, 2))
print(p_cal + ggtitle("B  Calibration Curves") +
        theme(plot.title = element_text(face = "bold", size = 12, hjust = 0)))

dev.off()
cat(sprintf("Combined figure saved: %s\n", combo_file))

# ============================================================================
#  MODEL PERFORMANCE SUMMARY
# ============================================================================
cat("\n=== Model Performance ===\n")
cat(sprintf("  C-index (apparent): %.3f\n", cox_fit$stats["C"]))
cat(sprintf("  n = %d, events = %d\n", nrow(df), sum(df$status)))

# Print calibration slopes
cat("\n=== Calibration (Brier scores) ===\n")
cat(sprintf("  12-month: see cal_12 object\n"))
cat(sprintf("  36-month: see cal_36 object\n"))
cat(sprintf("  60-month: see cal_60 object\n"))

cat("\nDone! All figures saved to:", out_dir, "\n")
