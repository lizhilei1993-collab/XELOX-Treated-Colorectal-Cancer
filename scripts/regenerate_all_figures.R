#!/usr/bin/env Rscript
# Regenerate all figures - uses ASCII paths only
.libPaths(c('/tmp/rlib', .libPaths()))
library(ggplot2)

OUT_DIR <- "/tmp/fig_new"
DATA <- "/tmp/scratch/fig_regen"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Inline theme definition (avoids Chinese path source issue)
theme_publication <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
  theme(text = element_text(family = "serif", size = base_size, color = "black"),
        axis.title = element_text(size = base_size + 1, face = "bold"),
        axis.text = element_text(size = base_size - 1, color = "black"),
        plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
        legend.title = element_text(size = base_size, face = "bold"),
        legend.text = element_text(size = base_size - 1),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.3))
}

W <- 6.5; DPI <- 300
sf <- function(p, fn, h = 4.5) {
  ggsave(file.path(OUT_DIR, fn), p, width = W, height = h, dpi = DPI, device = "png", bg = "white")
  cat("  [OK]", fn, "\n")
}

cat("=== Regenerating Figures ===\n")

# Fig2: SHAP bar + beeswarm
cat("[2/6] SHAP\n")
sh <- read.csv(file.path(DATA, "shap.csv"))
sh <- sh[order(-sh$MeanAbsSHAP), ]
top20 <- head(sh, 20); top20$Gene <- factor(top20$Gene, levels = rev(top20$Gene))
sf(ggplot(top20, aes(x = Gene, y = MeanAbsSHAP)) +
   geom_bar(stat = "identity", fill = "#378ADD") + coord_flip() +
   labs(x = NULL, y = "Mean |SHAP|") + theme_publication(), "Fig2_SHAP_bar.png", 4)

top10 <- head(sh, 10); top10$Gene <- factor(top10$Gene, levels = rev(top10$Gene))
sf(ggplot(top10, aes(x = Gene, y = MeanAbsSHAP)) +
   geom_point(aes(size = MeanAbsSHAP), color = "#378ADD", alpha = 0.7) + coord_flip() +
   labs(x = NULL, y = "SHAP Value") + theme_publication() + theme(legend.position = "none"),
   "Fig2_SHAP_beeswarm.png", 4)

# Fig3: Pathway volcano (using delta as effect size)
cat("[3/6] Volcano\n")
pd <- read.csv(file.path(DATA, "pathway_diff.csv"))
pd$sig <- ifelse(pd$p_value < 0.05 & abs(pd$delta) > 0.01, "Sig",
                 ifelse(pd$p_value < 0.1, "Borderline", "NS"))
sf(ggplot(pd, aes(x = delta, y = -log10(p_value), color = sig)) +
   geom_point(size = 2, alpha = 0.8) +
   scale_color_manual(values = c("Sig"="#D85A30","Borderline"="#F0997B","NS"="#B4B2A9")) +
   geom_hline(yintercept = -log10(0.1), linetype = "dashed") +
   geom_vline(xintercept = c(-0.01, 0.01), linetype = "dashed") +
   labs(x = "Delta (R - S)", y = "-log10(p)") +
   theme_publication() + theme(legend.position = "bottom"), "Fig3A_pathway_volcano.png", 4.5)

# Fig4: Cox forest
cat("[4/6] Cox forest\n")
cr <- read.csv(file.path(DATA, "cox_results.csv"))
cr$Variable <- factor(cr$Variable, levels = rev(cr$Variable))
sf(ggplot(cr, aes(x = HR, y = Variable)) +
   geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
   geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper), height = 0.2, color = "#378ADD") +
   geom_point(aes(color = p_value < 0.05), size = 3) +
   scale_color_manual(values = c("TRUE"="#D85A30","FALSE"="#378ADD")) +
   scale_x_log10() + labs(x = "HR (95% CI)", y = NULL) +
   theme_publication() + theme(legend.position = "none"), "Fig4C_cox_forest.png", 3.5)

# Fig6: Fair comparison
cat("[6/6] Fair comparison\n")
fc <- read.csv(file.path(DATA, "fair_compare.csv"))
c_gene <- as.numeric(fc[fc$Metric == "C-corrected", "Gene"])
c_path <- as.numeric(fc[fc$Metric == "C-corrected", "Pathway"])
fc_df <- data.frame(Level = c("Gene-Level", "Pathway-Level"), C = c(c_gene, c_path))
sf(ggplot(fc_df, aes(x = Level, y = C, fill = Level)) +
   geom_bar(stat = "identity") + geom_hline(yintercept = 0.5, linetype = "dashed") +
   scale_fill_manual(values = c("Gene-Level"="#B4B2A9","Pathway-Level"="#378ADD")) +
   labs(y = "Corrected C-index", x = NULL) + ylim(0, 0.6) +
   theme_publication() + theme(legend.position = "none"), "Fig6B_method_comparison.png", 3.5)

cat("\n=== Done ===\n")
cat("Output:", OUT_DIR, "\n")
cat("Copy to HD: cp", OUT_DIR, "/*.png <PROJECT>/results/figures_png/hd/\n")
