#!/usr/bin/env Rscript
.libPaths(c('C:/tmp/rlib', .libPaths()))
library(ggplot2)
# Unified ggplot2 theme for publication figures
# Run this BEFORE generating any figure:
#   source("scripts/unified_figure_theme.R")
#
# Then call: ggplot(...) + theme_publication()
#
# To regenerate all existing figures with this theme, 
# add + theme_publication() to each ggplot call in the R scripts.

library(ggplot2)

theme_publication <- function(base_size = 11, base_family = "serif") {
  theme_bw(base_size = base_size, base_family = base_family) +
  theme(
    # Text
    text = element_text(family = base_family, size = base_size, color = "black"),
    axis.title = element_text(size = base_size + 1, face = "bold"),
    axis.text = element_text(size = base_size - 1, color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = base_size, hjust = 0.5),
    legend.title = element_text(size = base_size, face = "bold"),
    legend.text = element_text(size = base_size - 1),
    strip.text = element_text(size = base_size, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "black"),
    
    # Grid
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    
    # Legend
    legend.position = "right",
    legend.key = element_blank(),
    
    # Axis
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    
    # Margins
    plot.margin = margin(5, 5, 5, 5)
  )
}

# Save file dimensions (for journals)
# Single column:  width = 3.5 in (89 mm)
# 1.5 column:     width = 5.0 in (127 mm)
# Full page:      width = 7.0 in (178 mm)
# Height:         typically 2/3 to 3/4 of width

message("Unified publication theme loaded. Use: ggplot(...) + theme_publication()")
message("Suggested dimensions for ggsave():")
message("  Single col: width=3.5, height=2.5-3.0")
message("  Double col: width=6.5, height=4.0-5.0")
message("  DPI: 300 (ggsave(dpi=300))")
