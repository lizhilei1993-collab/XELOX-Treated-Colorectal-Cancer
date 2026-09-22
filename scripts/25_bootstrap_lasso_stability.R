# ============================================================
# 25_BOOTSTRAP_LASSO_STABILITY.R
# Path A: Bootstrap Stability Analysis for LASSO-PRS Model
#
# Purpose:
#   Diagnose the stability of LASSO feature selection in the
#   44-pathway × 277-sample merged training set (GSE39582 +
#   GSE104645). B=1000 bootstrap resamples, each running the
#   exact same pipeline as Script 16 (univariate AUC pre-filter
#   + cv.glmnet lambda.1se).
#
# Key Question:
#   Are the 3 selected pathways (Chemokine/ECM/BER) consistently
#   chosen across bootstrap samples, or are they artifacts of a
#   single data split?
#
# Input:
#   GSE39582_pathway_scores.rds         — 44 pathways x 585 samples
#   GSE104645_pathway_scores.rds        — 44 pathways x 193 samples
#   GSE39582_xelox_groups.csv           — response labels
#   GSE104645_clinical_data.csv         — RECIST response
#
# Output (C:/xelox_work/results/tables/bootstrap_stability/):
#   bootstrap_selection_freq.csv        — selection frequency table
#   bootstrap_coefficient_dist.csv      — coefficient distribution
#   bootstrap_stability_report.txt      — text summary
#   bootstrap_corr_matrix.csv           — pairwise selection correlation
#
# Output (C:/xelox_work/results/figures/bootstrap_stability/):
#   bootstrap_selection_freq_bar.pdf    — selection frequency barplot
#   bootstrap_coef_heatmap.pdf          — coefficient heatmap (top 20)
#   bootstrap_stability_curve.pdf       — cumulative stability curve
# ============================================================

Sys.setenv(TMPDIR = "C:/temp", TMP = "C:/temp", TEMP = "C:/temp")
.libPaths(c("C:/Rlibs", .libPaths()))

PROJECT_ROOT_XELOX <- "/path/to/xelox_project"
RESULTS_TAB_DIR     <- file.path(PROJECT_ROOT_XELOX, "results", "tables")

OUT_TAB_DIR <- file.path(PROJECT_ROOT_XELOX, "results", "tables", "bootstrap_stability")
OUT_FIG_DIR <- file.path(PROJECT_ROOT_XELOX, "results", "figures", "bootstrap_stability")
dir.create(OUT_TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_FIG_DIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("Path A: Bootstrap LASSO Stability Diagnosis\n")
cat("  B = 1000 bootstrap resamples\n")
cat("  Pipeline: AUC pre-filter → cv.glmnet (lambda.1se)\n")
cat("============================================================\n\n")

# ============================================================
# 0. Load packages
# ============================================================
cat("=== [0] Loading packages ===\n\n")

required_pkgs <- c("glmnet", "pROC", "ggplot2", "data.table", "reshape2")

for (pkg in required_pkgs) {
  loaded <- require(pkg, lib.loc = "C:/Rlibs", character.only = TRUE, quietly = TRUE)
  if (!loaded) {
    loaded <- require(pkg, character.only = TRUE, quietly = TRUE)
  }
  if (loaded) {
    cat(sprintf("  + %s\n", pkg))
  } else {
    stop(sprintf("  x %s - required package not found", pkg))
  }
}
cat("\n")

# ============================================================
# 1. Load training data (exactly as Script 16)
# ============================================================
cat("=== [1] Loading training data ===\n\n")

# Pathway scores
pa_dir <- file.path(RESULTS_TAB_DIR, "pathway_activity")
ps_39582  <- readRDS(file.path(pa_dir, "GSE39582_pathway_scores.rds"))
ps_104645 <- readRDS(file.path(pa_dir, "GSE104645_pathway_scores.rds"))

cat(sprintf("  GSE39582 pathway scores: %d pathways x %d samples\n",
            nrow(ps_39582), ncol(ps_39582)))
cat(sprintf("  GSE104645 pathway scores: %d pathways x %d samples\n",
            nrow(ps_104645), ncol(ps_104645)))

# Response labels: GSE39582 (RFS proxy)
groups_39582 <- read.csv(file.path(RESULTS_TAB_DIR, "GSE39582_xelox_groups.csv"),
                         stringsAsFactors = FALSE)
rownames(groups_39582) <- groups_39582$sample_id

g39582_binary <- groups_39582$group
names(g39582_binary) <- groups_39582$sample_id
g39582_keep <- g39582_binary %in% c("sensitive", "resistant")
g39582_y <- ifelse(g39582_binary[g39582_keep] == "resistant", 1, 0)
g39582_common <- intersect(names(g39582_y), colnames(ps_39582))
g39582_y <- g39582_y[g39582_common]

cat(sprintf("  GSE39582 XELOX: Sens=%d, Res=%d (RFS proxy)\n",
            sum(g39582_y == 0), sum(g39582_y == 1)))

# Response labels: GSE104645 (RECIST)
clin_104645 <- read.csv(file.path(RESULTS_TAB_DIR, "GSE104645_clinical_data.csv"),
                         row.names = 1, stringsAsFactors = FALSE)

regimen_col <- "X1st.line.chemotherapy.regimens.ch1"
response_col <- "best.response.of.1st.line.chemotherapy.ch1"

is_oxali <- grepl("XELOX|FOLFOX|SOX|CAPOX",
                   clin_104645[[regimen_col]], ignore.case = TRUE)
oxali_samples <- rownames(clin_104645)[is_oxali]

resp <- clin_104645[[response_col]]
names(resp) <- rownames(clin_104645)

g104645_y <- rep(NA, length(oxali_samples))
names(g104645_y) <- oxali_samples
g104645_y[resp[oxali_samples] %in% c("Complete response", "Partial response")] <- 0
g104645_y[resp[oxali_samples] %in% c("Progressive disease", "Stable disease")] <- 1
g104645_y <- g104645_y[!is.na(g104645_y)]
g104645_common <- intersect(names(g104645_y), colnames(ps_104645))
g104645_y <- g104645_y[g104645_common]

cat(sprintf("  GSE104645 oxali: Sens=%d, Res=%d (RECIST)\n",
            sum(g104645_y == 0), sum(g104645_y == 1)))

# Build merged training matrix
train_ps <- cbind(
  ps_39582[, names(g39582_y), drop = FALSE],
  ps_104645[, names(g104645_y), drop = FALSE]
)
train_y <- c(g39582_y, g104645_y)
train_source <- c(rep("GSE39582", length(g39582_y)),
                  rep("GSE104645", length(g104645_y)))

pathway_names <- rownames(train_ps)
n_pathways <- length(pathway_names)
n_samples <- ncol(train_ps)
n_39582 <- length(g39582_y)

cat(sprintf("  Merged training: %d pathways x %d samples\n", n_pathways, n_samples))
cat(sprintf("  Resist/Sens: %d/%d (%.1f%% resistant)\n",
            sum(train_y == 1), sum(train_y == 0), mean(train_y) * 100))
cat(sprintf("  Source: GSE39582=%d, GSE104645=%d\n\n", n_39582, n_samples - n_39582))

# ============================================================
# 2. Helper: Univariate AUC computation
# ============================================================
compute_auc <- function(x, y) {
  tryCatch({
    roc_obj <- roc(y, x, direction = "<", quiet = TRUE)
    as.numeric(auc(roc_obj))
  }, error = function(e) NA)
}

# ============================================================
# 3. Bootstrap LASSO (B = 1000)
# ============================================================
cat("=== [2] Bootstrap LASSO (B = 1000) ===\n\n")

intermediate_file <- file.path(OUT_TAB_DIR, "bootstrap_intermediate.rds")

if (file.exists(intermediate_file)) {
  cat("  Loading saved intermediate results (skip bootstrap)...\n")
  inter_data <- readRDS(intermediate_file)
  selection_matrix <- inter_data$selection_matrix
  coef_matrix      <- inter_data$coef_matrix
  n_selected_vec   <- inter_data$n_selected_vec
  bootstrap_auc    <- inter_data$bootstrap_auc
  lambda_used      <- inter_data$lambda_used
  pathway_names    <- inter_data$pathway_names
  train_ps         <- inter_data$train_ps
  train_y          <- inter_data$train_y
  train_source     <- inter_data$train_source
  n_samples        <- ncol(train_ps)
  n_39582          <- sum(train_source == "GSE39582")
  B                <- nrow(selection_matrix)
  n_pathways       <- length(pathway_names)
  cat(sprintf("  Loaded B=%d, %d pathways x %d samples\n", B, n_pathways, n_samples))
} else {

B <- 1000
n_total <- n_samples

# Storage
selection_matrix <- matrix(0, nrow = B, ncol = n_pathways)
colnames(selection_matrix) <- pathway_names
coef_matrix      <- matrix(NA, nrow = B, ncol = n_pathways)
colnames(coef_matrix) <- pathway_names
bootstrap_auc    <- numeric(B)      # Training AUC per bootstrap
n_selected_vec   <- integer(B)      # Number of pathways selected per bootstrap
lambda_used      <- character(B)    # "lambda.1se" or "lambda.min"

set.seed(42)

for (b in seq_len(B)) {
  # --- 3a. Resample with replacement ---
  idx <- sample(n_total, size = n_total, replace = TRUE)
  
  boot_x_raw <- train_ps[, idx, drop = FALSE]
  boot_y     <- train_y[idx]
  boot_src   <- train_source[idx]
  
  n_boot_39582 <- sum(boot_src == "GSE39582")
  
  # --- 3b. Univariate AUC pre-filter on GSE39582 bootstrap subset ---
  g39582_in_boot <- which(boot_src == "GSE39582")
  
  if (length(g39582_in_boot) < 10) {
    # Too few GSE39582 samples for meaningful AUC
    next
  }
  
  auc_vec <- apply(boot_x_raw, 1, function(pw) {
    compute_auc(pw[g39582_in_boot], boot_y[g39582_in_boot])
  })
  names(auc_vec) <- pathway_names
  
  # Filter: |AUC - 0.5| > 0.05
  auc_filter <- abs(auc_vec - 0.5) > 0.05
  auc_filter[is.na(auc_filter)] <- FALSE
  
  if (sum(auc_filter) < 2) {
    # Too few pathways pass filter
    next
  }
  
  boot_x_filtered <- t(boot_x_raw[auc_filter, , drop = FALSE])
  
  # --- 3c. Scale and run cv.glmnet ---
  boot_x_scaled <- scale(boot_x_filtered)
  
  cv_fit <- tryCatch({
    cv.glmnet(
      x = boot_x_scaled,
      y = boot_y,
      family = "binomial",
      alpha = 1,
      nfolds = min(10, n_total),
      type.measure = "deviance"
    )
  }, error = function(e) NULL)
  
  if (is.null(cv_fit)) next
  
  # --- 3d. Select at lambda.1se ---
  coef_1se <- as.matrix(coef(cv_fit, s = "lambda.1se"))
  idx_1se <- which(coef_1se[-1, 1] != 0)
  
  if (length(idx_1se) == 0) {
    # Fallback to lambda.min
    coef_min <- as.matrix(coef(cv_fit, s = "lambda.min"))
    idx_sel <- which(coef_min[-1, 1] != 0)
    lambda_used[b] <- "lambda.min"
    coef_vals <- coef_min[idx_sel + 1, 1]
  } else {
    idx_sel <- idx_1se
    lambda_used[b] <- "lambda.1se"
    coef_vals <- coef_1se[idx_sel + 1, 1]
  }
  
  # Map back to full pathway list
  filtered_names <- pathway_names[auc_filter]
  selected_names <- filtered_names[idx_sel]
  
  selection_matrix[b, selected_names] <- 1
  coef_matrix[b, selected_names] <- coef_vals
  n_selected_vec[b] <- length(selected_names)
  
  # Compute training AUC
  if (length(selected_names) > 0) {
    pred <- predict(cv_fit, newx = boot_x_scaled, s = "lambda.1se", type = "response")
    bootstrap_auc[b] <- compute_auc(as.numeric(pred), boot_y)
  }
  
  # Progress
  if (b %% 100 == 0) {
    cat(sprintf("  Bootstrap %4d/%d | Avg selected: %.1f pathways | Fallback rate: %.1f%%\n",
                b, B, mean(n_selected_vec[1:b], na.rm = TRUE),
                mean(lambda_used[1:b] == "lambda.min", na.rm = TRUE) * 100))
  }
}

cat(sprintf("\n  Completed %d successful bootstrap iterations\n", B))
cat(sprintf("  Mean pathways selected: %.2f ± %.2f\n",
            mean(n_selected_vec, na.rm = TRUE), sd(n_selected_vec, na.rm = TRUE)))
cat(sprintf("  lambda.min fallback rate: %.1f%%\n",
            mean(lambda_used == "lambda.min", na.rm = TRUE) * 100))

# Save intermediate data to avoid re-running bootstrap
saveRDS(list(selection_matrix = selection_matrix, coef_matrix = coef_matrix,
             n_selected_vec = n_selected_vec, bootstrap_auc = bootstrap_auc,
             lambda_used = lambda_used, pathway_names = pathway_names,
             train_ps = train_ps, train_y = train_y, train_source = train_source),
        file.path(OUT_TAB_DIR, "bootstrap_intermediate.rds"))
cat("  -> Saved intermediate data: bootstrap_intermediate.rds\n")
}  # end else (bootstrap from scratch)

# ============================================================
# 4. Compute stability metrics
# ============================================================
cat("\n=== [3] Stability Metrics ===\n\n")

# 4a. Selection frequency
sel_freq <- colSums(selection_matrix) / B
sel_freq_df <- data.frame(
  Pathway = pathway_names,
  Selection_Frequency = round(sel_freq, 4),
  N_Selected = colSums(selection_matrix),
  stringsAsFactors = FALSE
)
sel_freq_df <- sel_freq_df[order(sel_freq_df$Selection_Frequency, decreasing = TRUE), ]

cat("  Top 15 pathways by selection frequency:\n")
for (i in 1:min(15, nrow(sel_freq_df))) {
  cat(sprintf("    %2d. %-50s %.1f%% (%d/%d)\n",
              i, sel_freq_df$Pathway[i],
              sel_freq_df$Selection_Frequency[i] * 100,
              sel_freq_df$N_Selected[i], B))
}

# 4b. Coefficient distribution for frequently selected pathways
coef_summary <- data.frame(
  Pathway = pathway_names,
  Mean_Coef = round(colMeans(coef_matrix, na.rm = TRUE), 4),
  Median_Coef = round(apply(coef_matrix, 2, median, na.rm = TRUE), 4),
  SD_Coef = round(apply(coef_matrix, 2, sd, na.rm = TRUE), 4),
  Q025_Coef = round(apply(coef_matrix, 2, quantile, 0.025, na.rm = TRUE), 4),
  Q975_Coef = round(apply(coef_matrix, 2, quantile, 0.975, na.rm = TRUE), 4),
  Selection_Freq = sel_freq,
  stringsAsFactors = FALSE
)
coef_summary <- coef_summary[order(coef_summary$Selection_Freq, decreasing = TRUE), ]

# 4c. Stability score (Jaccard-like pairwise similarity)
# Efficient implementation: sample random pairs instead of all O(B^2)
pairwise_jaccard <- function(mat, max_pairs = 50000) {
  # Matrix is B x n_pathways, binary
  B <- nrow(mat)
  if (B < 2) return(NA)
  
  n_pairs <- B * (B - 1) / 2
  if (n_pairs <= max_pairs) {
    # Full computation (pre-allocated)
    jac <- numeric(n_pairs)
    idx <- 0
    for (i in 1:(B - 1)) {
      row_i <- mat[i, ]
      for (j in (i + 1):B) {
        inter <- sum(row_i & mat[j, ])
        uni   <- sum(row_i | mat[j, ])
        idx <- idx + 1
        if (uni > 0) jac[idx] <- inter / uni
      }
    }
    return(mean(jac[1:idx], na.rm = TRUE))
  } else {
    # Random sample max_pairs pairs
    set.seed(123)
    pairs <- replicate(max_pairs, sample(B, 2, replace = FALSE))
    jac <- numeric(max_pairs)
    for (k in seq_len(max_pairs)) {
      i <- pairs[1, k]; j <- pairs[2, k]
      inter <- sum(mat[i, ] & mat[j, ])
      uni   <- sum(mat[i, ] | mat[j, ])
      if (uni > 0) jac[k] <- inter / uni
    }
    return(mean(jac, na.rm = TRUE))
  }
}

stability_jaccard <- pairwise_jaccard(selection_matrix)
cat(sprintf("\n  Mean pairwise Jaccard stability: %.4f\n", stability_jaccard))

# 4d. Number of unique pathways ever selected
n_unique <- sum(sel_freq > 0)
cat(sprintf("  Unique pathways ever selected: %d / %d (%.1f%%)\n",
            n_unique, n_pathways, n_unique / n_pathways * 100))

# 4e. Pathways with frequency > 50% (moderate stability)
stable_pathways <- sel_freq_df$Pathway[sel_freq_df$Selection_Frequency > 0.5]
cat(sprintf("  Pathways with >50%% selection frequency: %d\n", length(stable_pathways)))
if (length(stable_pathways) > 0) {
  for (pw in stable_pathways) {
    cat(sprintf("    - %s: %.1f%%\n", pw, sel_freq[sel_freq_df$Pathway == pw] * 100))
  }
}

# 4f. Check original LASSO 3 pathways
orig_pw <- c("KEGG_CHEMOKINE_SIGNALING_PATHWAY",
             "KEGG_ECM_RECEPTOR_INTERACTION",
             "KEGG_BASE_EXCISION_REPAIR")
cat("\n  Original LASSO (Script 16) 3 pathways stability:\n")
for (pw in orig_pw) {
  cat(sprintf("    %-50s: %.1f%%\n", pw, sel_freq[pw] * 100))
}

# ============================================================
# 5. Save tables
# ============================================================
cat("\n=== [4] Saving tables ===\n\n")

write.csv(sel_freq_df, file.path(OUT_TAB_DIR, "bootstrap_selection_freq.csv"),
          row.names = FALSE)
cat("  -> Saved: bootstrap_selection_freq.csv\n")

write.csv(coef_summary, file.path(OUT_TAB_DIR, "bootstrap_coefficient_dist.csv"),
          row.names = FALSE)
cat("  -> Saved: bootstrap_coefficient_dist.csv\n")

# Pairwise correlation of selections (for heatmap / network)
top_n <- min(20, sum(sel_freq > 0.01))
top_pathways <- sel_freq_df$Pathway[1:top_n]
corr_sel <- cor(selection_matrix[, top_pathways, drop = FALSE], use = "pairwise.complete.obs")
write.csv(corr_sel, file.path(OUT_TAB_DIR, "bootstrap_corr_matrix.csv"))
cat("  -> Saved: bootstrap_corr_matrix.csv\n")

# Text summary
sink(file.path(OUT_TAB_DIR, "bootstrap_stability_report.txt"))
cat("============================================================\n")
cat("Bootstrap LASSO Stability Report\n")
cat("============================================================\n\n")
cat(sprintf("Training set: %d samples (%d GSE39582 + %d GSE104645)\n",
            n_samples, n_39582, n_samples - n_39582))
cat(sprintf("Features: %d KEGG + HALLMARK pathways\n", n_pathways))
cat(sprintf("Bootstrap iterations: B = %d\n", B))
cat(sprintf("Pipeline: Univariate AUC pre-filter (|AUC-0.5|>0.05) + cv.glmnet(lambda.1se)\n\n"))

cat(sprintf("Mean pathways selected per bootstrap: %.2f +- %.2f\n",
            mean(n_selected_vec, na.rm = TRUE), sd(n_selected_vec, na.rm = TRUE)))
cat(sprintf("lambda.min fallback rate: %.1f%%\n",
            mean(lambda_used == "lambda.min", na.rm = TRUE) * 100))
cat(sprintf("Mean pairwise Jaccard stability: %.4f\n", stability_jaccard))
cat(sprintf("Unique pathways ever selected: %d / %d (%.1f%%)\n\n",
            n_unique, n_pathways, n_unique / n_pathways * 100))

cat("Pathways with >50% selection frequency:\n")
if (length(stable_pathways) > 0) {
  for (pw in stable_pathways) {
    cat(sprintf("  %s: %.1f%%\n", pw, sel_freq[pw] * 100))
  }
} else {
  cat("  NONE\n")
}

cat("\nOriginal Script 16 LASSO pathways:\n")
for (pw in orig_pw) {
  coef_row <- coef_summary[coef_summary$Pathway == pw, ]
  cat(sprintf("  %s: freq=%.1f%%, coef=%.4f [%.4f, %.4f]\n",
              pw, sel_freq[pw] * 100,
              coef_row$Median_Coef,
              coef_row$Q025_Coef, coef_row$Q975_Coef))
}

cat("\nTop 15 pathways by selection frequency:\n")
for (i in 1:min(15, nrow(sel_freq_df))) {
  coef_row <- coef_summary[i, ]
  cat(sprintf("  %2d. %-50s %.1f%%  coef=%.4f [%.4f,%.4f]\n",
              i, sel_freq_df$Pathway[i],
              sel_freq_df$Selection_Frequency[i] * 100,
              coef_row$Median_Coef,
              coef_row$Q025_Coef, coef_row$Q975_Coef))
}

cat("\nDiagnosis:\n")
if (length(stable_pathways) >= 3) {
  cat("  STABLE: >=3 pathways selected >50% of the time.\n")
  cat("  The LASSO signal appears robust. Consider Elastic Net for more conservative selection.\n")
} else if (length(stable_pathways) >= 1) {
  cat("  MODERATE: 1-2 pathways selected >50% of the time.\n")
  cat("  Limited stable signal. Consider feature space expansion or alternative modeling.\n")
} else {
  cat("  UNSTABLE: No pathway selected >50% of the time.\n")
  cat("  The LASSO signal is not stable. The original 3-pathway model likely\n")
  cat("  reflects a single random split result rather than a true underlying\n")
  cat("  biological signal. Consider switching to Elastic Net (alpha=0.5) or\n")
  cat("  gene-level features with Ridge/LASSO.\n")
}

sink()
cat("  -> Saved: bootstrap_stability_report.txt\n")

# ============================================================
# 6. Generate figures
# ============================================================
cat("\n=== [5] Generating figures ===\n\n")

# 6a. Selection frequency barplot (top 30)
n_show <- min(30, sum(sel_freq > 0))
if (n_show > 0) {
  top_sel <- head(sel_freq_df, n_show)
  top_sel$Pathway_short <- sapply(top_sel$Pathway, function(x) {
    s <- gsub("^KEGG_", "", gsub("^HALLMARK_", "", x))
    if (nchar(s) > 45) paste0(substr(s, 1, 42), "...") else s
  })
  top_sel$Pathway_short <- factor(top_sel$Pathway_short,
                                   levels = rev(top_sel$Pathway_short))
  
  # Color: original 3 pathways in red
  top_sel$is_original <- top_sel$Pathway %in% orig_pw
  
  pdf(file.path(OUT_FIG_DIR, "bootstrap_selection_freq_bar.pdf"),
      width = 12, height = 8)
  
  p <- ggplot(top_sel, aes(x = Selection_Frequency * 100, y = Pathway_short,
                            fill = is_original)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
                      labels = c("TRUE" = "Original LASSO", "FALSE" = "Other"),
                      name = NULL) +
    geom_vline(xintercept = 50, linetype = "dashed", color = "gray50") +
    annotate("text", x = 51, y = n_show + 0.5, label = "50% threshold",
             hjust = 0, size = 3, color = "gray50") +
    labs(title = "Bootstrap LASSO Pathway Selection Frequency",
         subtitle = sprintf("B=%d bootstrap | %d pathways ever selected | Jaccard stability=%.4f",
                            B, n_unique, stability_jaccard),
         x = "Selection Frequency (%)", y = "") +
    theme_bw(base_size = 11) +
    theme(legend.position = "top")
  
  print(p)
  dev.off()
  cat("  -> Saved: bootstrap_selection_freq_bar.pdf\n")
}

# 6b. Coefficient heatmap (top 20 most selected)
n_heat <- min(20, n_unique)
if (n_heat > 0) {
  # Extract coefficient matrix for top pathways
  top_heat_pathways <- sel_freq_df$Pathway[1:n_heat]
  coef_heat <- coef_matrix[, top_heat_pathways, drop = FALSE]
  
  # For display, sample 200 bootstrap iterations
  disp_n <- min(200, B)
  coef_disp <- coef_heat[1:disp_n, , drop = FALSE]
  
  # Melt for ggplot
  coef_long <- melt(coef_disp, varnames = c("Bootstrap", "Pathway"),
                    value.name = "Coefficient")
  coef_long <- coef_long[!is.na(coef_long$Coefficient), ]
  
  # Shorten pathway names
  coef_long$Pathway_short <- sapply(as.character(coef_long$Pathway), function(x) {
    s <- gsub("^KEGG_", "", gsub("^HALLMARK_", "", x))
    if (nchar(s) > 40) paste0(substr(s, 1, 37), "...") else s
  })
  
  if (nrow(coef_long) > 0) {
    pdf(file.path(OUT_FIG_DIR, "bootstrap_coef_heatmap.pdf"),
        width = 14, height = 8)
    
    p2 <- ggplot(coef_long, aes(x = Bootstrap, y = Pathway_short, fill = Coefficient)) +
      geom_tile() +
      scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",
                           midpoint = 0, name = "Coefficient") +
      labs(title = "Bootstrap LASSO Coefficient Stability (Top 20 Pathways)",
           subtitle = sprintf("B=%d iterations (%d displayed) | Red=Pro-Resistance, Blue=Pro-Sensitive",
                              B, disp_n),
           x = "Bootstrap Iteration", y = "") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
    
    print(p2)
    dev.off()
    cat("  -> Saved: bootstrap_coef_heatmap.pdf\n")
  }
}

# 6c. Cumulative stability curve (sampled at intervals)
cat("  Generating cumulative stability curve...\n")
cum_jaccard <- rep(NA_real_, B)
cum_steps <- unique(c(seq(10, B, by = 10), B))
if (B > 1) {
  for (k in cum_steps) {
    cum_jaccard[k] <- tryCatch(
      pairwise_jaccard(selection_matrix[1:k, ]),
      error = function(e) NA_real_
    )
  }
  # Interpolate for smooth curve (fill forward)
  last_val <- NA_real_
  for (i in seq_len(B)) {
    if (!is.na(cum_jaccard[i])) {
      last_val <- cum_jaccard[i]
    } else {
      cum_jaccard[i] <- last_val
    }
  }
  
  pdf(file.path(OUT_FIG_DIR, "bootstrap_stability_curve.pdf"),
      width = 8, height = 5)
  
  y_vals <- cum_jaccard[!is.na(cum_jaccard)]
  if (length(y_vals) > 1) {
    plot(cum_jaccard, type = "l", lwd = 2, col = "#4575b4",
         xlab = "Number of Bootstrap Iterations",
         ylab = "Mean Pairwise Jaccard Stability",
         main = "Cumulative LASSO Stability Over Bootstrap Iterations")
    abline(h = 0.5, lty = 2, col = "gray50")
    legend("bottomright", legend = c("Stability curve", "50% threshold"),
           col = c("#4575b4", "gray50"), lty = c(1, 2), lwd = c(2, 1))
  } else {
    plot(1, 1, type = "n", xlab = "", ylab = "",
         main = "Cumulative LASSO Stability - No data")
    text(1, 1, "Insufficient data for stability curve")
  }
  
  dev.off()
  cat("  -> Saved: bootstrap_stability_curve.pdf\n")
}

# ============================================================
# 7. Summary
# ============================================================
cat("\n=== [6] Summary ===\n\n")

cat(sprintf("Bootstrap iterations: %d\n", B))
cat(sprintf("Mean pathways selected: %.2f +- %.2f\n",
            mean(n_selected_vec, na.rm = TRUE), sd(n_selected_vec, na.rm = TRUE)))
cat(sprintf("lambda.min fallback: %.1f%%\n",
            mean(lambda_used == "lambda.min", na.rm = TRUE) * 100))
cat(sprintf("Jaccard stability: %.4f\n", stability_jaccard))
cat(sprintf("Unique pathways: %d / %d\n", n_unique, n_pathways))
cat(sprintf(">50%% stable pathways: %d\n", length(stable_pathways)))

cat("\n  Output tables:\n")
for (f in list.files(OUT_TAB_DIR)) {
  cat(sprintf("    %s\n", f))
}
cat("\n  Output figures:\n")
for (f in list.files(OUT_FIG_DIR)) {
  cat(sprintf("    %s\n", f))
}

cat("\n=== Bootstrap LASSO Stability Complete ===\n\n")
