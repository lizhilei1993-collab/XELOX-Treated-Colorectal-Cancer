# ============================================================
# run_13_combat.R - Wrapper for ComBat batch correction (v6.0 fix)
# Core fix: Drop NA samples instead of filling with 0.5
# Reads from /path/to/xelox_project, prints results to stdout
# ============================================================

Sys.setenv(TMPDIR = "/tmp", TMP = "/tmp", TEMP = "/tmp")
.libPaths(c("/path/to/Rlibs", .libPaths()))
if (.Platform$OS.type == "windows") {
  tryCatch(Sys.setlocale("LC_ALL", "Chinese"), error = function(e) {
    tryCatch(Sys.setlocale("LC_ALL", "chs"), error = function(e2) {})
  })
}

suppressPackageStartupMessages({
  library(Biobase)
  library(sva)
  library(limma)
})

READ_ROOT <- "/path/to/xelox_project"
DATA_GEO_DIR <- file.path(READ_ROOT, "data", "geo")
READ_TAB_DIR <- file.path(READ_ROOT, "results", "tables")

cat("=== ComBat Batch Correction (v6.0 wrapper) ===\n\n")
cat(sprintf("  READ_ROOT: %s (exists=%s)\n", READ_ROOT, file.exists(READ_ROOT)))

# ============================================================
# 1. Load expression matrices
# ============================================================
cat("\n[1] Loading expression matrices\n")

cat("  Loading GSE39582...\n")
eset_39582 <- readRDS(file.path(DATA_GEO_DIR, "GSE39582_eset.rds"))
expr_39582 <- exprs(eset_39582)
cat(sprintf("    GSE39582: %d probes x %d samples\n", nrow(expr_39582), ncol(expr_39582)))

cat("  Loading GSE72970...\n")
expr_72970 <- readRDS(file.path(DATA_GEO_DIR, "GSE72970_expression.rds"))
cat(sprintf("    GSE72970: %d probes x %d samples\n", nrow(expr_72970), ncol(expr_72970)))

cat("  Loading GSE69657...\n")
expr_69657 <- readRDS(file.path(DATA_GEO_DIR, "GSE69657_expression.rds"))
cat(sprintf("    GSE69657: %d probes x %d samples\n", nrow(expr_69657), ncol(expr_69657)))

# ============================================================
# 2. Build combined matrix
# ============================================================
cat("\n[2] Building combined matrix\n")

common_probes <- Reduce(intersect, list(
  rownames(expr_39582), rownames(expr_72970), rownames(expr_69657)
))
cat(sprintf("  Common probes: %d\n", length(common_probes)))

X_39582 <- expr_39582[common_probes, ]
X_72970 <- expr_72970[common_probes, ]
X_69657 <- expr_69657[common_probes, ]
X_combined <- cbind(X_39582, X_72970, X_69657)
cat(sprintf("  Combined: %d probes x %d samples\n", nrow(X_combined), ncol(X_combined)))

batch <- c(rep("GSE39582", ncol(X_39582)),
           rep("GSE72970", ncol(X_72970)),
           rep("GSE69657", ncol(X_69657)))

# ============================================================
# 3. Build response labels
# ============================================================
cat("\n[3] Building response labels\n")

clin_39582 <- read.csv(file.path(READ_TAB_DIR, "GSE39582_clinical_data.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
clin_72970 <- read.csv(file.path(READ_TAB_DIR, "GSE72970_clinical_data.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
clin_69657 <- read.csv(file.path(READ_TAB_DIR, "GSE69657_clinical_data.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)

# GSE39582 response
if ("response_status" %in% colnames(clin_39582)) {
  resp_39582 <- ifelse(clin_39582$response_status == "resistant", 1,
                       ifelse(clin_39582$response_status == "sensitive", 0, NA))
} else if ("Relapse" %in% colnames(clin_39582)) {
  resp_39582 <- clin_39582$Relapse
} else {
  resp_39582 <- rep(NA, nrow(clin_39582))
}

# GSE72970 response
resp_col_72970 <- grep("response status", names(clin_72970), value = TRUE)[1]
resp_72970 <- ifelse(clin_72970[[resp_col_72970]] == "NR", 1, 0)

# GSE69657 response
resp_col_69657 <- grep("chemoresponse", names(clin_69657), value = TRUE)[1]
if (length(resp_col_69657) == 0) resp_col_69657 <- grep("response", names(clin_69657), value = TRUE)[1]
resp_69657 <- ifelse(clin_69657[[resp_col_69657]] == "noresponder", 1,
                     ifelse(clin_69657[[resp_col_69657]] == "responder", 0, NA))

# Align
idx_39582 <- match(colnames(X_39582), clin_39582$geo_accession)
idx_72970 <- match(colnames(X_72970), clin_72970$geo_accession)
idx_69657 <- match(colnames(X_69657), clin_69657$geo_accession)

all_resp <- c(resp_39582[idx_39582], resp_72970[idx_72970], resp_69657[idx_69657])
cat(sprintf("  Response: R=%d, S=%d, NA=%d\n",
            sum(all_resp == 1, na.rm = TRUE),
            sum(all_resp == 0, na.rm = TRUE),
            sum(is.na(all_resp))))

# ============================================================
# 4. ComBat - OLD vs NEW approach
# ============================================================
cat("\n[4] ComBat batch correction\n")

# ---- OLD approach: fill NA with 0.5 ----
cat("\n  --- OLD approach (NA filled with 0.5) ---\n")
resp_filled <- all_resp
resp_filled[is.na(resp_filled)] <- 0.5
mod_old <- model.matrix(~ resp_filled)
cat(sprintf("  OLD mod: %d samples, mean_resp=%.3f\n", nrow(mod_old), mean(resp_filled)))

X_combat_old <- tryCatch(
  ComBat(dat = as.matrix(X_combined), batch = as.factor(batch),
         mod = mod_old, par.prior = TRUE, prior.plots = FALSE),
  error = function(e) { cat(sprintf("  OLD ComBat failed: %s\n", e$message)); NULL }
)

if (!is.null(X_combat_old)) {
  cat(sprintf("  OLD ComBat output: %d x %d\n", nrow(X_combat_old), ncol(X_combat_old)))
}

# ---- NEW approach (v6.0): drop NA samples ----
cat("\n  --- NEW approach (v6.0: drop NA samples) ---\n")
complete_idx <- !is.na(all_resp)
n_na <- sum(!complete_idx)
cat(sprintf("  Removing %d samples with NA response\n", n_na))

X_combined_clean <- X_combined[, complete_idx]
batch_clean <- batch[complete_idx]
all_resp_clean <- all_resp[complete_idx]
mod_new <- model.matrix(~ all_resp_clean)
cat(sprintf("  NEW mod: %d samples, mean_resp=%.3f\n", nrow(mod_new), mean(all_resp_clean)))

X_combat_new <- tryCatch(
  ComBat(dat = as.matrix(X_combined_clean), batch = as.factor(batch_clean),
         mod = mod_new, par.prior = TRUE, prior.plots = FALSE),
  error = function(e) { cat(sprintf("  NEW ComBat failed: %s\n", e$message)); NULL }
)

if (!is.null(X_combat_new)) {
  cat(sprintf("  NEW ComBat output: %d x %d\n", nrow(X_combat_new), ncol(X_combat_new)))
}

# ============================================================
# 5. Compare: PCA batch separability
# ============================================================
cat("\n[5] PCA batch separability comparison\n")

compare_pca <- function(X_mat, batch_vec, label) {
  probe_var <- apply(X_mat, 1, var, na.rm = TRUE)
  top5k <- head(order(probe_var, decreasing = TRUE), 5000)
  X_pca <- X_mat[top5k, ]
  X_scaled <- t(scale(t(X_pca), center = TRUE, scale = TRUE))
  X_scaled[is.na(X_scaled)] <- 0
  pca <- prcomp(t(X_scaled), center = FALSE, scale. = FALSE)
  var_explained <- summary(pca)$importance[2, 1:3] * 100
  
  # LDA batch classification
  pca_df <- data.frame(PC1 = pca$x[, 1], Batch = batch_vec)
  lda_acc <- tryCatch({
    lda_fit <- MASS::lda(PC1 ~ Batch, data = pca_df)
    lda_pred <- predict(lda_fit)
    mean(lda_pred$class == pca_df$Batch)
  }, error = function(e) NA)
  
  cat(sprintf("  %s: PC1=%.1f%%, PC2=%.1f%%, LDA batch acc=%.1f%%\n",
              label, var_explained[1], var_explained[2], lda_acc * 100))
  return(lda_acc)
}

# Pre-correction
cat("  Pre-correction:\n")
lda_pre <- compare_pca(X_combined, batch, "Pre-ComBat")

# Post OLD
if (!is.null(X_combat_old)) {
  cat("  Post OLD (NA=0.5):\n")
  lda_post_old <- compare_pca(X_combat_old, batch, "Post-OLD")
}

# Post NEW
if (!is.null(X_combat_new)) {
  cat("  Post NEW (drop NA):\n")
  lda_post_new <- compare_pca(X_combat_new, batch_clean, "Post-NEW")
}

# ============================================================
# 6. Compare: response signal preservation
# ============================================================
cat("\n[6] Response signal preservation\n")

check_response_signal <- function(X_mat, resp_vec, batch_vec, label) {
  # For each dataset, compute mean expression of known resistant/sensitive marker genes
  # Just check if the response labels correlate with PC1
  probe_var <- apply(X_mat, 1, var, na.rm = TRUE)
  top5k <- head(order(probe_var, decreasing = TRUE), 5000)
  X_pca <- X_mat[top5k, ]
  X_scaled <- t(scale(t(X_pca), center = TRUE, scale = TRUE))
  X_scaled[is.na(X_scaled)] <- 0
  pca <- prcomp(t(X_scaled), center = FALSE, scale. = FALSE)
  
  # Correlation between PC1 and response
  valid <- !is.na(resp_vec)
  if (sum(valid) < 10) {
    cat(sprintf("  %s: too few valid samples\n", label))
    return(invisible(NULL))
  }
  cor_val <- cor(pca$x[valid, 1], resp_vec[valid], use = "complete.obs")
  cat(sprintf("  %s: cor(PC1, response) = %.4f\n", label, cor_val))
  return(cor_val)
}

cat("  Pre-correction:\n")
cor_pre <- check_response_signal(X_combined, all_resp, batch, "Pre-ComBat")

if (!is.null(X_combat_old)) {
  cat("  Post OLD (NA=0.5):\n")
  cor_post_old <- check_response_signal(X_combat_old, resp_filled, batch, "Post-OLD")
}

if (!is.null(X_combat_new)) {
  cat("  Post NEW (drop NA):\n")
  cor_post_new <- check_response_signal(X_combat_new, all_resp_clean, batch_clean, "Post-NEW")
}

# ============================================================
# 7. Summary
# ============================================================
cat("\n=== ComBat v6.0 Summary ===\n")
cat(sprintf("  Total samples: %d\n", ncol(X_combined)))
cat(sprintf("  NA samples removed: %d\n", n_na))
cat(sprintf("  Clean samples: %d\n", sum(complete_idx)))

cat("\n  Batch separability (LDA accuracy, lower = better):\n")
cat(sprintf("    Pre-ComBat:        %.1f%%\n", lda_pre * 100))
if (!is.null(X_combat_old)) {
  cat(sprintf("    Post-OLD (NA=0.5): %.1f%%\n", lda_post_old * 100))
}
if (!is.null(X_combat_new)) {
  cat(sprintf("    Post-NEW (drop NA):%.1f%%\n", lda_post_new * 100))
}

cat("\n  Response signal (cor with PC1, higher = more biology preserved):\n")
cat(sprintf("    Pre-ComBat:        %.4f\n", cor_pre))
if (!is.null(X_combat_old)) {
  cat(sprintf("    Post-OLD (NA=0.5): %.4f\n", cor_post_old))
}
if (!is.null(X_combat_new)) {
  cat(sprintf("    Post-NEW (drop NA):%.4f\n", cor_post_new))
}

cat("\n=== Script 13 Complete ===\n")
