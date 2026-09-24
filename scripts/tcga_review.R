# ============================================================
# TCGA分析审查：患者分组 + DEGs + 批次效应
# ============================================================
Sys.setenv(HOME = "/tmp", TMPDIR = "/tmp", R_LIBS = "/path/to/Rlibs")
.libPaths(c("/path/to/Rlibs", .libPaths()))

PROJ <- "/path/to/xelox_project基于可解释性机器学习的XELOX耐药分子指纹研究"

# ============================================================
# 任务1: 患者分组审查
# ============================================================
cat("========================================\n")
cat("任务1: TCGA患者分组审查\n")
cat("========================================\n\n")

resp <- read.csv(file.path(PROJ, "data/tcga/TCGA_COAD_READ_xelox_response.csv"), 
                 stringsAsFactors = FALSE)
cat(sprintf("总患者数: %d\n\n", nrow(resp)))

# 分组统计
tb <- table(resp$response_group)
cat("--- response_group 分布 ---\n")
print(tb)
cat("\n")

pct <- round(prop.table(tb) * 100, 2)
cat("--- 比例 (%) ---\n")
print(pct)
cat("\n")

cat(sprintf("Unknown比例: %.1f%%\n", 
    sum(resp$response_group == "unknown") / nrow(resp) * 100))
cat(sprintf("可分析: sensitive=%d, resistant=%d\n", 
    sum(resp$response_group == "sensitive"), 
    sum(resp$response_group == "resistant")))
cat(sprintf("mixed (排除): %d\n\n", 
    sum(resp$response_group == "mixed")))

# 癌症类型交叉
origins <- read.csv(file.path(PROJ, "data/tcga/TCGA_COAD_READ_sample_origins.csv"), 
                    stringsAsFactors = FALSE)
origins$patient_barcode <- substr(origins$sample_barcode, 1, 12)
origins <- origins[!duplicated(origins$patient_barcode), ]

cat(sprintf("患者来源(去重): COAD=%d, READ=%d\n\n", 
    sum(origins$cancer_type == "COAD"), 
    sum(origins$cancer_type == "READ")))

merged <- merge(resp, origins[, c("patient_barcode", "cancer_type")], 
                by = "patient_barcode", all.x = TRUE)
cat("--- 按癌症类型×response 交叉表 ---\n")
print(table(merged$cancer_type, merged$response_group))
cat("\n")

# XELOX用药确认
drug <- read.csv(file.path(PROJ, "data/tcga/TCGA_COAD_READ_drug_data.csv"), 
                 stringsAsFactors = FALSE)
cat(sprintf("药物记录总数: %d\n", nrow(drug)))
cat(sprintf("含奥沙利铂(Oxali)记录: %d\n", sum(drug$is_oxali == TRUE, na.rm = TRUE)))
cat(sprintf("含氟尿嘧啶(Fluoro)记录: %d\n", sum(drug$is_fluoro == TRUE, na.rm = TRUE)))
xelox_pts <- unique(drug$bcr_patient_barcode[drug$is_oxali == TRUE & drug$is_fluoro == TRUE])
cat(sprintf("同时使用奥沙利铂+氟尿嘧啶的患者数: %d\n", length(xelox_pts)))
cat("\n")

# ============================================================
# 任务2: TCGA DEGs审查
# ============================================================
cat("========================================\n")
cat("任务2: TCGA DESeq2 DEGs审查\n")
cat("========================================\n\n")

deg <- read.csv(file.path(PROJ, "results/tables/TCGA_DESeq2_results.csv"), 
                stringsAsFactors = FALSE)
cat(sprintf("总基因数: %d\n", nrow(deg)))
cat(sprintf("列名: %s\n\n", paste(colnames(deg), collapse = ", ")))

# 显著DEG (padj < 0.05)
sig <- deg[deg$padj < 0.05, ]
cat(sprintf("显著DEG (padj<0.05): %d\n", nrow(sig)))

sig_up <- sig[sig$log2FoldChange > 0, ]
sig_down <- sig[sig$log2FoldChange < 0, ]
cat(sprintf("  上调 (log2FC>0): %d\n", nrow(sig_up)))
cat(sprintf("  下调 (log2FC<0): %d\n", nrow(sig_down)))

# 严格DEG (padj<0.05 & |log2FC|>1)
sig_strict <- sig[abs(sig$log2FoldChange) > 1, ]
cat(sprintf("\n严格DEG (padj<0.05 & |log2FC|>1): %d\n", nrow(sig_strict)))
cat(sprintf("  上调: %d\n", sum(sig_strict$log2FoldChange > 0)))
cat(sprintf("  下调: %d\n", sum(sig_strict$log2FoldChange < 0)))

# 极显著DEG (padj<1e-10)
sig_top <- sig[sig$padj < 1e-10, ]
cat(sprintf("\n极显著DEG (padj<1e-10): %d\n", nrow(sig_top)))

# Top10
cat("\n--- Top10 DEG (按padj排序) ---\n")
top10 <- sig[order(sig$padj), ][1:10, ]
for (i in 1:nrow(top10)) {
  direction <- ifelse(top10$log2FoldChange[i] > 0, "UP", "DOWN")
  cat(sprintf("  %s | log2FC=%.2f | pval=%.2e | padj=%.2e | %s\n",
      top10$gene[i], top10$log2FoldChange[i], top10$pvalue[i], 
      top10$padj[i], direction))
}

# ENSG版本号去除
deg$gene_clean <- gsub("\\..*$", "", deg$gene)
cat(sprintf("\n唯一ENSG ID数: %d\n", length(unique(deg$gene_clean))))

# 保存严格DEG列表用于后续分析
strict_out <- file.path(PROJ, "results/tables/TCGA_strict_DEGs.csv")
write.csv(sig_strict, strict_out, row.names = FALSE)
cat(sprintf("\n严格DEG已保存: %s\n", strict_out))

# ============================================================
# 任务3: COAD/READ批次效应初步评估
# ============================================================
cat("\n========================================\n")
cat("任务3: COAD/READ批次效应评估\n")
cat("========================================\n\n")

# 加载gene_pool_expression (已整合的表达矩阵)
gp_expr <- read.csv(file.path(PROJ, "data/tcga/TCGA_COAD_READ_gene_pool_expression.csv"), 
                    stringsAsFactors = FALSE)
cat(sprintf("基因池表达矩阵: %d 基因 × %d 样本\n", nrow(gp_expr), ncol(gp_expr)))

# 提取样本barcode
sample_cols <- colnames(gp_expr)[-1]  # 第一列是gene

# 获取样本的cancer_type
sample_info <- read.csv(file.path(PROJ, "data/tcga/TCGA_COAD_READ_sample_origins.csv"), 
                        stringsAsFactors = FALSE)

# 匹配样本
# 样本barcode格式: TCGA-XX-XXXX-XX...
matched <- sample_info[sample_info$sample_barcode %in% sample_cols, ]
cat(sprintf("匹配的样本: COAD=%d, READ=%d\n", 
    sum(matched$cancer_type == "COAD"), 
    sum(matched$cancer_type == "READ")))

# 简单PCA评估批次效应
gene_mat <- as.matrix(gp_expr[, -1])
rownames(gene_mat) <- gp_expr[, 1]

# 转置: 样本为行
expr_t <- t(gene_mat)
expr_t <- expr_t[, colSums(is.na(expr_t)) == 0]  # 去掉NA列

# PCA
pca <- prcomp(expr_t, scale. = TRUE, center = TRUE)
pca_var <- summary(pca)$importance[2, 1:3] * 100

cat(sprintf("\nPCA: PC1=%.1f%%, PC2=%.1f%%, PC3=%.1f%%\n", 
    pca_var[1], pca_var[2], pca_var[3]))

# 获取批次信息
batch_info <- data.frame(
  sample = rownames(expr_t),
  cancer_type = matched$cancer_type[match(rownames(expr_t), matched$sample_barcode)],
  stringsAsFactors = FALSE
)
batch_info$cancer_type[is.na(batch_info$cancer_type)] <- "UNKNOWN"

# PC1/PC2 by cancer type
pca_data <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  cancer_type = batch_info$cancer_type
)

cat("\n--- 按cancer_type的PC1描述性统计 ---\n")
for (ct in unique(pca_data$cancer_type)) {
  sub <- pca_data$PC1[pca_data$cancer_type == ct]
  cat(sprintf("  %s: mean=%.2f, sd=%.2f, n=%d\n", ct, mean(sub), sd(sub), length(sub)))
}

# 如果R包sva可用，尝试ComBat
if (requireNamespace("sva", quietly = TRUE)) {
  library(sva)
  
  # 准备ComBat输入
  # 需要: 完整表达矩阵, 批次向量
  combat_expr <- gene_mat[, colnames(gene_mat) %in% matched$sample_barcode]
  combat_samples <- colnames(combat_expr)
  combat_batch <- matched$cancer_type[match(combat_samples, matched$sample_barcode)]
  
  # 检查每组的样本数
  batch_tab <- table(combat_batch)
  cat(sprintf("\nComBat批次信息: %s\n", 
      paste(names(batch_tab), batch_tab, sep = "=", collapse = ", ")))
  
  if (length(unique(combat_batch)) > 1 && min(batch_tab) >= 3) {
    cat("\n运行 sva::ComBat() 批次校正...\n")
    
    # 检查是否有零方差基因
    gene_vars <- apply(combat_expr, 1, var, na.rm = TRUE)
    combat_expr_filt <- combat_expr[gene_vars > 0, ]
    cat(sprintf("  零方差基因移除后: %d 基因\n", nrow(combat_expr_filt)))
    
    # 执行ComBat
    combat_corrected <- tryCatch({
      ComBat(dat = as.matrix(combat_expr_filt), 
             batch = combat_batch,
             mod = NULL,
             par.prior = TRUE,
             prior.plots = FALSE)
    }, error = function(e) {
      cat(sprintf("  ComBat出错: %s\n", conditionMessage(e)))
      cat("  尝试使用par.prior=FALSE...\n")
      ComBat(dat = as.matrix(combat_expr_filt), 
             batch = combat_batch,
             mod = NULL,
             par.prior = FALSE,
             prior.plots = FALSE)
    })
    
    # 校正后PCA
    combat_t <- t(combat_corrected)
    pca2 <- prcomp(combat_t, scale. = TRUE, center = TRUE)
    pca2_var <- summary(pca2)$importance[2, 1:3] * 100
    
    cat(sprintf("校正后PCA: PC1=%.1f%%, PC2=%.1f%%, PC3=%.1f%%\n", 
        pca2_var[1], pca2_var[2], pca2_var[3]))
    
    # 保存校正后表达矩阵
    combat_out <- file.path(PROJ, "data/tcga/TCGA_COAD_READ_expr_combat.csv")
    combat_df <- data.frame(gene = rownames(combat_corrected), combat_corrected, 
                            stringsAsFactors = FALSE, check.names = FALSE)
    write.csv(combat_df, combat_out, row.names = FALSE)
    cat(sprintf("校正后表达矩阵已保存: %s\n\n", combat_out))
    
  } else {
    cat("  SKIP: 批次效应校正条件不满足（需至少2个批次，每批>=3样本）\n\n")
  }
} else {
  cat("\n  sva包未安装，跳过ComBat批次校正\n\n")
}

cat("\n========================================\n")
cat("TCGA分析审查完成\n")
cat("========================================\n")
