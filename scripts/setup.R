# ============================================================
# XELOX 耐药分子指纹研究 - 环境初始化脚本
# ============================================================
# 使用方法: Rscript scripts/setup.R
# 或在 R 中: source("scripts/setup.R")

# --- 设置 R 库路径 ---
lib_path <- "/path/to/Rlibs"
.libPaths(c(lib_path, .libPaths()))

# --- 设置临时目录（避免中文路径问题）---
Sys.setenv(TMPDIR = "/tmp", TMP = "/tmp", TEMP = "/tmp")

# --- 加载必需包 ---
required_pkgs <- c(
  "limma", "DESeq2", "WGCNA", "GEOquery", 
  "biomaRt", "preprocessCore", "impute", "sva"
)

cat("=== XELOX 耐药分子指纹 环境检查 ===\n")
for (pkg in required_pkgs) {
  status <- tryCatch({
    library(pkg, character.only = TRUE, quietly = TRUE)
    sprintf("  ✓ %s v%s", pkg, packageVersion(pkg))
  }, error = function(e) {
    sprintf("  ✗ %s - %s", pkg, conditionMessage(e))
  })
  cat(status, "\n")
}

# --- 项目根目录 ---
project_root <- normalizePath(file.path(getwd(), ".."))
cat(sprintf("\n项目根目录: %s\n", project_root))

# --- 数据目录 ---
data_geo_dir    <- file.path(project_root, "data", "geo")
data_tcga_dir   <- file.path(project_root, "data", "tcga")
data_proc_dir   <- file.path(project_root, "data", "processed")
results_fig_dir <- file.path(project_root, "results", "figures")
results_tab_dir <- file.path(project_root, "results", "tables")

cat("\n数据目录:\n")
cat(sprintf("  GEO:    %s\n", data_geo_dir))
cat(sprintf("  TCGA:   %s\n", data_tcga_dir))
cat(sprintf("  处理后: %s\n", data_proc_dir))
cat(sprintf("  图表:   %s\n", results_fig_dir))
cat(sprintf("  表格:   %s\n", results_tab_dir))

cat("\n=== 环境初始化完成 ===\n")
