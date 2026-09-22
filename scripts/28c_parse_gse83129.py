#!/usr/bin/env python3
"""
28c_parse_gse83129.py
直接解析 GSE83129 系列矩阵文件，提取表达矩阵和临床数据
"""
import gzip
import re
import pandas as pd
import numpy as np
import os

GEO_DIR = "C:/xelox_research/data/geo"
GZ_FILE = os.path.join(GEO_DIR, "GSE83129_series_matrix.txt.gz")

# ============================================================
# 1. 解析系列矩阵头部的临床信息
# ============================================================
print("=" * 60)
print("解析临床信息...")
print("=" * 60)

# 读取所有行
with gzip.open(GZ_FILE, "rt") as f:
    lines = f.readlines()

# 分离头部和数据
header_lines = []
data_start = 0
for i, line in enumerate(lines):
    if line.startswith("!series_matrix_table_begin"):
        data_start = i + 1
        break
    header_lines.append(line)

# 解析样本级信息（以 !Sample_ 开头的行）
# 这些行格式为：!Sample_xxx\tval1\tval2\t...
sample_data = {}
for line in header_lines:
    if line.startswith("!Sample_"):
        parts = line.strip().split("\t")
        key = parts[0].replace("!Sample_", "")
        values = [p.strip('"') for p in parts[1:]]
        sample_data[key] = values

sample_names = sample_data.get("geo_accession", [])
print(f"样本数: {len(sample_names)}")
print(f"样本ID: {sample_names}")

# 提取样本标题（包含组织类型信息）
titles = sample_data.get("title", [])
print(f"\n样本标题: {titles}")

# 提取所有临床特征
print("\n--- 所有样本级字段 ---")
for key, values in sample_data.items():
    if len(values) < 30:
        continue  # 跳过非样本级字段
    unique_vals = set(values)
    print(f"\n{key} (N={len(unique_vals)} unique):")
    for v in sorted(unique_vals)[:20]:
        print(f"  {v}")

# ============================================================
# 2. 提取表达矩阵
# ============================================================
print("\n" + "=" * 60)
print("提取表达矩阵...")
print("=" * 60)

# 读取数据行
data_lines = []
for line in lines[data_start:]:
    stripped = line.strip()
    if stripped == "!series_matrix_table_end":
        break
    data_lines.append(stripped)

print(f"数据行数: {len(data_lines)}")

# 解析第一行（表头）
header = data_lines[0].split("\t")
col_names = [h.strip('"') for h in header]
print(f"列名 (前5): {col_names[:5]}")
print(f"总列数: {len(col_names)}")

# 解析数据
data_rows = []
probe_ids = []
for line in data_lines[1:]:
    parts = line.split("\t")
    probe_ids.append(parts[0].strip('"'))
    data_rows.append([float(x) for x in parts[1:]])

expr_df = pd.DataFrame(data_rows, index=probe_ids, columns=col_names[1:])
print(f"\n表达矩阵: {expr_df.shape[0]} probes x {expr_df.shape[1]} samples")

# 保存
expr_df.to_pickle(os.path.join(GEO_DIR, "GSE83129_expression.pkl"))
print(f"已保存: {os.path.join(GEO_DIR, 'GSE83129_expression.pkl')}")

# 保存临床信息
# 样本标题中包含组织类型信息（"M226-01" 等，有些带 "L" 后缀可能表示正常/肿瘤？）
# 需要从标题中的 "L" 后缀识别正常 vs 肿瘤组织
sample_info = pd.DataFrame({
    "geo_accession": sample_names,
    "title": titles
})

# 根据样本名称中的 "L" 判断组织类型
# 从论文来看，"L" 可能代表 "normal epithelium" 激光显微切割的正常组织
is_normal = [t.endswith("L") for t in titles]
sample_info["tissue_type"] = ["normal" if n else "tumor" for n in is_normal]

print(f"\n组织类型分布:")
print(sample_info["tissue_type"].value_counts())

# 检查样本标题中是否包含响应信息
# 从论文 miR-625-3p 来看，26 例患者应该是 mCRC 肿瘤组织样本
# 其中 9 NR (non-responder) / 17 R (responder)
# 不对应信息不在系列矩阵中，需要从论文补充

sample_info.to_csv(os.path.join(GEO_DIR, "GSE83129_sample_info.csv"), index=False)
print(f"已保存: {os.path.join(GEO_DIR, 'GSE83129_sample_info.csv')}")

print("\nDone.")
