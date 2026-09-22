"""Convert RDS files to CSV using pyreadr"""
import pyreadr
import pandas as pd
import os

PROJECT = r"/path/to/xelox_project"

files = [
    (os.path.join(PROJECT, "results", "tables", "pathway_activity", "GSE39582_pathway_scores.rds"),
     os.path.join(PROJECT, "results", "tables", "pathway_activity", "GSE39582_pathway_scores.csv")),
    (os.path.join(PROJECT, "results", "tables", "pathway_activity", "GSE39582_gene_expression.rds"),
     os.path.join(PROJECT, "results", "tables", "pathway_activity", "GSE39582_gene_expression.csv")),
    (os.path.join(PROJECT, "data", "processed", "GSE39582_clinical_processed.rds"),
     os.path.join(PROJECT, "data", "processed", "GSE39582_clinical_processed.csv")),
]

for rds_path, csv_path in files:
    print(f"Converting: {os.path.basename(rds_path)}")
    result = pyreadr.read_r(rds_path)
    if isinstance(result, dict):
        df = list(result.values())[0]
    else:
        df = result
    # Transpose if it looks like a matrix (pathway scores / expression)
    if isinstance(df, pd.DataFrame):
        print(f"  Shape: {df.shape}, cols sample: {list(df.columns[:3])}")
    df.to_csv(csv_path)
    print(f"  -> Saved: {csv_path}")

print("\nDone!")
