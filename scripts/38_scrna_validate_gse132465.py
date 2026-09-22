"""Script 38: GSE132465 scRNA-seq Cell Type Validation
==========
XELOX resistance gene validation in CRC scRNA-seq data.

Pipeline:
  1. Load UMI count matrix + cell annotations
  2. Create AnnData object
  3. QC filtering
  4. Normalization, HVG selection, PCA, UMAP
  5. Cell type annotation (using provided labels)
  6. XELOX fingerprint gene expression by cell type
  7. Capecitabine metabolism enzyme cell-type specificity
  8. Save results and generate figures
"""

import os, sys, warnings, gzip, zlib
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import scanpy as sc
from scipy.sparse import csr_matrix
warnings.filterwarnings('ignore')

PROJECT = r"/path/to/xelox_project"
DATA_DIR = os.path.join(PROJECT, "data", "scrna")
FIG_DIR = os.path.join(PROJECT, "results", "figures", "scrna_seq")
os.makedirs(FIG_DIR, exist_ok=True)

print("=" * 70)
print("GSE132465 CRC scRNA-seq — XELOX Resistance Gene Validation")
print("=" * 70)

# ============================================================
# 1. Load data
# ============================================================
print("\n[1] Loading data...")

# Load cell annotations
ann_path = os.path.join(DATA_DIR, "GSE132465_cell_annotation.txt.gz")
print(f"  Reading cell annotations: {ann_path}")
ann = pd.read_csv(ann_path, sep='\t', compression='gzip')
ann = ann.set_index('Index')
print(f"  Cells: {ann.shape[0]}")
print(f"  Cell types: {ann['Cell_type'].value_counts().to_dict()}")

# Load UMI count matrix (gzipped TSV, may be truncated)
umi_path = os.path.join(DATA_DIR, "GSE132465_raw_UMI_count_matrix.txt.gz")
print(f"  Reading UMI matrix (large file, may take a while)...")

# Use zlib for streaming decompress (handles truncated gzip)
with open(umi_path, 'rb') as f:
    compressed = f.read()
decompressor = zlib.decompressobj(15 + 32)
raw_data = decompressor.decompress(compressed)
raw_text = raw_data.decode('utf-8', errors='replace')
all_lines = raw_text.split('\n')

# Parse header
header = all_lines[0].strip().split('\t')
cell_barcodes = header[1:]  # first col is gene index
print(f"  Header: {len(cell_barcodes)} cells")

# Parse gene-by-gene and build matrix
n_cells = len(cell_barcodes)
n_genes = len(all_lines) - 1  # subtract header
print(f"  Genes to read: {n_genes}")

gene_names = []
umi_values = np.zeros((n_genes, n_cells), dtype=np.float32)

for i, line in enumerate(all_lines[1:]):
    if not line.strip():
        n_genes = i
        break
    parts = line.strip().split('\t')
    # Handle truncated rows - only take complete rows
    if len(parts) < 1 + n_cells:
        n_genes = i
        gene_names = gene_names[:i]
        umi_values = umi_values[:i]
        print(f"\n    Truncated at row {i+1}, using {i} complete genes")
        break
    gene_names.append(parts[0])
    vals = np.array([float(x) if x else 0.0 for x in parts[1:1+n_cells]], dtype=np.float32)
    umi_values[i] = vals
    if (i+1) % 5000 == 0:
        print(f"    {i+1}/{n_genes} genes parsed", end='\r')

# Trim
umi_values = umi_values[:len(gene_names)]
print(f"  Matrix: {len(gene_names)} genes x {n_cells} cells")

# Create AnnData - cell barcodes in matrix match annotation index
barcodes = pd.Index(cell_barcodes)
common_idx = barcodes.intersection(ann.index)
n_matched = len(common_idx)
mask = barcodes.isin(common_idx)
print(f"  Matching barcodes: {n_matched}/{len(cell_barcodes)}")

# Subset to matched cells and transpose for AnnData (obs x var)
matched_umi = csr_matrix(umi_values[:, mask].T)  # shape: (n_matched_cells, n_genes), sparse
barcodes_matched = barcodes[mask]
print(f"  AnnData input shape: {matched_umi.shape}")

adata = sc.AnnData(
    X=matched_umi,
    obs=ann.loc[barcodes_matched],
    var=pd.DataFrame(index=gene_names)
)
# Check matching
common_cells = adata.obs_names.intersection(ann.index)
print(f"  Matched cells with annotation: {len(common_cells)}")

# ============================================================
# 2. Basic QC
# ============================================================
print("\n[2] QC filtering...")
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

print(f"  Before QC: {adata.n_obs} cells")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.pct_counts_mt < 25, :].copy()
print(f"  After QC: {adata.n_obs} cells, {adata.n_vars} genes")

# QC figures
fig, axes = plt.subplots(2, 3, figsize=(14, 8))
axes[0,0].hist(adata.obs.n_genes_by_counts, bins=50, color='#378ADD')
axes[0,0].set_xlabel('Genes per cell'); axes[0,0].set_ylabel('Count')
axes[0,1].hist(adata.obs.total_counts, bins=50, color='#378ADD')
axes[0,1].set_xlabel('UMI counts per cell')
axes[0,2].hist(adata.obs.pct_counts_mt, bins=50, color='#E64B35')
axes[0,2].set_xlabel('% mitochondrial genes')
sc.pl.violin(adata, 'n_genes_by_counts', groupby='Cell_type', ax=axes[1,0], show=False)
sc.pl.violin(adata, 'total_counts', groupby='Cell_type', ax=axes[1,1], show=False)
sc.pl.violin(adata, 'pct_counts_mt', groupby='Cell_type', ax=axes[1,2], show=False)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'QC_metrics.png'), dpi=150, bbox_inches='tight')
plt.close()
print("  Saved QC_metrics.png")

# ============================================================
# 3. Normalization & Dimensionality Reduction
# ============================================================
print("\n[3] Normalization & clustering...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# HVG
sc.pp.highly_variable_genes(adata, n_top_genes=3000)
n_hvg = adata.var.highly_variable.sum()
print(f"  HVG: {n_hvg}")

# Scale + PCA
adata_hvg = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=50, svd_solver='arpack')
print(f"  PCA done: {adata_hvg.obsm['X_pca'].shape}")

# Neighbors + UMAP
sc.pp.neighbors(adata_hvg, n_pcs=30, n_neighbors=15)
sc.tl.umap(adata_hvg, min_dist=0.3)
print("  UMAP done")

# Copy to main adata
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
adata.obsm['X_umap'] = adata_hvg.obsm['X_umap']
adata.uns['neighbors'] = adata_hvg.uns['neighbors']
adata.obsp['distances'] = adata_hvg.obsp['distances']
adata.obsp['connectivities'] = adata_hvg.obsp['connectivities']
sc.tl.leiden(adata, resolution=0.5)

# ============================================================
# 4. UMAP Visualization
# ============================================================
print("\n[4] Visualization...")

# Cell type UMAP
fig, ax = plt.subplots(figsize=(9, 7))
cell_types = adata.obs['Cell_type'].unique()
colors = plt.cm.tab20(np.linspace(0, 1, len(cell_types)))
for ct, c in zip(cell_types, colors):
    mask = adata.obs['Cell_type'] == ct
    ax.scatter(adata.obsm['X_umap'][mask, 0], adata.obsm['X_umap'][mask, 1],
              s=2, c=[c], label=ct, alpha=0.5)
ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
ax.set_title('GSE132465 CRC: Cell Types (provided annotation)', fontweight='bold')
ax.legend(fontsize=7, markerscale=5, loc='best')
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'UMAP_cell_types.png'), dpi=200, bbox_inches='tight')
plt.close()
print("  Saved UMAP_cell_types.png")

# Tumor vs Normal UMAP
fig, ax = plt.subplots(figsize=(9, 7))
for cls, c in zip(['Tumor', 'Normal'], ['#E64B35', '#4DBBD5']):
    mask = adata.obs['Class'] == cls
    n = mask.sum()
    ax.scatter(adata.obsm['X_umap'][mask, 0], adata.obsm['X_umap'][mask, 1],
              s=2, c=c, label=f'{cls} (n={n})', alpha=0.4)
ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
ax.set_title('GSE132465: Tumor vs Normal Tissue', fontweight='bold')
ax.legend(fontsize=9, markerscale=5)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'UMAP_tumor_vs_normal.png'), dpi=200, bbox_inches='tight')
plt.close()
print("  Saved UMAP_tumor_vs_normal.png")

# Cell subtype UMAP (finer classification)
fig, ax = plt.subplots(figsize=(10, 8))
subtypes = adata.obs['Cell_subtype'].unique()
colors2 = plt.cm.tab20(np.linspace(0, 1, len(subtypes)))
for st, c in zip(subtypes, colors2):
    mask = adata.obs['Cell_subtype'] == st
    n = mask.sum()
    ax.scatter(adata.obsm['X_umap'][mask, 0], adata.obsm['X_umap'][mask, 1],
              s=2, c=[c], label=f'{st} (n={n})', alpha=0.5)
ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
ax.set_title('GSE132465: Cell Subtypes', fontweight='bold')
ax.legend(fontsize=6, markerscale=4, loc='best')
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'UMAP_cell_subtypes.png'), dpi=200, bbox_inches='tight')
plt.close()
print("  Saved UMAP_cell_subtypes.png")

# ============================================================
# 5. XELOX Resistance Gene Expression
# ============================================================
print("\n[5] XELOX resistance gene analysis...")

# Load gene lists
meta_df = pd.read_csv(os.path.join(PROJECT, "results/tables/drug_prediction/253_meta_genes.csv"))
genes_up = set(meta_df[meta_df['direction_consensus']=='up']['gene_symbol'])
genes_dn = set(meta_df[meta_df['direction_consensus']=='down']['gene_symbol'])

fingerprint = ['CSNK1G2', 'C5ORF42', 'KAZN', 'KLK6', 'MGA', 'MID2', 'SOX11', 'TAS2R40', 'ZNF451']
cape_enz = ['CES1', 'CES2', 'TYMS', 'DPYD', 'DCK', 'CDA', 'TP']

# Find available genes
fp_avail = [g for g in fingerprint if g in adata.var_names]
meta_up_avail = [g for g in genes_up if g in adata.var_names]
meta_dn_avail = [g for g in genes_dn if g in adata.var_names]
cape_avail = [g for g in cape_enz if g in adata.var_names]

print(f"  Fingerprint genes: {len(fp_avail)}/{len(fingerprint)} available")
print(f"  Meta up-genes: {len(meta_up_avail)}/{len(genes_up)} available")
print(f"  Meta down-genes: {len(meta_dn_avail)}/{len(genes_dn)} available")
print(f"  Capecitabine enzymes: {len(cape_avail)}/{len(cape_enz)} available")

if fp_avail:
    print(f"    Available: {', '.join(fp_avail)}")

# Compute per-cell expression scores
adata.obs['fp_score'] = adata[:, [g for g in fp_avail if g in adata.var_names]].X.toarray().mean(axis=1) if fp_avail else 0
adata.obs['meta_up_score'] = adata[:, [g for g in meta_up_avail if g in adata.var_names]].X.toarray().mean(axis=1) if meta_up_avail else 0
adata.obs['meta_dn_score'] = adata[:, [g for g in meta_dn_avail if g in adata.var_names]].X.toarray().mean(axis=1) if meta_dn_avail else 0
adata.obs['cape_enz_score'] = adata[:, [g for g in cape_avail if g in adata.var_names]].X.toarray().mean(axis=1) if cape_avail else 0

# Figure: Fingerprint gene expression on UMAP
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
score_cols = ['fp_score', 'meta_up_score', 'meta_dn_score', 'cape_enz_score']
titles = ['Fingerprint Genes (avg log-expr)',
          'Meta Up-in-Resistant Genes (avg log-expr)',
          'Meta Down-in-Resistant Genes (avg log-expr)',
          'Capecitabine Metabolism Enzymes (avg log-expr)']

for ax, col, title in zip(axes.flat[:4], score_cols, titles):
    scatter = ax.scatter(adata.obsm['X_umap'][:, 0], adata.obsm['X_umap'][:, 1],
                   c=adata.obs[col], s=1, cmap='Reds', alpha=0.5)
    ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
    ax.set_title(title, fontsize=10, fontweight='bold')
    plt.colorbar(scatter, ax=ax, shrink=0.6)

# Highlight tumor vs normal
for cls, c in zip(['Tumor', 'Normal'], ['#E64B35', '#4DBBD5']):
    mask = adata.obs['Class'] == cls
    axes.flat[4].scatter(adata.obsm['X_umap'][mask, 0], adata.obsm['X_umap'][mask, 1],
                        s=2, c=c, label=cls, alpha=0.3)
axes.flat[4].set_title('Tumor vs Normal', fontweight='bold')
axes.flat[4].legend(fontsize=8, markerscale=5)
axes.flat[4].set_xlabel('UMAP1'); axes.flat[4].set_ylabel('UMAP2')

axes.flat[5].axis('off')
plt.suptitle('XELOX Resistance Gene Expression in CRC scRNA-seq', fontsize=13, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'resistance_genes_UMAP.png'), dpi=200, bbox_inches='tight')
plt.close()
print("  Saved resistance_genes_UMAP.png")

# ============================================================
# 6. Cell-type Specificity Analysis
# ============================================================
print("\n[6] Cell-type specificity analysis...")

# Expression by cell type
ct_order = adata.obs.groupby('Cell_type').size().sort_values(ascending=False).index

# Bar plot: fingerprint score by cell type
fig, ax = plt.subplots(figsize=(10, 5))
ct_scores = adata.obs.groupby('Cell_type')[['fp_score', 'meta_up_score', 'meta_dn_score']].mean()
ct_scores.loc[ct_order].plot(kind='bar', ax=ax, color=['#7F77DD', '#E64B35', '#4DBBD5'])
ax.set_ylabel('Mean log-expression')
ax.set_title('XELOX Resistance Gene Scores by Cell Type', fontweight='bold')
ax.legend(fontsize=8)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'gene_scores_by_celltype.png'), dpi=200, bbox_inches='tight')
plt.close()
print("  Saved gene_scores_by_celltype.png")

# Capecitabine enzymes by cell type
if cape_avail:
    fig, ax = plt.subplots(figsize=(10, 5))
    cape_ct = adata.obs.groupby('Cell_type')['cape_enz_score'].mean()
    cape_ct.loc[ct_order].plot(kind='bar', ax=ax, color='#378ADD')
    ax.set_ylabel('Mean log-expression')
    ax.set_title('Capecitabine Metabolism Enzyme Expression by Cell Type', fontweight='bold')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'cape_enz_by_celltype.png'), dpi=200, bbox_inches='tight')
    plt.close()
    print("  Saved cape_enz_by_celltype.png")

# Per-gene Dotplot for key genes
all_key_genes = fp_avail + meta_up_avail[:5] + meta_dn_avail[:5] + cape_avail
all_key_genes = list(dict.fromkeys(all_key_genes))  # deduplicate preserving order
all_key_genes = [g for g in all_key_genes if g in adata.var_names]

if all_key_genes:
    import scanpy.pl as splt
    try:
        splt.dotplot(adata, var_names=all_key_genes[:20], groupby='Cell_type',
                     save='_xelox_key_genes.png', figsize=(8, 5),
                     dendrogram=False)
        print("  Saved dotplot_xelox_key_genes.png")
    except Exception as e:
        print(f"  Dotplot failed (non-critical): {e}")
        # Fallback: create simple heatmap
        fig, ax = plt.subplots(figsize=(8, 4))
        ct_order = adata.obs.groupby('Cell_type').size().sort_values(ascending=False).index
        plot_data = pd.DataFrame(
            {g: adata[:, g].X.toarray().mean(axis=0) for g in all_key_genes[:15]},
            index=adata.obs_names
        )
        sns.heatmap(plot_data.groupby(adata.obs['Cell_type']).mean().T,
                   cmap='viridis', ax=ax, cbar_kws={'label': 'Mean log-expr'})
        ax.set_title('Key Gene Expression by Cell Type', fontweight='bold')
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, 'key_gene_heatmap.png'), dpi=200, bbox_inches='tight')
        plt.close()
        print("  Saved key_gene_heatmap.png (fallback)")

# ============================================================
# 7. Tumor-specific Analysis
# ============================================================
print("\n[7] Tumor-specific analysis...")
tumor_adata = adata[adata.obs['Class'] == 'Tumor'].copy()
if tumor_adata.n_obs > 0:
    print(f"  Tumor cells: {tumor_adata.n_obs}")
    
    # Tumor epithelial vs stromal fingerprint expression
    tumor_ct = tumor_adata.obs.groupby('Cell_type').size().sort_values(ascending=False)
    print(f"  Tumor cell types: {tumor_ct.to_dict()}")
    
    # Tumor vs normal comparison
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for idx, (col, title) in enumerate([('fp_score', 'Fingerprint Score'), ('cape_enz_score', 'Capecitabine Enz Score')]):
        tumor_vals = adata[adata.obs['Class']=='Tumor'].obs[col]
        normal_vals = adata[adata.obs['Class']=='Normal'].obs[col]
        axes[idx].boxplot([tumor_vals, normal_vals], labels=['Tumor', 'Normal'], 
                          patch_artist=True,
                          boxprops=dict(facecolor='#E64B35' if idx==0 else '#378ADD', alpha=0.6),
                          medianprops=dict(color='black'))
        axes[idx].set_ylabel('Score')
        axes[idx].set_title(title, fontweight='bold')
    plt.suptitle('Tumor vs Normal: Resistance Gene Expression', fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'tumor_vs_normal_scores.png'), dpi=200, bbox_inches='tight')
    plt.close()
    print("  Saved tumor_vs_normal_scores.png")

# ============================================================
# 8. Summary Statistics
# ============================================================
print("\n[8] Summary...")

# Cell type counts
n_tumor = (adata.obs['Class'] == 'Tumor').sum()
n_normal = (adata.obs['Class'] == 'Normal').sum()
top_ct = adata.obs['Cell_type'].value_counts().head(5)

print(f"""
=== GSE132465 scRNA-seq Summary ===
Dataset: 23 CRC patients + 10 normal mucosa
Total cells after QC: {adata.n_obs}
  Tumor: {n_tumor} ({n_tumor/adata.n_obs*100:.0f}%)
  Normal: {n_normal} ({n_normal/adata.n_obs*100:.0f}%)

Top 5 cell types:
{top_ct.to_string()}

XELOX resistance genes:
  Fingerprint: {len(fp_avail)}/{len(fingerprint)} found
  Meta up-in-res: {len(meta_up_avail)}/{len(genes_up)} found
  Meta down-in-res: {len(meta_dn_avail)}/{len(genes_dn)} found
  Capecitabine enz: {len(cape_avail)}/{len(cape_enz)} found

Top fingerprint-expressing cell types:
{adata.obs.groupby('Cell_type')['fp_score'].mean().sort_values(ascending=False).head(3).to_string()}

Highest capecitabine enzyme-expressing cell types:
{adata.obs.groupby('Cell_type')['cape_enz_score'].mean().sort_values(ascending=False).head(3).to_string()}
""")

# Save
adata.write(os.path.join(FIG_DIR, 'GSE132465_annotated.h5ad'))
print(f"  Saved annotated h5ad")
print(f"  Figures: {FIG_DIR}")
print("\nDONE!")
