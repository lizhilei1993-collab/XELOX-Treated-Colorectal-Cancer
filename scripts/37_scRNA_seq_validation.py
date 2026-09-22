"""Script 37: scRNA-seq Analysis for XELOX Resistance Validation
==========
Target: Validate 10-gene fingerprint + 253 meta-genes in scRNA-seq data.

Approach:
  1. Load CRC scRNA-seq dataset (GSE132465 or other available)
  2. Basic QC, normalization, clustering (scanpy)
  3. Annotate cell types (epithelial, CAF, immune, etc.)
  4. Check 10-fingerprint genes and 253-genes distribution across cell types
  5. Check capecitabine activation enzymes (CES2, TP, DPYD) cell-type specificity

NOTE: This is a template script. To run:
  1. Download scRNA-seq data (e.g., from GEO or Zenodo)
  2. Update DATA_PATH below
  3. Run: python scripts/37_scRNA_seq_validation.py

Prerequisites: pip install scanpy leidenalg
"""

import os, sys, warnings
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
warnings.filterwarnings('ignore')

PROJECT = r"/path/to/xelox_project"
TAB_DIR = os.path.join(PROJECT, "results", "tables", "drug_prediction")
FIG_DIR = os.path.join(PROJECT, "results", "figures", "scrna_seq")
os.makedirs(FIG_DIR, exist_ok=True)

print("=" * 70)
print("XELOX RESISTANCE — scRNA-seq VALIDATION")
print("=" * 70)

# ============================================================
# DATA: Need to download first
# ============================================================
# Option 1: CRC scRNA-seq from GSE132465 (treatment-naive CRC, 23 patients)
#   Download: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465
#   Or from Zenodo: https://zenodo.org/records/3967538
#
# Option 2: CRC scRNA-seq from GSE146771 (another CRC atlas)
#
# Option 3: Oxaliplatin-resistant organoid data (if available)
#   Paper: Cellular Oncology 2022 (PMID: 36136268) 
#   Check Supplementary Materials for data access

# ============================================================
# ANALYSIS PIPELINE (TO BE RUN AFTER DATA DOWNLOAD)
# ============================================================
def run_scrna_pipeline(h5ad_path, output_prefix="CRC"):
    """
    Complete scRNA-seq analysis pipeline.
    
    Parameters
    ----------
    h5ad_path : str
        Path to scanpy h5ad file or 10x count matrix directory
    output_prefix : str
        Prefix for output files
    """
    
    print(f"\n[1] Loading data from {h5ad_path}...")
    adata = sc.read_h5ad(h5ad_path)
    print(f"  Cells: {adata.n_obs}, Genes: {adata.n_vars}")
    
    # QC
    print("\n[2] QC filtering...")
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Filter
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs.pct_counts_mt < 20, :].copy()
    print(f"  After QC: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # Normalize
    print("\n[3] Normalization...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # HVG
    sc.pp.highly_variable_genes(adata, n_top_genes=3000)
    adata_hvg = adata[:, adata.var.highly_variable].copy()
    
    # Scale + PCA
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=50, svd_solver='arpack')
    
    # Neighbors + UMAP
    sc.pp.neighbors(adata_hvg, n_pcs=30)
    sc.tl.umap(adata_hvg)
    
    # Clustering
    sc.tl.leiden(adata_hvg, resolution=0.8)
    
    # Copy annotations back
    adata.obs['leiden'] = adata_hvg.obs['leiden']
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
    adata.obsm['X_umap'] = adata_hvg.obsm['X_umap']
    adata.obs['n_genes_by_counts'] = adata_hvg.obs['n_genes_by_counts']
    
    # ============================================================
    # Cell Type Annotation
    # ============================================================
    print("\n[4] Cell type annotation...")
    
    # Known markers for CRC
    cell_markers = {
        'Epithelial': ['EPCAM', 'KRT19', 'KRT18', 'CDH1'],
        'Fibroblast_CAF': ['COL1A1', 'COL1A2', 'FAP', 'ACTA2', 'PDGFRA', 'PDGFRB'],
        'T_cell': ['CD3D', 'CD3E', 'CD8A', 'CD4'],
        'NK_cell': ['NKG7', 'GNLY', 'KLRD1'],
        'B_cell': ['MS4A1', 'CD79A', 'CD19'],
        'Myeloid_Macrophage': ['CD68', 'CD14', 'FCGR3A', 'CSF1R', 'ITGAM'],
        'Endothelial': ['PECAM1', 'CDH5', 'VWF', 'ENG'],
        'Mast_cell': ['KIT', 'TPSAB1', 'CPA3'],
        'Plasma_cell': ['MZB1', 'SDC1', 'JCHAIN', 'IGHG1'],
    }
    
    scores = {}
    for ct, genes in cell_markers.items():
        avail = [g for g in genes if g in adata.var_names]
        if avail:
            scores[ct] = adata[:, avail].X.toarray().mean(axis=1)
    
    score_df = pd.DataFrame(scores, index=adata.obs_names)
    predicted = score_df.idxmax(axis=1)
    adata.obs['cell_type'] = predicted.values
    
    print("  Cell type distribution:")
    print(adata.obs['cell_type'].value_counts())
    
    # Save UMAP plot
    fig, ax = plt.subplots(figsize=(8, 6))
    for ct in adata.obs['cell_type'].unique():
        mask = adata.obs['cell_type'] == ct
        ax.scatter(adata.obsm['X_umap'][mask, 0], adata.obsm['X_umap'][mask, 1],
                  s=3, label=ct, alpha=0.6)
    ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
    ax.set_title(f'{output_prefix}: Cell Types', fontweight='bold')
    ax.legend(fontsize=7, loc='best', markerscale=3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, f'{output_prefix}_cell_types.png'), dpi=200)
    plt.close()
    print(f"\n  Saved {output_prefix}_cell_types.png")
    
    # ============================================================
    # XELOX Resistance Genes Expression
    # ============================================================
    print("\n[5] XELOX resistance genes distribution...")
    
    # 10-gene fingerprint
    fingerprint = ['CSNK1G2', 'C5ORF42', 'KAZN', 'KLK6', 'MGA', 'MID2', 'SOX11', 'TAS2R40', 'ZNF451', 'GGT1']
    fp_avail = [g for g in fingerprint if g in adata.var_names]
    print(f"  Fingerprint genes found: {len(fp_avail)}/10")
    
    # 253 meta-analysis genes
    genes_253 = pd.read_csv(os.path.join(TAB_DIR, "253_meta_genes.csv"))
    meta_avail = [g for g in genes_253['gene_symbol'] if g in adata.var_names]
    print(f"  Meta-genes found: {len(meta_avail)}/253")
    
    # Capecitabine activation enzymes
    cape_enz = ['CES1', 'CES2', 'TYMS', 'DPYD', 'DCK', 'CDA']
    cape_avail = [g for g in cape_enz if g in adata.var_names]
    print(f"  Capecitabine enzymes found: {len(cape_avail)}/{len(cape_enz)}")
    
    # Expression by cell type for key genes
    if meta_avail:
        # Top 20 most expressed meta-genes
        meta_expr = adata[:, meta_avail].X.toarray()
        mean_expr = meta_expr.mean(axis=0)
        top_idx = np.argsort(mean_expr)[-20:]
        top_meta = [meta_avail[i] for i in top_idx]
        
        # Dotplot
        sc.pl.dotplot(adata, var_names=top_meta, groupby='cell_type',
                      save=f'_{output_prefix}_top_meta_genes.png',
                      figsize=(10, 5))
        
        # Capecitabine enzyme dotplot
        if cape_avail:
            sc.pl.dotplot(adata, var_names=cape_avail, groupby='cell_type',
                         save=f'_{output_prefix}_cape_enz.png',
                         figsize=(8, 4))
        
        # Fingerprint dotplot
        if fp_avail:
            sc.pl.dotplot(adata, var_names=fp_avail, groupby='cell_type',
                         save=f'_{output_prefix}_fingerprint.png',
                         figsize=(8, 4))
    
    # ============================================================
    # CellChat analysis (requires R)
    # ============================================================
    print("\n[6] CellChat communication analysis (requires R)...")
    print("  To run CellChat, export data and use CellChat R package:")
    print(f"    adata.write('{output_prefix}_for_cellchat.h5ad')")
    
    # Save results
    adata.write(os.path.join(FIG_DIR, f'{output_prefix}_annotated.h5ad'))
    
    print(f"\n[{output_prefix}] Analysis complete!")
    return adata

# ============================================================
# SCRIPT
# ============================================================
if __name__ == '__main__':
    print("""
    ========================================
    scRNA-seq ANALYSIS — READY TO RUN
    ========================================
    
    This script is pre-configured for XELOX resistance gene validation.
    
    TO USE:
    1. Download CRC scRNA-seq data from one of:
       - GSE132465 (primary CRC, 23 patients)
         https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465
       - GSE146771 (CRC atlas)
         https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146771
       
    2. Convert to h5ad format if needed, or update the load path
       
    3. Run:
       python scripts/37_scRNA_seq_validation.py
    
    KEY QUESTIONS TO ANSWER:
    1. Are any of the 10-fingerprint genes enriched in specific cell types?
    2. Are capecitabine metabolism enzymes (CES2, TP, DPYD) enriched in
       tumor epithelium vs stroma?
    3. Are meta-analysis genes differentially expressed between
       resistant vs sensitive organoid clones?
    
    For the oxaliplatin-resistant organoid data (Cellular Oncology 2022,
    PMID: 36136268), data access is through the corresponding authors.
    The scRNA-seq data from this paper uses 10x Genomics.
    """)
    
    # Example usage (commented out)
    # run_scrna_pipeline("path/to/GSE132465_CRC.h5ad", output_prefix="GSE132465")
