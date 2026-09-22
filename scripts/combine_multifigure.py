#!/usr/bin/env python3
"""
Combine individual figures into publication-quality multi-panel figures.
Following journal guidelines: 17.5cm width (two-column), 300 DPI, Arial font.

Output multi-panel figures:
  - Figure 1 (merged): Pipeline + SHAP interpretability (2 panels: A, B)
  - Figure 2 (merged): Pathway-level analysis (3 panels: A, B, C)  
  - Figure 3 (merged): Nomogram performance (4 panels: A, B, C, D)
  - Figure 4 (merged): Model comparison (2 panels: A, B)
  - Figure 5 (merged): Meta-analysis enrichment (2 panels: A, B)
"""
import os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
import numpy as np

# Paths
BASE = r'/path/to/xelox_project'
SRC = os.path.join(BASE, 'results', 'figures_png', 'hd')
OUT = os.path.join(BASE, 'results', 'figures_png', 'multifig')
os.makedirs(OUT, exist_ok=True)

# Publication style constants
DPI = 300
FONT_FAMILY = 'Arial'
LABEL_SIZE = 10
TICK_SIZE = 8
TITLE_SIZE = 11

plt.rcParams.update({
    'font.family': FONT_FAMILY,
    'font.size': TICK_SIZE,
    'axes.titlesize': TITLE_SIZE,
    'axes.labelsize': LABEL_SIZE,
    'xtick.labelsize': TICK_SIZE,
    'ytick.labelsize': TICK_SIZE,
    'legend.fontsize': TICK_SIZE,
    'figure.dpi': DPI,
    'savefig.dpi': DPI,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
})


def add_panel_label(ax, label, x=-0.08, y=1.08):
    """Add bold uppercase panel label (A, B, C, ...) to axes."""
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=14, fontweight='bold', va='top', ha='left',
            fontfamily=FONT_FAMILY)


def load_img(filename):
    """Load image from hd directory."""
    path = os.path.join(SRC, filename)
    if not os.path.exists(path):
        print(f'  [WARN] Missing: {filename}')
        return None
    return mpimg.imread(path)


def create_figure3_nomogram():
    """
    Figure 3: Nomogram performance panels
    A: ROC curves (PRS validation)
    B: KM curves
    C: Forest plot (8-variable Cox)
    D: Nomogram
    """
    print('Creating Figure 3: Nomogram Performance...')
    
    imgs = {
        'A': load_img('Fig4A_PRS_ROC.png'),
        'B': load_img('Fig4B_prs_km_curve.png'),
        'C': load_img('Fig4C_cox_forest.png'),
        'D': load_img('Fig4D_nomogram.png'),
    }
    
    if any(v is None for v in imgs.values()):
        print('  [SKIP] Missing source images')
        return
    
    # 2x2 layout, 17.5cm width ≈ 6.9 inches
    fig = plt.figure(figsize=(6.9, 7.5))
    gs = GridSpec(2, 2, figure=fig, hspace=0.25, wspace=0.20)
    
    panels = [
        ('A', imgs['A'], 'LASSO-PRS ROC Curves'),
        ('B', imgs['B'], 'Kaplan-Meier Risk Stratification'),
        ('C', imgs['C'], '8-Variable Cox Forest Plot'),
        ('D', imgs['D'], 'Prognostic Nomogram'),
    ]
    
    for idx, (label, img, title) in enumerate(panels):
        ax = fig.add_subplot(gs[idx // 2, idx % 2])
        ax.imshow(img, aspect='auto')
        ax.set_title(title, fontsize=TITLE_SIZE, fontweight='bold', pad=6)
        ax.axis('off')
        add_panel_label(ax, label)
    
    out_path = os.path.join(OUT, 'Fig3_nomogram_performance.png')
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)
    print(f'  Saved: {out_path}')


def create_figure4_model_comparison():
    """
    Figure 4: Gene-level vs Pathway-level model comparison
    A: Fair comparison forest (gene vs pathway C-index)
    B: Method comparison ablation
    """
    print('Creating Figure 4: Model Comparison...')
    
    imgs = {
        'A': load_img('Fig6A_ablation_forest.png'),
        'B': load_img('Fig6B_method_comparison.png'),
    }
    
    if any(v is None for v in imgs.values()):
        print('  [SKIP] Missing source images')
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(6.9, 3.8))
    
    for ax, (label, img, title) in zip(axes, [
        ('A', imgs['A'], 'Gene-Level vs Pathway-Level C-index'),
        ('B', imgs['B'], 'Fair Comparison: Feature Type × Modeling Strategy'),
    ]):
        ax.imshow(img, aspect='auto')
        ax.set_title(title, fontsize=TITLE_SIZE, fontweight='bold', pad=6)
        ax.axis('off')
        add_panel_label(ax, label)
    
    out_path = os.path.join(OUT, 'Fig4_model_comparison.png')
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)
    print(f'  Saved: {out_path}')


def create_figure2_pathway():
    """
    Figure 2: Pathway-level analysis
    A: Volcano plot
    B: Consensus heatmap
    C: Cross-cohort boxplots
    """
    print('Creating Figure 2: Pathway Analysis...')
    
    imgs = {
        'A': load_img('Fig3A_pathway_volcano.png'),
        'B': load_img('Fig3B_pathway_heatmap.png'),
        'C': load_img('Fig3C_pathway_boxplot.png'),
    }
    
    if any(v is None for v in imgs.values()):
        print('  [SKIP] Missing source images')
        return
    
    fig = plt.figure(figsize=(6.9, 8.0))
    gs = GridSpec(2, 2, figure=fig, hspace=0.25, wspace=0.15,
                  height_ratios=[1, 1.2])
    
    # Panel A: Volcano (top-left)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.imshow(imgs['A'], aspect='auto')
    ax_a.set_title('Pathway Volcano Plot', fontsize=TITLE_SIZE, fontweight='bold', pad=6)
    ax_a.axis('off')
    add_panel_label(ax_a, 'A')
    
    # Panel B: Heatmap (top-right)
    ax_b = fig.add_subplot(gs[0, 1])
    ax_b.imshow(imgs['B'], aspect='auto')
    ax_b.set_title('Cross-Cohort Consensus Heatmap', fontsize=TITLE_SIZE, fontweight='bold', pad=6)
    ax_b.axis('off')
    add_panel_label(ax_b, 'B')
    
    # Panel C: Boxplots (bottom, full width)
    ax_c = fig.add_subplot(gs[1, :])
    ax_c.imshow(imgs['C'], aspect='auto')
    ax_c.set_title('Cross-Cohort Pathway Activity Boxplots', fontsize=TITLE_SIZE, fontweight='bold', pad=6)
    ax_c.axis('off')
    add_panel_label(ax_c, 'C')
    
    out_path = os.path.join(OUT, 'Fig2_pathway_analysis.png')
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)
    print(f'  Saved: {out_path}')


def create_figure5_enrichment():
    """
    Figure 5: Meta-analysis enrichment
    A: BP dotplot
    B: KEGG barplot
    """
    print('Creating Figure 5: Enrichment Analysis...')
    
    imgs = {
        'A': load_img('Fig5_BP_dotplot.png'),
        'B': load_img('Fig5_KEGG_barplot.png'),
    }
    
    if any(v is None for v in imgs.values()):
        print('  [SKIP] Missing source images')
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(6.9, 4.5))
    
    for ax, (label, img, title) in zip(axes, [
        ('A', imgs['A'], 'GO Biological Process Enrichment'),
        ('B', imgs['B'], 'KEGG Pathway Enrichment'),
    ]):
        ax.imshow(img, aspect='auto')
        ax.set_title(title, fontsize=TITLE_SIZE, fontweight='bold', pad=6)
        ax.axis('off')
        add_panel_label(ax, label)
    
    out_path = os.path.join(OUT, 'Fig5_enrichment.png')
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)
    print(f'  Saved: {out_path}')


def create_calibration_panel():
    """
    Figure 3E-F: Calibration curves (3 timepoints)
    """
    print('Creating Figure 3EF: Calibration Curves...')
    
    imgs = {
        'E': load_img('Fig4E_calibration_12mo.png'),
        'F': load_img('Fig4E_calibration_36mo.png'),
        'G': load_img('Fig4F_calibration_60mo.png'),
    }
    
    if any(v is None for v in imgs.values()):
        print('  [SKIP] Missing source images')
        return
    
    fig, axes = plt.subplots(1, 3, figsize=(6.9, 2.8))
    
    for ax, (label, img, title) in zip(axes, [
        ('E', imgs['E'], '12-Month Calibration'),
        ('F', imgs['F'], '36-Month Calibration'),
        ('G', imgs['G'], '60-Month Calibration'),
    ]):
        ax.imshow(img, aspect='auto')
        ax.set_title(title, fontsize=TITLE_SIZE-1, fontweight='bold', pad=4)
        ax.axis('off')
        add_panel_label(ax, label, y=1.05)
    
    out_path = os.path.join(OUT, 'Fig3EFG_calibration.png')
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)
    print(f'  Saved: {out_path}')


if __name__ == '__main__':
    print('=' * 60)
    print('Multi-Panel Figure Generator')
    print(f'Source: {SRC}')
    print(f'Output: {OUT}')
    print('=' * 60)
    print()
    
    create_figure3_nomogram()
    create_figure4_model_comparison()
    create_figure2_pathway()
    create_figure5_enrichment()
    create_calibration_panel()
    
    print()
    print('=' * 60)
    print('All multi-panel figures generated.')
    print('=' * 60)
