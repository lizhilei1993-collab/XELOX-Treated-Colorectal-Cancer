#!/usr/bin/env python3
"""
Regenerate ALL main figures for V5.8 manuscript.
Following BMC/Springer Nature Journal of Translational Medicine guidelines:
  - Size: 170mm width (single, landscape) or 85mm (single, portrait)
  - Resolution: 300 DPI minimum
  - Font: Arial, 8-10pt, bold panel labels (A, B, C)
  - Format: High-quality PNG 300 DPI
  - No titles on figures themselves (pure data presentation)
  - Panel labels: bold uppercase, top-left corner
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch
from PIL import Image

# ── Global Settings ────────────────────────────────────────────
DPI = 300
FAMILY = 'Arial'
LABEL_SIZE = 9
TICK_SIZE = 7
# Page width = 170mm ≈ 6.7 inches (BMC double-column)
FULL_WIDTH_INCHES = 6.7
SINGLE_WIDTH_MM = 85
FULL_WIDTH_MM = 170

plt.rcParams.update({
    'font.family': FAMILY,
    'font.size': LABEL_SIZE,
    'axes.labelsize': LABEL_SIZE,
    'axes.titlesize': 0,          # no title
    'xtick.labelsize': TICK_SIZE,
    'ytick.labelsize': TICK_SIZE,
    'legend.fontsize': 7,
    'figure.dpi': DPI,
    'savefig.dpi': DPI,
    'savefig.bbox': 'tight',
    'savefig.facecolor': 'white',
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'axes.grid': False,
})

# Paths
BASE = r'/path/to/xelox_project'
SRC_HD = os.path.join(BASE, 'results', 'figures_png', 'hd')
SRC_RAW = os.path.join(BASE, 'results', 'figures_png')
OUT = os.path.join(BASE, 'results', 'figures_png', 'v5.9')
os.makedirs(OUT, exist_ok=True)

TABLES_DIR = os.path.join(BASE, 'results', 'tables')

def add_label(ax, label, x=-0.06, y=1.06):
    """Add bold uppercase panel label."""
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=11, fontweight='bold', va='top', ha='left',
            fontfamily=FAMILY)

def add_stats(ax, text, x=0.98, y=0.02):
    """Add small statistics annotation bottom-right."""
    ax.text(x, y, text, transform=ax.transAxes,
            fontsize=6.5, va='bottom', ha='right', color='#555',
            fontfamily=FAMILY,
            bbox=dict(boxstyle='round,pad=0.15', facecolor='white', alpha=0.85, edgecolor='#ccc'))

def load_and_embed(ax, filename, src_dir=SRC_HD):
    """Load PNG image and embed into axes, removing any existing titles."""
    fp = os.path.join(src_dir, filename)
    if not os.path.exists(fp):
        fp = os.path.join(SRC_RAW, filename)
    if not os.path.exists(fp):
        print(f'  [WARN] Missing: {filename}')
        return False
    img = mpimg.imread(fp)
    ax.imshow(img, aspect='auto')
    ax.axis('off')
    return True


# ================================================================
# Figure 1: Study Pipeline (single column, 170mm full width)
# ================================================================
def create_figure1():
    print('Creating Figure 1: Study Pipeline...')
    fig, ax = plt.subplots(figsize=(FULL_WIDTH_INCHES, 9.5))
    ax.set_xlim(0, 7)
    ax.set_ylim(0, 9.5)
    ax.axis('off')

    CX, BW, BH = 2.0, 3.0, 0.6

    def rbox(x, y, w, h, title, sub='', fc='white', ec='#333', lw=1.5, ts=10, ss=7.5):
        box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.03",
                             facecolor=fc, edgecolor=ec, linewidth=lw)
        ax.add_patch(box)
        if sub:
            ax.text(x+w/2, y+h/2+0.06, title, ha='center', va='center',
                    fontsize=ts, fontweight='bold', color='#1a1a2e')
            ax.text(x+w/2, y+h/2-0.22, sub, ha='center', va='center',
                    fontsize=ss, color='#666')
        else:
            ax.text(x+w/2, y+h/2, title, ha='center', va='center',
                    fontsize=ts, fontweight='bold', color='#1a1a2e')

    def arrow_v(x, y1, y2):
        ax.annotate('', xy=(x, y2+0.02), xytext=(x, y1-0.02),
                    arrowprops=dict(arrowstyle='->', color='#555', lw=1.5))

    def right_label(x, y, text, color='#555', fs=8):
        ax.text(x, y, text, fontsize=fs, color=color, ha='left', va='center')

    colors = {
        'geo': '#E64B35', 'qc': '#F9A825', 'sel': '#F57C00',
        'gse': '#1565C0', 'cla': '#00ACC1', 'sur': '#43A047', 'val': '#7B1FA2'
    }

    # Row 1
    y = 8.6
    rbox(CX, y, BW, BH, 'GEO Data Retrieval', '5 Discovery Cohorts', ts=11)
    right_label(CX+BW+0.15, y+BH/2, 'n = 1,015', colors['geo'], 8)

    # Row 2
    y = 7.5
    arrow_v(CX+BW/2, y+BH+0.7, y+BH)
    rbox(CX, y, BW, BH, 'Quality Control & Preprocessing', '', fc='#FFF8E1', ec=colors['qc'], ts=11)
    right_label(CX+BW+0.15, y+BH/2, 'Excluded: 539', '#999', 8)

    # Row 3
    y = 6.4
    arrow_v(CX+BW/2, y+BH+0.7, y+BH)
    rbox(CX, y, BW, BH, 'XELOX / FOLFOX-treated mCRC', '5 Cohorts Combined', fc='#FFF3E0', ec=colors['sel'], ts=11)
    right_label(CX+BW+0.15, y+BH/2, 'n = 476', colors['sel'], 8)

    # Row 4
    y = 5.3
    arrow_v(CX+BW/2, y+BH+0.7, y+BH)
    rbox(CX, y, BW, BH, 'GSE39582 Discovery Cohort', 'Primary Analysis', fc='#E3F2FD', ec=colors['gse'], lw=2, ts=11)
    right_label(CX+BW+0.15, y+BH/2, 'n = 239 (79 events)', colors['gse'], 8)

    # Row 5: Split
    y = 3.8
    arrow_v(CX+BW/2, y+BH+0.7, y+BH)

    # Left
    rbox(0.3, y, 2.8, BH, 'Binary Classification', 'Resistant (45) vs Sensitive (119)', fc='#E0F7FA', ec=colors['cla'], ts=10, ss=7)
    right_label(3.4, y+BH/2, '[EX] Excluded: 75', '#999', 7)

    # Right
    rbox(3.9, y, 2.8, BH, 'Survival Analysis', 'RFS Endpoint, 79 Events', fc='#E8F5E9', ec=colors['sur'], ts=10, ss=7)
    right_label(7.0, y+BH/2, '[IN] Included: 66', colors['sur'], 7)

    # Row 6
    y = 2.7
    arrow_v(0.3+1.4, y+BH+0.7, y+BH)
    arrow_v(3.9+1.4, y+BH+0.7, y+BH)

    rbox(0.3, y, 2.8, BH, 'XGBoost + SHAP', '10-Gene Fingerprint (Table 1)', fc='#B2EBF2', ec=colors['cla'], ts=10, ss=7)
    rbox(3.9, y, 2.8, BH, 'ssGSEA + AIC-Cox', '8-Variable Nomogram (Table 2)', fc='#C8E6C9', ec=colors['sur'], ts=10, ss=7)

    # Row 7
    y = 1.5
    arrow_v(3.9+1.4, y+BH+0.7, y+BH)
    rbox(0.5, y, 6.0, BH, 'External Validation Cohorts', 'TCGA-COAD/READ (n=140) | GSE104645 (n=113) | GSE72970 (n=32) | GSE83129 (n=26)',
         fc='#F3E5F5', ec=colors['val'], ts=11, ss=7.5)
    right_label(0.5, y-0.15, 'Total: n = 311', colors['val'], 7)

    # Summary box
    sx, sy = 5.3, 5.5
    sw, sh = 1.55, 1.6
    box = FancyBboxPatch((sx, sy), sw, sh, boxstyle="round,pad=0.03",
                         facecolor='#FFF8E1', edgecolor='#F9A825', linewidth=1.5)
    ax.add_patch(box)
    ax.annotate('', xy=(sx, sy+sh/2), xytext=(CX+BW, sy+sh/2),
                arrowprops=dict(arrowstyle='->', color='#999', lw=1, linestyle='dashed'))
    lines = [
        ('Sample Flow', 9, 'bold', sy+sh-0.12), ('────────', 6, 'normal', sy+sh-0.30),
        ('Pool:    1,015', 7.5, 'normal', sy+sh-0.48), ('QC:       476', 7.5, 'normal', sy+sh-0.66),
        ('GSE39582:  239', 7.5, 'bold', sy+sh-0.84), ('  Class: 164', 7, 'normal', sy+sh-1.00),
        ('  Surv:  227', 7, 'normal', sy+sh-1.14), ('External:  311', 7.5, 'bold', sy+sh-1.32),
    ]
    for text, fs, wt, ly in lines:
        ax.text(sx+sw/2, ly, text, ha='center', va='center', fontsize=fs, fontweight=wt, color='#1a1a2e')

    add_label(ax, 'A')
    return fig


# ================================================================
# Figure 2: SHAP + Pathway (merged panels)
# ================================================================
def create_figure2():
    print('Creating Figure 2: SHAP Interpretability...')
    fig = plt.figure(figsize=(FULL_WIDTH_INCHES, 4.5))
    gs = GridSpec(1, 2, figure=fig, wspace=0.15)

    ax_a = fig.add_subplot(gs[0, 0])
    load_and_embed(ax_a, 'Fig2_SHAP_beeswarm.png')
    add_label(ax_a, 'A')

    ax_b = fig.add_subplot(gs[0, 1])
    load_and_embed(ax_b, 'Fig2_SHAP_bar.png')
    add_label(ax_b, 'B')

    return fig


# ================================================================
# Figure 3: Pathway Analysis (3 panels)
# ================================================================
def create_figure3():
    print('Creating Figure 3: Pathway Analysis...')
    fig = plt.figure(figsize=(FULL_WIDTH_INCHES, 7.0))
    gs = GridSpec(2, 2, figure=fig, hspace=0.2, wspace=0.15, height_ratios=[1, 1.1])

    ax_a = fig.add_subplot(gs[0, 0])
    load_and_embed(ax_a, 'Fig3A_pathway_volcano.png')
    add_label(ax_a, 'A')

    ax_b = fig.add_subplot(gs[0, 1])
    load_and_embed(ax_b, 'Fig3B_pathway_heatmap.png')
    add_label(ax_b, 'B')

    ax_c = fig.add_subplot(gs[1, :])
    load_and_embed(ax_c, 'Fig3C_pathway_boxplot.png')
    add_label(ax_c, 'C')

    return fig


# ================================================================
# Figure 4: Nomogram (4 panels: ROC + KM + Forest + Nomogram)
# ================================================================
def create_figure4():
    print('Creating Figure 4: Nomogram Performance...')
    fig = plt.figure(figsize=(FULL_WIDTH_INCHES, 6.5))
    gs = GridSpec(2, 2, figure=fig, hspace=0.18, wspace=0.12)

    panels = [
        (0, 0, 'Fig4A_PRS_ROC.png', 'A'),
        (0, 1, 'Fig4B_prs_km_curve.png', 'B'),
        (1, 0, 'Fig4C_cox_forest.png', 'C'),
        (1, 1, 'Fig4D_nomogram.png', 'D'),
    ]
    for r, c, fn, lbl in panels:
        ax = fig.add_subplot(gs[r, c])
        load_and_embed(ax, fn)
        add_label(ax, lbl)

    return fig


# ================================================================
# Figure 5: Calibration (3 panels)
# ================================================================
def create_figure5():
    print('Creating Figure 5: Calibration Curves...')
    fig, axes = plt.subplots(1, 3, figsize=(FULL_WIDTH_INCHES, 3.0))
    plt.subplots_adjust(wspace=0.12)

    panels = [
        (0, 'Fig4E_calibration_12mo.png', 'A'),
        (1, 'Fig4E_calibration_36mo.png', 'B'),
        (2, 'Fig4F_calibration_60mo.png', 'C'),
    ]
    metrics = {
        0: 'Slope=0.817\nBrier=0.142',
        1: 'Slope=0.817\nBrier=0.201',
        2: 'Slope=0.817\nBrier=0.218',
    }
    for idx, fn, lbl in panels:
        ax = axes[idx]
        load_and_embed(ax, fn)
        add_label(ax, lbl)
        add_stats(ax, metrics[idx])

    return fig


# ================================================================
# Figure 6: Fair Comparison (2 panels)
# ================================================================
def create_figure6():
    print('Creating Figure 6: Fair Comparison...')
    fig, axes = plt.subplots(1, 2, figsize=(FULL_WIDTH_INCHES, 3.5))
    plt.subplots_adjust(wspace=0.12)

    for ax, fn, lbl in zip(axes,
                           ['Fig6A_ablation_forest.png', 'Fig6B_method_comparison.png'],
                           ['A', 'B']):
        load_and_embed(ax, fn)
        add_label(ax, lbl)

    return fig


# ================================================================
# Figure 7: Enrichment (2 panels)
# ================================================================
def create_figure7():
    print('Creating Figure 7: Enrichment Analysis...')
    fig, axes = plt.subplots(1, 2, figsize=(FULL_WIDTH_INCHES, 4.0))
    plt.subplots_adjust(wspace=0.12)

    for ax, fn, lbl in zip(axes,
                           ['Fig5_BP_dotplot.png', 'Fig5_KEGG_barplot.png'],
                           ['A', 'B']):
        load_and_embed(ax, fn)
        add_label(ax, lbl)

    return fig


# ================================================================
# Main
# ================================================================
def save_figure(fig, name):
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    kb = os.path.getsize(path) / 1024
    print(f'  -> {name} ({kb:.0f} KB)')
    return path


if __name__ == '__main__':
    print('=' * 60)
    print('V5.9 Full Figure Regeneration')
    print('Standard: BMC / Journal of Translational Medicine')
    print(f'Width: {FULL_WIDTH_MM}mm, DPI: {DPI}, Font: {FAMILY} {LABEL_SIZE}pt')
    print('=' * 60)
    print()

    save_figure(create_figure1(), 'Fig1_pipeline.png')
    save_figure(create_figure2(), 'Fig2_SHAP_interpretability.png')
    save_figure(create_figure3(), 'Fig3_pathway_analysis.png')
    save_figure(create_figure4(), 'Fig4_nomogram_performance.png')
    save_figure(create_figure5(), 'Fig5_calibration.png')
    save_figure(create_figure6(), 'Fig6_fair_comparison.png')
    save_figure(create_figure7(), 'Fig7_enrichment.png')

    print()
    print('=' * 60)
    print(f'All figures saved to: {OUT}')
    print('Done!')
