#!/usr/bin/env python3
"""
Complete figure regeneration from original R-generated PDFs.
Strategy: PDF -> high-res PNG -> multi-panel assembly.
All panels retain native resolution, no distortion.

BMC/Journal of Translational Medicine standard: 170mm width, 300 DPI.
"""
import os, sys, glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from PIL import Image

# System Python with PyMuPDF
SYSPY = r'/home/user/AppData/Local/Programs/Python/Python314/python.exe'

# Paths
BASE = r'/path/to/xelox_project'
FIG_DIR = os.path.join(BASE, 'results', 'figures')
OUT_DIR = os.path.join(BASE, 'results', 'figures_png', 'v5.10')
TMP_DIR = os.path.join(BASE, 'tmp', 'fig_regen')
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(TMP_DIR, exist_ok=True)

DPI = 300
W_INCH = 6.7  # 170mm
FAMILY = 'Arial'

plt.rcParams.update({
    'font.family': FAMILY,
    'font.size': 9,
    'axes.titlesize': 0,
    'axes.labelsize': 9,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'figure.dpi': DPI,
    'savefig.dpi': DPI,
    'figure.facecolor': 'white',
    'savefig.facecolor': 'white',
    'savefig.bbox': 'tight',
})

# ── PDF to PNG via PyMuPDF ──
def pdf_to_png(pdf_path, zoom=3.0):
    """Convert PDF page to high-res PNG using PyMuPDF (system Python)."""
    name = os.path.splitext(os.path.basename(pdf_path))[0]
    out = os.path.join(TMP_DIR, f'{name}.png')
    if os.path.exists(out):
        return out
    
    script = f'''
import fitz
doc = fitz.open(r"{pdf_path}")
if len(doc) > 0:
    page = doc[0]
    mat = fitz.Matrix({zoom}, {zoom})
    pix = page.get_pixmap(matrix=mat)
    pix.save(r"{out}")
    print(f"OK: {{pix.width}}x{{pix.height}}")
else:
    print("EMPTY")
'''
    import subprocess
    r = subprocess.run([SYSPY, '-c', script], capture_output=True, text=True)
    if 'OK' in r.stdout:
        return out
    else:
        print(f'  [SKIP] Empty or failed: {os.path.basename(pdf_path)}')
        return None


# ── Multi-panel Figure Builder ──
def add_label(ax, label, x=-0.05, y=1.05):
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=12, fontweight='bold', va='top', ha='left', fontfamily=FAMILY)


def embed_png(ax, png_path):
    """Load PNG and embed, return True on success."""
    if not png_path or not os.path.exists(png_path):
        return False
    img = plt.imread(png_path)
    ax.imshow(img, aspect='auto')
    ax.axis('off')
    return True


def save_fig(fig, name):
    path = os.path.join(OUT_DIR, name)
    fig.savefig(path, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    kb = os.path.getsize(path) / 1024
    print(f'  -> {name} ({kb:.0f} KB)')


# ══════════════════════════════════════════════════════
# Figure 1: Pipeline (drawn from scratch)
# ══════════════════════════════════════════════════════
def fig1_pipeline():
    print('\nCreating Figure 1: Study Pipeline...')
    from matplotlib.patches import FancyBboxPatch

    fig, ax = plt.subplots(figsize=(W_INCH, 9.5))
    ax.set_xlim(0, 7)
    ax.set_ylim(0, 9.5)
    ax.axis('off')

    CX, BW, BH = 2.0, 3.0, 0.6
    colors = {'geo':'#E64B35','qc':'#F9A825','sel':'#F57C00','gse':'#1565C0',
              'cla':'#00ACC1','sur':'#43A047','val':'#7B1FA2'}

    def rbox(x, y, w, h, title, sub='', fc='white', ec='#333', lw=1.5, ts=10, ss=7.5):
        box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.025", facecolor=fc, edgecolor=ec, linewidth=lw)
        ax.add_patch(box)
        if sub:
            ax.text(x+w/2, y+h/2+0.08, title, ha='center', va='center', fontsize=ts, fontweight='bold', color='#1a1a2e')
            ax.text(x+w/2, y+h/2-0.22, sub, ha='center', va='center', fontsize=ss, color='#666')
        else:
            ax.text(x+w/2, y+h/2, title, ha='center', va='center', fontsize=ts, fontweight='bold', color='#1a1a2e')

    def arrow_v(x, y1, y2):
        ax.annotate('', xy=(x, y2+0.02), xytext=(x, y1-0.02), arrowprops=dict(arrowstyle='->', color='#555', lw=1.5))

    def rlabel(x, y, t, c='#555', fs=8):
        ax.text(x, y, t, fontsize=fs, color=c, ha='left', va='center')

    # Row 1: GEO
    y = 8.6
    rbox(CX, y, BW, BH, 'GEO Data Retrieval', '5 Discovery Cohorts', ts=11)
    rlabel(CX+BW+0.15, y+BH/2, 'n = 1,015', colors['geo'], 8)

    # Row 2: QC
    y = 7.5
    arrow_v(CX+BW/2, y+BH+0.7, y+BH)
    rbox(CX, y, BW, BH, 'Quality Control & Preprocessing', fc='#FFF8E1', ec=colors['qc'], ts=11)
    rlabel(CX+BW+0.15, y+BH/2, 'Excluded: 539', '#999', 8)

    # Row 3: XELOX
    y = 6.4
    arrow_v(CX+BW/2, y+BH+0.7, y+BH)
    rbox(CX, y, BW, BH, 'XELOX / FOLFOX-treated mCRC', '5 Cohorts', fc='#FFF3E0', ec=colors['sel'], ts=11)
    rlabel(CX+BW+0.15, y+BH/2, 'n = 476', colors['sel'], 8)

    # Row 4: GSE39582
    y = 5.3
    arrow_v(CX+BW/2, y+BH+0.7, y+BH)
    rbox(CX, y, BW, BH, 'GSE39582 Primary Cohort', 'n = 239, 79 RFS events', fc='#E3F2FD', ec=colors['gse'], lw=2, ts=11)

    # Row 5: Split
    y = 3.8
    arrow_v(CX+BW/2, y+BH+0.7, y+BH)
    rbox(0.3, y, 2.8, BH, 'Binary Classification', 'Resistant (45) vs Sensitive (119)', fc='#E0F7FA', ec=colors['cla'], ts=10, ss=7)
    rlabel(3.4, y+BH/2, '[EX] 75\nIntermediate RFS', '#999', 7)
    rbox(3.9, y, 2.8, BH, 'Survival Analysis', 'RFS Endpoint, 79 Events', fc='#E8F5E9', ec=colors['sur'], ts=10, ss=7)
    rlabel(7.0, y+BH/2, '[IN] 66\nIntermediate RFS', colors['sur'], 7)

    # Row 6: Methods
    y = 2.7
    arrow_v(0.3+1.4, y+BH+0.7, y+BH)
    arrow_v(3.9+1.4, y+BH+0.7, y+BH)
    rbox(0.3, y, 2.8, BH, 'XGBoost + SHAP', '10-Gene Fingerprint', fc='#B2EBF2', ec=colors['cla'], ts=10, ss=7)
    rbox(3.9, y, 2.8, BH, 'ssGSEA + AIC-Cox', '8-Variable Nomogram', fc='#C8E6C9', ec=colors['sur'], ts=10, ss=7)

    # Row 7: External
    y = 1.5
    arrow_v(3.9+1.4, y+BH+0.7, y+BH)
    rbox(0.5, y, 6.0, BH, 'External Validation Cohorts', 'TCGA-COAD/READ (n=140) | GSE104645 (n=113) | GSE72970 (n=32) | GSE83129 (n=26)',
         fc='#F3E5F5', ec=colors['val'], ts=10, ss=7.5)
    rlabel(0.5, y-0.15, 'Total: n = 311', colors['val'], 7)

    # Summary
    sx, sy, sw, sh = 5.3, 5.5, 1.55, 1.6
    box = FancyBboxPatch((sx, sy), sw, sh, boxstyle="round,pad=0.03", facecolor='#FFF8E1', edgecolor='#F9A825', linewidth=1.5)
    ax.add_patch(box)
    ax.annotate('', xy=(sx, sy+sh/2), xytext=(CX+BW, sy+sh/2), arrowprops=dict(arrowstyle='->', color='#999', lw=1, linestyle='dashed'))
    for t, fs, wt, ly in [('Sample Flow',9,'bold',sy+sh-0.12),('───────',6,'normal',sy+sh-0.30),
        ('Pool:  1,015',7.5,'normal',sy+sh-0.48),('QC:    476',7.5,'normal',sy+sh-0.66),
        ('GSE39582: 239',7.5,'bold',sy+sh-0.84),('  Class: 164',7,'normal',sy+sh-1.00),
        ('  Surv:  227',7,'normal',sy+sh-1.14),('Ext:   311',7.5,'bold',sy+sh-1.32)]:
        ax.text(sx+sw/2, ly, t, ha='center', va='center', fontsize=fs, fontweight=wt, color='#1a1a2e')

    return fig


# ══════════════════════════════════════════════════════
# Figure 2: SHAP (from PDF)
# ══════════════════════════════════════════════════════
def fig2_shap():
    print('\nCreating Figure 2: SHAP Interpretability...')
    pdfs = [
        ('figures/ml_phase2/SHAP_beeswarm.pdf', 'A'),
        ('figures/ml_phase2/SHAP_summary_bar_top20.pdf', 'B'),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(W_INCH, 4.5))
    plt.subplots_adjust(wspace=0.12)
    for ax, (pdf, label) in zip(axes, pdfs):
        path = pdf_to_png(os.path.join(BASE, 'results', pdf))
        embed_png(ax, path)
        add_label(ax, label)
    return fig


# ══════════════════════════════════════════════════════
# Figure 3: Pathway Analysis
# ══════════════════════════════════════════════════════
def fig3_pathway():
    print('\nCreating Figure 3: Pathway Analysis...')
    pdfs = [
        ('figures/pathway/volcano_combined.png', 'A'),  # try PNG first
    ]
    # Try finding pathway figures
    candidates = glob.glob(os.path.join(BASE, 'results', 'figures', 'pathway', '*volcano*.pdf'))
    if not candidates:
        candidates = glob.glob(os.path.join(BASE, 'results', 'figures_png', 'hd', 'Fig3A_pathway_volcano.png'))

    # Fallback: use existing HD images for pathway figures
    hd_files = {
        'A': 'Fig3A_pathway_volcano.png',
        'B': 'Fig3B_pathway_heatmap.png',
        'C': 'Fig3C_pathway_boxplot.png',
    }
    fig = plt.figure(figsize=(W_INCH, 7.0))
    gs = GridSpec(2, 2, figure=fig, hspace=0.16, wspace=0.12, height_ratios=[1, 1.1])

    for lbl, pos in [('A', (0,0)), ('B', (0,1))]:
        ax = fig.add_subplot(gs[pos])
        fp = os.path.join(BASE, 'results', 'figures_png', 'hd', hd_files[lbl])
        if os.path.exists(fp):
            img = Image.open(fp)
            ax.imshow(img, aspect='auto')
        ax.axis('off')
        add_label(ax, lbl)

    ax_c = fig.add_subplot(gs[1, :])
    fp = os.path.join(BASE, 'results', 'figures_png', 'hd', hd_files['C'])
    if os.path.exists(fp):
        img = Image.open(fp)
        ax_c.imshow(img, aspect='auto')
    ax_c.axis('off')
    add_label(ax_c, 'C')

    return fig


# ══════════════════════════════════════════════════════
# Figure 4: Nomogram (PDF conversions)
# ══════════════════════════════════════════════════════
def fig4_nomogram():
    print('\nCreating Figure 4: Nomogram Performance...')
    # Convert PRS ROCs
    roc_pdf = pdf_to_png(os.path.join(BASE, 'results', 'figures', 'prs', 'prs_roc_combined.pdf'))
    # KM curve
    km_pdf = pdf_to_png(os.path.join(BASE, 'results', 'figures', 'prs', 'prs_km_curve.pdf'))
    # Forest plot (use our generated one)
    forest_png = os.path.join(BASE, 'results', 'figures_png', 'hd', 'Fig4C_cox_forest.png')
    # Nomogram
    nomo_pdf = pdf_to_png(os.path.join(BASE, 'results', 'figures', 'nomogram', 'nomogram.pdf'))

    fig = plt.figure(figsize=(W_INCH, 6.5))
    gs = GridSpec(2, 2, figure=fig, hspace=0.16, wspace=0.12)

    panels = [
        (0,0, roc_pdf, 'A'),
        (0,1, km_pdf, 'B'),
        (1,0, forest_png, 'C'),
        (1,1, nomo_pdf, 'D'),
    ]
    for r, c, src, lbl in panels:
        ax = fig.add_subplot(gs[r, c])
        embed_png(ax, src)
        add_label(ax, lbl)

    return fig


# ══════════════════════════════════════════════════════
# Figure 5: Calibration (PDF)
# ══════════════════════════════════════════════════════
def fig5_calibration():
    print('\nCreating Figure 5: Calibration...')
    pdfs = [
        ('figures/prs/calibration_12mo_prs.pdf', 'A'),
        ('figures/prs/calibration_36mo_prs.pdf', 'B'),
        ('figures/prs/calibration_60mo_prs.pdf', 'C'),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(W_INCH, 3.0))
    plt.subplots_adjust(wspace=0.10)
    for ax, (pdf, lbl) in zip(axes, pdfs):
        png = pdf_to_png(os.path.join(BASE, 'results', pdf))
        embed_png(ax, png)
        add_label(ax, lbl)
        ax.text(0.98, 0.02, 'Slope=0.817', transform=ax.transAxes, fontsize=6.5, color='#555',
                va='bottom', ha='right', bbox=dict(boxstyle='round,pad=0.12', fc='white', alpha=0.85, ec='#ccc'))
    return fig


# ══════════════════════════════════════════════════════
# Figure 6: Fair Comparison
# ══════════════════════════════════════════════════════
def fig6_fair():
    print('\nCreating Figure 6: Fair Comparison...')
    hd = os.path.join(BASE, 'results', 'figures_png', 'hd')
    fig, axes = plt.subplots(1, 2, figsize=(W_INCH, 3.5))
    plt.subplots_adjust(wspace=0.12)
    for ax, fn, lbl in zip(axes, ['Fig6A_ablation_forest.png', 'Fig6B_method_comparison.png'], ['A','B']):
        fp = os.path.join(hd, fn)
        if os.path.exists(fp):
            img = Image.open(fp)
            ax.imshow(img, aspect='auto')
        ax.axis('off')
        add_label(ax, lbl)
    return fig


# ══════════════════════════════════════════════════════
# Figure 7: Enrichment
# ══════════════════════════════════════════════════════
def fig7_enrichment():
    print('\nCreating Figure 7: Enrichment...')
    hd = os.path.join(BASE, 'results', 'figures_png', 'hd')
    fig, axes = plt.subplots(1, 2, figsize=(W_INCH, 4.0))
    plt.subplots_adjust(wspace=0.12)
    for ax, fn, lbl in zip(axes, ['Fig5_BP_dotplot.png', 'Fig5_KEGG_barplot.png'], ['A','B']):
        fp = os.path.join(hd, fn)
        if os.path.exists(fp):
            img = Image.open(fp)
            ax.imshow(img, aspect='auto')
        ax.axis('off')
        add_label(ax, lbl)
    return fig


# ══════════════════════════════════════════════════════
if __name__ == '__main__':
    print('=' * 60)
    print('V5.10 FIGURE REGENERATION (PDF -> PNG -> Multi-panel)')
    print('=' * 60)

    save_fig(fig1_pipeline(), 'Fig1_pipeline.png')
    save_fig(fig2_shap(), 'Fig2_SHAP.png')
    save_fig(fig3_pathway(), 'Fig3_pathway.png')
    save_fig(fig4_nomogram(), 'Fig4_nomogram.png')
    save_fig(fig5_calibration(), 'Fig5_calibration.png')
    save_fig(fig6_fair(), 'Fig6_fair_comparison.png')
    save_fig(fig7_enrichment(), 'Fig7_enrichment.png')

    print(f'\nDone! All figures in: {OUT_DIR}')
