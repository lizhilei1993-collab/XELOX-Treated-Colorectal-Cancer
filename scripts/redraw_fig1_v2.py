#!/usr/bin/env python3
"""
Redraw Figure 1: Study Design & Analytical Pipeline (v2)
Professional flowchart following TRIPOD/REMARK guidelines.
Clear sample flow with exclusion counts at each stage.
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

# Publication settings
DPI = 300
FONT = 'Arial'

plt.rcParams.update({
    'font.family': FONT,
    'font.size': 9,
    'figure.dpi': DPI,
    'savefig.dpi': DPI,
    'savefig.bbox': 'tight',
})

# Colors (colorblind-friendly)
C = {
    'discovery': '#E64B35',    # Red
    'qc': '#F9A825',           # Amber
    'primary': '#1565C0',      # Blue
    'path_class': '#4DBBD5',   # Teal
    'path_surv': '#00A087',    # Green
    'external': '#7B1FA2',     # Purple
    'exclude': '#9E9E9E',      # Gray
    'box_bg': '#FAFAFA',
    'border': '#333333',
    'text': '#212121',
    'light_bg': '#FFF8E1',    # Light yellow
}


def rounded_box(ax, x, y, w, h, text, fc='white', ec='#333', lw=1.2, fontsize=8, bold=False, alpha=1.0):
    box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.015", facecolor=fc, edgecolor=ec, linewidth=lw, alpha=alpha)
    ax.add_patch(box)
    weight = 'bold' if bold else 'normal'
    ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=fontsize, fontweight=weight, color=C['text'])


def arrow(ax, x1, y1, x2, y2, color='#333', lw=1.5, style='->', connectionstyle='arc3,rad=0'):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle=style, color=color, lw=lw, connectionstyle=connectionstyle))


def label_box(ax, x, y, text, color='#333', fontsize=7, bg='white'):
    ax.text(x, y, text, ha='center', va='center', fontsize=fontsize, color=color,
            bbox=dict(boxstyle='round,pad=0.2', facecolor=bg, edgecolor='none', alpha=0.9))


def create_figure():
    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.set_xlim(0, 8.5)
    ax.set_ylim(0, 11)
    ax.axis('off')

    # Title
    ax.text(4.25, 10.6, 'Figure 1. Study Design and Analytical Pipeline',
            ha='center', va='center', fontsize=13, fontweight='bold', color=C['primary'])
    ax.text(4.25, 10.35, 'Following TRIPOD Statement for Transparent Reporting',
            ha='center', va='center', fontsize=9, color='gray', fontstyle='italic')

    # ── Row 1: GEO Data Retrieval ──
    y1 = 9.7
    rounded_box(ax, 2.5, y1, 3.5, 0.5, 'GEO Data Retrieval\n5 Discovery Cohorts',
                fc=C['discovery'], ec=C['discovery'], fontsize=10, bold=True)
    label_box(ax, 2.0, y1 + 0.25, 'n = 1,015\nsamples', color=C['discovery'], fontsize=8)

    # ── Row 2: QC ──
    y2 = 8.7
    rounded_box(ax, 2.5, y2, 3.5, 0.5, 'Quality Control &\nPreprocessing',
                fc=C['qc'], ec=C['qc'], fontsize=9, bold=True)
    arrow(ax, 4.25, y1, 4.25, y2 + 0.5, color=C['border'], lw=2)
    label_box(ax, 6.5, 9.2, 'Excluded: 539\n• Non-CRC\n• No chemo info\n• Low quality',
              color=C['exclude'], fontsize=7)

    # ── Row 3: XELOX Subgroup ──
    y3 = 7.7
    rounded_box(ax, 2.5, y3, 3.5, 0.5, 'XELOX/FOLFOX-treated\nmCRC Subgroup',
                fc='#FFE0B2', ec=C['qc'], fontsize=9, bold=True)
    arrow(ax, 4.25, y2, 4.25, y3 + 0.5, color=C['border'], lw=2)
    label_box(ax, 2.0, y3 + 0.25, 'n = 476\n(47%)', color=C['qc'], fontsize=8)

    # ── Row 4: Cohort Distribution ──
    y4 = 6.7
    rounded_box(ax, 0.3, y4, 1.6, 0.5, 'GSE39582\n(n = 239)',
                fc='#FFECB3', ec=C['primary'], fontsize=8, bold=True, lw=1.5)
    rounded_box(ax, 2.2, y4, 1.3, 0.5, 'GSE104645\n(n = 193)',
                fc=C['box_bg'], ec=C['border'], fontsize=7)
    rounded_box(ax, 3.8, y4, 1.1, 0.5, 'GSE28702\n(n = 83)',
                fc=C['box_bg'], ec=C['border'], fontsize=7)
    rounded_box(ax, 5.1, y4, 1.1, 0.5, 'GSE72970\n(n = 124)',
                fc=C['box_bg'], ec=C['border'], fontsize=7)
    rounded_box(ax, 6.4, y4, 1.1, 0.5, 'GSE69657\n(n = 30)',
                fc=C['box_bg'], ec=C['border'], fontsize=7)
    arrow(ax, 4.25, y3, 4.25, y4 + 0.5, color=C['border'], lw=2)

    # ── Row 5: GSE39582 Split ──
    y5 = 5.5
    rounded_box(ax, 1.5, y5, 3.0, 0.7,
                'GSE39582 XELOX Cohort\n(n = 239)',
                fc=C['primary'], ec=C['primary'], fontsize=10, bold=True, lw=2)
    arrow(ax, 1.1, y4 + 0.25, 1.5, y5 + 0.35, color=C['primary'], lw=2)

    # ── Row 6: Two Analysis Paths ──
    y6_class = 4.3
    y6_surv = 2.5

    # Classification path (left)
    rounded_box(ax, 0.2, y6_class, 3.8, 0.8,
                'Binary Classification Task\n(Resistant vs. Sensitive)',
                fc=C['path_class'], ec=C['path_class'], fontsize=9, bold=True)
    arrow(ax, 2.0, y5, 2.1, y6_class + 0.8, color=C['path_class'], lw=2)

    # Classification details
    y7_class = 3.3
    rounded_box(ax, 0.2, y7_class, 3.8, 0.8,
                'GSE39582 (n = 164)\nR = 45  |  S = 119\nIntermediate (n = 75): Excluded',
                fc='#E0F7FA', ec=C['path_class'], fontsize=8, lw=1)
    arrow(ax, 2.1, y6_class, 2.1, y7_class + 0.8, color=C['path_class'], lw=1.5)
    label_box(ax, 4.5, 4.5, 'Excluded: 75 intermediate\n(RFS 12-36 months)', color=C['exclude'], fontsize=7)

    # Survival path (right)
    rounded_box(ax, 4.5, y6_class, 3.8, 0.8,
                'Survival Analysis Task\n(RFS Endpoint)',
                fc=C['path_surv'], ec=C['path_surv'], fontsize=9, bold=True)
    arrow(ax, 6.5, y5, 6.4, y6_class + 0.8, color=C['path_surv'], lw=2)

    # Survival details
    y7_surv = 3.3
    rounded_box(ax, 4.5, y7_surv, 3.8, 0.8,
                'GSE39582 (n = 227)\nEvents: 79 (34.8%)\nIntermediate: Included (n = 66)',
                fc='#E8F5E9', ec=C['path_surv'], fontsize=8, lw=1)
    arrow(ax, 6.4, y6_class, 6.4, y7_surv + 0.8, color=C['path_surv'], lw=1.5)

    # ── Row 7: Methods ──
    y8 = 2.2
    rounded_box(ax, 0.2, y8, 3.8, 0.6,
                'XGBoost → SHAP → 10-gene\nFingerprint (Table 1)',
                fc='#B2EBF2', ec=C['path_class'], fontsize=8, lw=1)
    arrow(ax, 2.1, y7_class, 2.1, y8 + 0.6, color=C['path_class'], lw=1.5)

    rounded_box(ax, 4.5, y8, 3.8, 0.6,
                'ssGSEA → AIC-Cox →\n8-variable Nomogram (Table 2)',
                fc='#C8E6C9', ec=C['path_surv'], fontsize=8, lw=1)
    arrow(ax, 6.4, y7_surv, 6.4, y8 + 0.6, color=C['path_surv'], lw=1.5)

    # ── Row 8: External Validation ──
    y9 = 1.0
    rounded_box(ax, 0.5, y9, 7.5, 0.7,
                'External Validation Cohorts',
                fc=C['external'], ec=C['external'], fontsize=10, bold=True)
    arrow(ax, 4.25, y8, 4.25, y9 + 0.7, color=C['external'], lw=2)

    # Validation details
    y10 = 0.2
    val_text = (
        'TCGA-COAD/READ (Ox, n=140)   •   '
        'GSE104645 (Ox, n=113)   •   '
        'GSE72970 (FOLFOX, n=32)   •   '
        'GSE83129 (FOLFOX, n=26)'
    )
    rounded_box(ax, 0.5, y10, 7.5, 0.5, val_text,
                fc='#F3E5F5', ec=C['external'], fontsize=8)

    # ── Sample Flow Summary Box ──
    summary = (
        'Sample Flow Summary\n'
        '───────────────────────\n'
        'GEO pool: 1,015\n'
        'After QC: 476 (47%)\n'
        'GSE39582: 239\n'
        '  ├ Classification: 164 (R vs S)\n'
        '  └ Survival: 227 (79 events)\n'
        'External: 311'
    )
    rounded_box(ax, 6.8, 6.0, 1.8, 1.5, summary, fc='#FFF9C4', ec=C['border'], fontsize=7, lw=1)

    plt.tight_layout()
    return fig


if __name__ == '__main__':
    print('Creating Figure 1 v2...')

    fig = create_figure()

    import os
    base = r'/path/to/xelox_project'
    hd_path = os.path.join(base, 'results', 'figures_png', 'hd', 'Fig1_pipeline.png')
    fig_path = os.path.join(base, 'results', 'figures_png', 'Fig1_pipeline.png')

    fig.savefig(hd_path, dpi=DPI, bbox_inches='tight', facecolor='white')
    fig.savefig(fig_path, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)

    print(f'  Saved: {hd_path}')
    print('Done!')
