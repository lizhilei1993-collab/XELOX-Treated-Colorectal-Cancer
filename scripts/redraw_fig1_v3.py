#!/usr/bin/env python3
"""
Figure 1 v3: Clean, professional study pipeline flowchart.
Key improvements: balanced layout, consistent font sizes, clear hierarchy.
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np

DPI = 300
FONT = 'Arial'

plt.rcParams.update({
    'font.family': FONT,
    'font.size': 10,
    'figure.dpi': DPI,
    'savefig.dpi': DPI,
    'savefig.bbox': 'tight',
    'figure.facecolor': 'white',
})

BW = 0.8   # box width
BH = 0.4   # base box height
GAP = 0.6  # gap between rows

COLORS = {
    'title': '#1a1a2e',
    'geo': '#E64B35',
    'filter': '#F9A825',
    'split': '#1565C0',
    'class': '#00ACC1',
    'surv': '#43A047',
    'valid': '#7B1FA2',
    'summary': '#FFF9C4',
    'text': '#212121',
    'muted': '#757575',
    'line': '#555555',
    'bg': '#FAFAFA',
}


def box(ax, x, y, w, h, lines, fc='white', ec='#333', lw=1.2,
        fontsize=10, bold_first=False, align='center'):
    """Draw rounded box with centered text, supporting multi-line."""
    box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02",
                         facecolor=fc, edgecolor=ec, linewidth=lw, zorder=2)
    ax.add_patch(box)

    n = len(lines)
    for i, line in enumerate(lines):
        weight = 'bold' if (bold_first and i == 0) else 'normal'
        fs = fontsize - 1 if not bold_first else fontsize
        if bold_first and i == 0:
            fs = fontsize + 1
        va = 'center' if n == 1 else ('bottom' if i == 0 else 'top')
        offset = 0 if n == 1 else (BH * 0.2 * (n / 2 - i - 0.5))
        ax.text(x + w/2, y + h/2 + offset, line,
                ha='center', va='center', fontsize=fs, fontweight=weight,
                color=COLORS['text'])


def arrow_v(ax, x, y1, y2, color='#555', lw=1.5, label=''):
    """Vertical arrow from (x,y1) to (x,y2)."""
    ax.annotate('', xy=(x, y2), xytext=(x, y1 + 0.02),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw))
    if label:
        ax.text(x + 0.15, (y1 + y2) / 2, label, fontsize=7,
                color=COLORS['muted'], ha='left', va='center',
                style='italic')


def arrow_h(ax, x1, x2, y, color='#555', lw=1.5):
    """Horizontal arrow."""
    ax.annotate('', xy=(x2, y), xytext=(x1, y + 0.02),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw))


def side_label(ax, x, y, text, color='#888', fs=7, align='left'):
    """Small label next to a box."""
    ax.text(x, y, text, fontsize=fs, color=color, ha=align, va='center', style='italic')


def create():
    fig, ax = plt.subplots(figsize=(8, 10))
    ax.set_xlim(0, 8)
    ax.set_ylim(0, 10)
    ax.axis('off')

    # ── Title ──
    ax.text(4, 9.7, 'Study Design and Analytical Pipeline',
            fontsize=16, fontweight='bold', ha='center', color=COLORS['title'])
    ax.text(4, 9.4, 'TRIPOD-compliant reporting | XELOX Resistance Molecular Fingerprint',
            fontsize=9, ha='center', color=COLORS['muted'], style='italic')

    # ── 1. GEO Discovery ──
    y = 8.5
    box(ax, 3, y, 2, BH, ['GEO Data Retrieval', '5 Discovery Cohorts'],
        fc=COLORS['geo'], ec=COLORS['geo'], bold_first=True, fontsize=10)
    side_label(ax, 5.3, y + BH/2, 'n = 1,015', COLORS['geo'], fs=8)

    # ── 2. QC ──
    y2 = 7.4
    arrow_v(ax, 4, y - 0.02, y2 + BH)
    box(ax, 3, y2, 2, BH, ['Quality Control & Preprocessing'],
        fc='#FFE0B2', ec='#F9A825', fontsize=10)
    side_label(ax, 5.3, y2 + BH/2, 'Excluded: 539', '#999', fs=7,
               align='left')
    ax.text(5.5, y2 + BH/2 - 0.25, '(non-CRC, no chemo)', fontsize=6,
            color='#999', ha='left', va='center')

    # ── 3. XELOX Subgroup ──
    y3 = 6.3
    arrow_v(ax, 4, y2 - 0.02, y3 + BH)
    box(ax, 3, y3, 2, BH, ['XELOX / FOLFOX mCRC', 'Subgroup (5 Cohorts)'],
        fc='#FFCC80', ec='#F57C00', bold_first=True, fontsize=10)
    side_label(ax, 5.3, y3 + BH/2, 'n = 476', '#E65100', fs=8)

    # ── 4. Cohort details row ──
    y4 = 5.2
    arrow_v(ax, 4, y3 - 0.02, y4 + BH)
    cohorts = [
        ('GSE39582', 0.8, 'n=239', '#1565C0', 2.0),
        ('GSE104645', 0.6, 'n=193', '#666', 1.0),
        ('GSE28702', 0.6, 'n=83', '#666', 1.0),
        ('GSE72970', 0.6, 'n=124', '#666', 1.0),
        ('GSE69657', 0.6, 'n=30', '#666', 1.0),
    ]
    cx = 1.0
    for name, w, n, ec, lw in cohorts:
        box(ax, cx, y4, w, BH, [f'{name}', n], fc='white', ec=ec, lw=lw, fontsize=8)
        cx += w + 0.2

    # ── 5. GSE39582 Focus ──
    y5 = 4.0
    arrow_v(ax, 1.3, y4 - 0.02, y5 + BH, color=COLORS['split'])
    box(ax, 0.5, y5, 1.6, BH, ['GSE39582', 'Primary Cohort'],
        fc='#E3F2FD', ec=COLORS['split'], bold_first=True,
        fontsize=10, lw=2)
    side_label(ax, 0.5, y5 - 0.05, 'n = 239', COLORS['split'], fs=7, align='left')

    # ── 6. Two Analysis Paths ──

    # Path A: Classification (left)
    y6a = 2.8
    arrow_v(ax, 0.5 + 0.8, y5 - 0.02, y6a + BH, color=COLORS['class'])
    box(ax, 0.2, y6a, 2.2, 0.5,
        ['Binary Classification', 'Resistant (45) vs Sensitive (119)'],
        fc='#E0F7FA', ec=COLORS['class'], bold_first=True, fontsize=9)

    y7a = 1.9
    arrow_v(ax, 1.3, y6a - 0.02, y7a + 0.5)
    box(ax, 0.2, y7a, 2.2, 0.5,
        ['XGBoost + SHAP', '10-gene Fingerprint (Table 1)'],
        fc='#B2EBF2', ec=COLORS['class'], fontsize=9)

    # Intermediate exclusion label
    ax.text(2.7, y6a + 0.25, '✕ Excluded: 75', fontsize=8, color=COLORS['geo'],
            ha='left', va='center', fontweight='bold')
    ax.text(2.7, y6a - 0.05, '(Intermediate RFS)', fontsize=6.5, color='#999',
            ha='left', va='center')

    # Path B: Survival (right)
    y6b = 2.8
    arrow_v(ax, 3 + 1.0, y4 - 0.02, y6b + BH, color=COLORS['surv'])
    box(ax, 3.8, y6b, 2.2, 0.5,
        ['Survival Analysis', 'RFS endpoint, 79 events (34.8%)'],
        fc='#E8F5E9', ec=COLORS['surv'], bold_first=True, fontsize=9)

    y7b = 1.9
    arrow_v(ax, 4.9, y6b - 0.02, y7b + 0.5)
    box(ax, 3.8, y7b, 2.2, 0.5,
        ['ssGSEA + AIC-Cox', '8-variable Nomogram (Table 2)'],
        fc='#C8E6C9', ec=COLORS['surv'], fontsize=9)

    # Intermediate included label
    ax.text(6.3, y6b + 0.25, '✓ Included: 66', fontsize=8, color=COLORS['surv'],
            ha='left', va='center', fontweight='bold')
    ax.text(6.3, y6b - 0.05, '(Intermediate RFS)', fontsize=6.5, color='#999',
            ha='left', va='center')

    # ── 7. External Validation ──
    y8 = 0.7
    arrow_v(ax, 4.9, y7b - 0.02, y8 + 0.5, color=COLORS['valid'])
    box(ax, 0.3, y8, 3.4, 0.5,
        ['External Validation Cohorts'],
        fc='#F3E5F5', ec=COLORS['valid'], bold_first=True, fontsize=10)
    box(ax, 3.9, y8, 3.8, 0.5,
        ['TCGA (n=140) | GSE104645 (n=113) | GSE72970 (n=32) | GSE83129 (n=26)'],
        fc='#F3E5F5', ec=COLORS['valid'], fontsize=8)

    side_label(ax, 0.3, y8 - 0.05, 'n = 311 total', COLORS['valid'], fs=7, align='left')

    # ── Summary box ──
    sum_text = (
        'Sample Flow Summary\n'
        '─────────────────\n'
        'Discovery pool:  1,015\n'
        'After QC:        476 (47%)\n'
        'GSE39582 Focus:  239\n'
        '  Classification: 164\n'
        '  Survival:       227\n'
        'External valid:   311'
    )
    box(ax, 6.5, 5.5, 1.3, 2.0, [sum_text],
        fc=COLORS['summary'], ec='#F9A825', fontsize=7)

    # ── Arrow from GSE39582 to summary ──
    arrow_h(ax, 3.8, 6.5, 5.8, color='#999', lw=1)

    plt.tight_layout()
    return fig


if __name__ == '__main__':
    print('Creating Figure 1 v3...')
    fig = create()
    import os
    base = r'/path/to/xelox_project'
    for p in [
        os.path.join(base, 'results', 'figures_png', 'hd', 'Fig1_pipeline.png'),
        os.path.join(base, 'results', 'figures_png', 'Fig1_pipeline.png'),
    ]:
        fig.savefig(p, dpi=DPI, bbox_inches='tight', facecolor='white')
        print(f'  Saved: {p}')
    plt.close(fig)
    print('Done!')
