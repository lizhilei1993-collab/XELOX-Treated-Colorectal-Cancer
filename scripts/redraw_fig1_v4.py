#!/usr/bin/env python3
"""
Figure 1 v4: Clean layout, unified font sizes, no overlap.
Each box has adequate padding. Simple vertical flow with two side branches.
"""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np

DPI = 300
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 10,
    'figure.dpi': DPI,
    'savefig.dpi': DPI,
    'savefig.facecolor': 'white',
    'figure.facecolor': 'white',
})

# Box style helpers
BW = 3.0   # standard box width
BH = 0.6   # standard box height

C = {
    'geo': '#E64B35',
    'qc': '#F9A825',
    'sel': '#F57C00',
    'gse': '#1565C0',
    'cla': '#00ACC1',
    'sur': '#43A047',
    'val': '#7B1FA2',
    'txt': '#212121',
    'mut': '#757575',
    'bg': '#FAFAFA',
}


def rbox(ax, x, y, w, h, title, subtitle='', color='#333', lw=1.2, ts=10, ss=8):
    """Draw a rounded box with title and optional subtitle."""
    box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.03",
                         facecolor='white', edgecolor=color, linewidth=lw)
    ax.add_patch(box)
    ax.text(x + w/2, y + h/2 + 0.05, title, ha='center', va='center',
            fontsize=ts, fontweight='bold', color=C['txt'])
    if subtitle:
        ax.text(x + w/2, y + h/2 - 0.25, subtitle, ha='center', va='center',
                fontsize=ss, color=C['mut'])


def rbox_filled(ax, x, y, w, h, title, subtitle='', fc='white', ec='#333',
                lw=1.2, ts=10, ss=8):
    """Rounded box with filled background."""
    box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.03",
                         facecolor=fc, edgecolor=ec, linewidth=lw)
    ax.add_patch(box)
    yc = y + h/2
    ax.text(x + w/2, yc + (0.05 if subtitle else 0), title,
            ha='center', va='center', fontsize=ts, fontweight='bold', color=C['txt'])
    if subtitle:
        ax.text(x + w/2, yc - 0.25, subtitle, ha='center', va='center',
                fontsize=ss, color=C['mut'])


def arrow(ax, x, y1, y2, lw=1.5):
    """Vertical arrow."""
    ax.annotate('', xy=(x, y2), xytext=(x, y1),
                arrowprops=dict(arrowstyle='->', color='#555', lw=lw))


def count_label(ax, x, y, text, color='#555', fs=8):
    """Count label next to a box."""
    ax.text(x, y, text, fontsize=fs, color=color, ha='left', va='center')


def main():
    fig, ax = plt.subplots(figsize=(7, 9.5))
    ax.set_xlim(0, 7)
    ax.set_ylim(0, 9.5)
    ax.axis('off')

    CX = 2  # center x

    # ── Title ──
    ax.text(3.5, 9.2, 'Study Design and Analytical Pipeline', fontsize=13,
            fontweight='bold', ha='center', color='#1a1a2e')
    ax.text(3.5, 8.9, 'TRIPOD-compliant reporting for transparent model development and validation',
            fontsize=8, ha='center', color=C['mut'], style='italic')

    # ── Row 1: GEO Retrieval ──
    y1 = 8.0
    rbox_filled(ax, CX, y1, BW, BH, 'GEO Data Retrieval',
                '5 Discovery Cohorts', fc='#FFEBEE', ec=C['geo'], ts=11)
    count_label(ax, CX + BW + 0.15, y1 + BH/2, 'n = 1,015', C['geo'], 9)

    # ── Row 2: Quality Control ──
    y2 = 6.9
    arrow(ax, CX + BW/2, y1, y2 + BH)
    rbox_filled(ax, CX, y2, BW, BH, 'Quality Control & Preprocessing', '',
                fc='#FFF8E1', ec=C['qc'], ts=11)
    count_label(ax, CX + BW + 0.15, y2 + BH/2, 'Excluded: 539', '#999', 8)
    ax.text(CX + BW + 0.15, y2 + BH/2 - 0.2, '(non-CRC, no chemo info)', fontsize=6, color='#bbb', ha='left')

    # ── Row 3: XELOX Subgroup ──
    y3 = 5.8
    arrow(ax, CX + BW/2, y2, y3 + BH)
    rbox_filled(ax, CX, y3, BW, BH, 'XELOX / FOLFOX-treated mCRC',
                '5 Cohorts Combined', fc='#FFF3E0', ec=C['sel'], ts=11)
    count_label(ax, CX + BW + 0.15, y3 + BH/2, 'n = 476', C['sel'], 9)

    # ── Row 4: GSE39582 Focus ──
    y4 = 4.7
    arrow(ax, CX + BW/2, y3, y4 + BH)
    rbox_filled(ax, CX, y4, BW, BH, 'GSE39582 Discovery Cohort',
                'Primary Analysis', fc='#E3F2FD', ec=C['gse'], ts=11, lw=2)
    count_label(ax, CX + BW + 0.15, y4 + BH/2, 'n = 239  (79 events)', C['gse'], 9)

    # ── Row 5: Two Paths ──
    y5 = 3.3
    arrow(ax, CX + BW/2, y4, y5 + BH)

    # Left: Classification
    rbox_filled(ax, 0.3, y5, 2.8, BH, 'Binary Classification',
                'Resistant (45) vs Sensitive (119)', fc='#E0F7FA', ec=C['cla'], ts=10, ss=8)
    ax.text(3.4, y5 + BH/2 + 0.1, '[EX]', fontsize=8, color=C['geo'],
            fontweight='bold', ha='center', va='center')
    ax.text(3.5, y5 + BH/2 + 0.1, 'Excluded: 75', fontsize=8, color='#999',
            ha='left', va='center')
    ax.text(3.5, y5 + BH/2 - 0.2, '(Intermediate RFS)', fontsize=6.5, color='#bbb',
            ha='left', va='center')

    # Right: Survival
    rbox_filled(ax, 3.9, y5, 2.8, BH, 'Survival Analysis',
                'RFS Endpoint, 79 Events', fc='#E8F5E9', ec=C['sur'], ts=10, ss=8)
    ax.text(7.0, y5 + BH/2 + 0.1, '[IN]', fontsize=8, color=C['sur'],
            fontweight='bold', ha='center', va='center')
    ax.text(7.1, y5 + BH/2 + 0.1, 'Included: 66', fontsize=8, color=C['sur'],
            ha='left', va='center')
    ax.text(7.1, y5 + BH/2 - 0.2, '(Intermediate RFS)', fontsize=6.5, color='#999',
            ha='left', va='center')

    # ── Row 6: Methods ──
    y6 = 2.2
    arrow(ax, 0.3 + 1.4, y5, y6 + BH)
    arrow(ax, 3.9 + 1.4, y5, y6 + BH)

    rbox_filled(ax, 0.3, y6, 2.8, BH, 'XGBoost + SHAP',
                '10-Gene Fingerprint (Table 1)', fc='#B2EBF2', ec=C['cla'], ts=10, ss=8)

    rbox_filled(ax, 3.9, y6, 2.8, BH, 'ssGSEA + AIC-Cox',
                '8-Variable Nomogram (Table 2)', fc='#C8E6C9', ec=C['sur'], ts=10, ss=8)

    # ── Row 7: External Validation ──
    y7 = 1.0
    arrow(ax, 3.9 + 1.4, y6, y7 + BH)
    rbox_filled(ax, 0.5, y7, 6.0, BH, 'External Validation Cohorts',
                'TCGA-COAD/READ (n=140)  |  GSE104645 (n=113)  |  GSE72970 (n=32)  |  GSE83129 (n=26)',
                fc='#F3E5F5', ec=C['val'], ts=11, ss=8)
    count_label(ax, 0.5, y7 - 0.15, 'Total: n = 311', C['val'], 8)

    # ── Summary Box ──
    sx, sy = 5.2, 5.8
    sh = 1.7
    sw = 1.6
    summary = FancyBboxPatch((sx, sy), sw, sh, boxstyle="round,pad=0.03",
                             facecolor='#FFF8E1', edgecolor='#F9A825', linewidth=1.5)
    ax.add_patch(summary)

    lines = [
        ('Sample Flow', 10, 'bold', sy + sh - 0.15),
        ('───────────', 7, 'normal', sy + sh - 0.35),
        ('Pool:  1,015', 8, 'normal', sy + sh - 0.55),
        ('QC:    476 (47%)', 8, 'normal', sy + sh - 0.75),
        ('GSE39582:  239', 8, 'bold', sy + sh - 0.95),
        ('   Class: 164', 7, 'normal', sy + sh - 1.12),
        ('   Surv: 227', 7, 'normal', sy + sh - 1.27),
        ('External:  311', 8, 'bold', sy + sh - 1.47),
    ]
    for text, fs, weight, ly in lines:
        ax.text(sx + sw/2, ly, text, ha='center', va='center',
                fontsize=fs, fontweight=weight, color=C['txt'])

    # Arrow to summary
    ax.annotate('', xy=(sx, sy + sh/2), xytext=(CX + BW, sy + sh/2),
                arrowprops=dict(arrowstyle='->', color='#999', lw=1, linestyle='dashed'))

    plt.tight_layout()
    return fig


if __name__ == '__main__':
    print('Creating Figure 1 v4...')
    fig = main()
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
