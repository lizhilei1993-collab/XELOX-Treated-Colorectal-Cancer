"""
Figure 5 Redesign — Reframe from "absolute C-index" to "overfitting resilience"
Story: pathway aggregation prevents overfitting collapse
Data: identical LASSO-Cox fair comparison (GSE39582, n=227, 79 events)
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

# ── Style ──
GENE_COLOR = '#E64B35'      # red (muted)
PATH_COLOR = '#4DBBD5'      # teal
BG_COLOR = '#FAFAFA'
GRID_COLOR = '#E0E0E0'
FONT_FAMILY = 'Arial'
DPI = 300

plt.rcParams.update({
    'font.family': FONT_FAMILY,
    'font.size': 8,
    'axes.linewidth': 0.6,
    'figure.dpi': DPI,
    'savefig.dpi': DPI,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
})

# ── Data (from fair comparison) ──
models = ['Gene-level\nLASSO-Cox', 'Pathway-level\nLASSO-Cox']
apparent = [0.500, 0.500]
corrected = [0.406, 0.489]
optimism = [0.0945, 0.011]
epv = [1.1, 11.3]
features = [74, 7]

BASE = r'/path/to/xelox_project'
OUT = os.path.join(BASE, 'results', 'figures_png', 'v5.11', 'Fig5_fair_comparison_v2.png')


def add_panel_label(ax, label, x=-0.12, y=1.05):
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=12, fontweight='bold', va='top', ha='left',
            color='#2C3E50')


def fig5_redesign():
    fig = plt.figure(figsize=(6.8, 3.2))
    gs = fig.add_gridspec(1, 2, wspace=0.45, left=0.10, right=0.96,
                          top=0.88, bottom=0.15)

    # ═══════════════════════════════════════════
    # Panel A: Slope chart — Apparent → Corrected
    # ═══════════════════════════════════════════
    ax = fig.add_subplot(gs[0, 0])
    ax.set_facecolor(BG_COLOR)

    x_app = 0  # Apparent
    x_cor = 1  # Corrected
    colors = [GENE_COLOR, PATH_COLOR]
    lw = 2.2

    for i, (c, label_short) in enumerate(zip(colors, ['Gene', 'Pathway'])):
        # Connecting line
        ax.plot([x_app, x_cor], [apparent[i], corrected[i]],
                color=c, linewidth=lw, alpha=0.9, zorder=3,
                solid_capstyle='round')
        # Dots
        ax.scatter([x_app], [apparent[i]], s=60, color=c,
                   edgecolors='white', linewidths=0.8, zorder=4)
        ax.scatter([x_cor], [corrected[i]], s=60, color=c,
                   edgecolors='white', linewidths=0.8, zorder=4)
        # Value labels
        ax.text(x_app - 0.06, apparent[i], f'{apparent[i]:.3f}',
                ha='right', va='center', fontsize=7, color=c, fontweight='bold')
        # For corrected, show corrected value + EPV annotation
        corrected_label = f'{corrected[i]:.3f}'
        ax.text(x_cor + 0.06, corrected[i], corrected_label,
                ha='left', va='center', fontsize=7, color=c, fontweight='bold')
        # EPV annotation near the corrected dot
        epv_y_offset = 0.025 if i == 0 else -0.025
        ax.text(x_cor + 0.06, corrected[i] + epv_y_offset,
                f'EPV={epv[i]:.1f}', ha='left', va='center',
                fontsize=5.5, color=c, alpha=0.8, style='italic')

    # 0.500 reference line
    ax.axhline(0.500, color='#999999', linestyle='--', linewidth=0.7,
               alpha=0.7, zorder=1)
    ax.text(1.02, 0.500, 'Random\n(C=0.500)', transform=ax.get_yaxis_transform(),
            fontsize=5, color='#999999', va='center', ha='left')

    # Collapse arrows (annotate the drop)
    # Gene drop
    mid_gene = (apparent[0] + corrected[0]) / 2
    ax.annotate('', xy=(0.38, corrected[0] + 0.015), xytext=(0.38, apparent[0] - 0.015),
                arrowprops=dict(arrowstyle='->', color=GENE_COLOR, lw=1.2,
                                connectionstyle='arc3,rad=0'))
    ax.text(0.42, mid_gene, f'-{optimism[0]:.3f}', fontsize=5.5,
            color=GENE_COLOR, va='center', fontweight='bold')

    # Pathway drop
    mid_path = (apparent[1] + corrected[1]) / 2
    ax.annotate('', xy=(0.62, corrected[1] + 0.008), xytext=(0.62, apparent[1] - 0.008),
                arrowprops=dict(arrowstyle='->', color=PATH_COLOR, lw=1.2,
                                connectionstyle='arc3,rad=0'))
    ax.text(0.66, mid_path, f'-{optimism[1]:.3f}', fontsize=5.5,
            color=PATH_COLOR, va='center', fontweight='bold')

    ax.set_xlim(-0.15, 1.35)
    ax.set_ylim(0.35, 0.55)
    ax.set_xticks([x_app, x_cor])
    ax.set_xticklabels(['Apparent\nC-index', 'Corrected\nC-index'], fontsize=7.5)
    ax.set_ylabel('C-index', fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('Overfitting collapse under identical\nLASSO-Cox methodology',
                 fontsize=8.5, fontweight='bold', color='#2C3E50', pad=6)
    add_panel_label(ax, 'A')

    # Legend
    handles = [mpatches.Patch(color=GENE_COLOR, label='Gene-level (74 features)'),
               mpatches.Patch(color=PATH_COLOR, label='Pathway-level (7 features)')]
    ax.legend(handles=handles, fontsize=6, loc='lower left',
              framealpha=0.9, edgecolor=GRID_COLOR, handlelength=1.2)

    # ═══════════════════════════════════════════
    # Panel B: Optimism comparison (horizontal bars)
    # ═══════════════════════════════════════════
    ax = fig.add_subplot(gs[0, 1])
    ax.set_facecolor(BG_COLOR)

    y_pos = [1, 0]
    bar_height = 0.55
    opt_vals = [optimism[0], optimism[1]]  # Gene first, then Pathway

    bars = ax.barh(y_pos, opt_vals, height=bar_height, color=colors,
                   edgecolor='white', linewidth=0.6, alpha=0.85, zorder=3)

    # Value labels
    for i, (val, yp) in enumerate(zip(opt_vals, y_pos)):
        ax.text(val + 0.003, yp, f'{val:.4f}', va='center', ha='left',
                fontsize=7, fontweight='bold', color=colors[i])

    # 88% reduction annotation
    reduction_pct = (1 - optimism[1] / optimism[0]) * 100
    ax.annotate(f'88% reduction\nin overfitting',
                xy=(optimism[1], 0), xytext=(0.06, -0.45),
                fontsize=6.5, fontweight='bold', color='#2C3E50',
                ha='center', va='top',
                arrowprops=dict(arrowstyle='->', color='#2C3E50', lw=0.8,
                                connectionstyle='arc3,rad=-0.2'))

    ax.set_yticks(y_pos)
    ax.set_yticklabels(['Gene-level', 'Pathway-level'], fontsize=7.5)
    ax.set_xlabel('Bootstrap Optimism', fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(0, 0.12)
    ax.set_title('Bootstrap optimism\n(lower = less overfitting)',
                 fontsize=8.5, fontweight='bold', color='#2C3E50', pad=6)
    add_panel_label(ax, 'B')

    # Bottom annotation
    fig.text(0.5, 0.01,
             'GSE39582 XELOX cohort (n=227, 79 events) · λ.1se, 10-fold CV · B=200 bootstrap',
             ha='center', fontsize=5.5, color='#888888', style='italic')

    plt.savefig(OUT, dpi=DPI, bbox_inches='tight', pad_inches=0.05)
    print(f'Saved: {OUT}')
    return fig


if __name__ == '__main__':
    fig5_redesign()
    print('Done.')
