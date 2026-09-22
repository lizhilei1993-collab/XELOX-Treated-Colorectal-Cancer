"""
Figure 5 Redesign v3 — Fix layout issues
Story: pathway aggregation prevents overfitting collapse
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
import numpy as np
import os

# ── Style ──
GENE_COLOR = '#D9534F'      # softer red
PATH_COLOR = '#5BC0DE'      # softer teal
BG_COLOR = '#FAFAFA'
FONT_FAMILY = 'Arial'
DPI = 300

plt.rcParams.update({
    'font.family': FONT_FAMILY,
    'font.size': 8,
    'axes.linewidth': 0.6,
    'figure.dpi': DPI,
    'savefig.dpi': DPI,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.08,
})

# ── Data ──
apparent = [0.500, 0.500]
corrected = [0.406, 0.489]
optimism = [0.0945, 0.011]
epv = [1.1, 11.3]

BASE = r'/path/to/xelox_project'
OUT = os.path.join(BASE, 'results', 'figures_png', 'v5.11', 'Fig5_fair_comparison_v3.png')


def add_panel_label(ax, label, x=-0.14, y=1.02):
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=13, fontweight='bold', va='top', ha='left',
            color='#2C3E50')


def fig5_v3():
    fig = plt.figure(figsize=(7.0, 3.4))
    gs = fig.add_gridspec(1, 2, wspace=0.50, left=0.10, right=0.96,
                          top=0.86, bottom=0.18)

    # ═══════════════════════════════════════════
    # Panel A: Slope chart
    # ═══════════════════════════════════════════
    ax = fig.add_subplot(gs[0, 0])
    ax.set_facecolor(BG_COLOR)

    x_app, x_cor = 0, 1
    colors = [GENE_COLOR, PATH_COLOR]

    # Draw lines and dots
    for i, c in enumerate(colors):
        ax.plot([x_app, x_cor], [apparent[i], corrected[i]],
                color=c, linewidth=2.5, alpha=0.9, zorder=3,
                solid_capstyle='round')
        ax.scatter([x_app, x_cor], [apparent[i], corrected[i]],
                   s=70, color=c, edgecolors='white', linewidths=1.0, zorder=4)

    # ── Labels: Apparent (left side, staggered) ──
    # Only one apparent label since both are 0.500
    ax.text(x_app - 0.08, apparent[0] + 0.008, '0.500',
            ha='right', va='bottom', fontsize=7.5, color='#555555', fontweight='bold')

    # ── Labels: Corrected (right side) ──
    ax.text(x_cor + 0.08, corrected[0] - 0.005, f'{corrected[0]:.3f}',
            ha='left', va='top', fontsize=8, color=GENE_COLOR, fontweight='bold')
    ax.text(x_cor + 0.08, corrected[1] + 0.005, f'{corrected[1]:.3f}',
            ha='left', va='bottom', fontsize=8, color=PATH_COLOR, fontweight='bold')

    # ── EPV annotations (right side, below/above corrected values) ──
    ax.text(x_cor + 0.08, corrected[0] - 0.022, f'EPV = {epv[0]:.1f}',
            ha='left', va='top', fontsize=5.5, color=GENE_COLOR, alpha=0.75, style='italic')
    ax.text(x_cor + 0.08, corrected[1] + 0.022, f'EPV = {epv[1]:.1f}',
            ha='left', va='bottom', fontsize=5.5, color=PATH_COLOR, alpha=0.75, style='italic')

    # ── Random reference line + label ──
    ax.axhline(0.500, color='#AAAAAA', linestyle='--', linewidth=0.8, alpha=0.6, zorder=1)
    # Place "Random" label on the right, aligned with the line, below EPV text
    ax.text(1.32, 0.495, 'Random\n(C = 0.500)',
            fontsize=5.5, color='#999999', va='top', ha='left',
            linespacing=0.9)

    # ── Optimism drop annotations (vertical arrows in middle) ──
    # Gene: big drop
    ax.annotate('', xy=(0.45, corrected[0] + 0.01), xytext=(0.45, apparent[0] - 0.01),
                arrowprops=dict(arrowstyle='<->', color=GENE_COLOR, lw=1.0,
                                connectionstyle='arc3,rad=0'))
    ax.text(0.50, (apparent[0] + corrected[0]) / 2, f'$\\Delta$ = {optimism[0]:.3f}',
            ha='left', va='center', fontsize=6, color=GENE_COLOR, fontweight='bold')

    # Pathway: small drop
    ax.annotate('', xy=(0.55, corrected[1] - 0.005), xytext=(0.55, apparent[1] + 0.005),
                arrowprops=dict(arrowstyle='<->', color=PATH_COLOR, lw=1.0,
                                connectionstyle='arc3,rad=0'))
    ax.text(0.60, (apparent[1] + corrected[1]) / 2, f'$\\Delta$ = {optimism[1]:.3f}',
            ha='left', va='center', fontsize=6, color=PATH_COLOR, fontweight='bold')

    # ── Axis settings ──
    ax.set_xlim(-0.20, 1.55)
    ax.set_ylim(0.37, 0.535)
    ax.set_xticks([x_app, x_cor])
    ax.set_xticklabels(['Apparent\nC-index', 'Corrected\nC-index'], fontsize=8)
    ax.set_ylabel('C-index', fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('Overfitting collapse under identical LASSO-Cox',
                 fontsize=9, fontweight='bold', color='#2C3E50', pad=8)
    add_panel_label(ax, 'A')

    # Legend
    handles = [mpatches.Patch(color=GENE_COLOR, label='Gene-level (74 features)'),
               mpatches.Patch(color=PATH_COLOR, label='Pathway-level (7 features)')]
    ax.legend(handles=handles, fontsize=6.5, loc='lower left',
              framealpha=0.9, edgecolor='#CCCCCC', handlelength=1.2)

    # ═══════════════════════════════════════════
    # Panel B: Optimism comparison
    # ═══════════════════════════════════════════
    ax = fig.add_subplot(gs[0, 1])
    ax.set_facecolor(BG_COLOR)

    y_pos = [1, 0]
    bar_height = 0.50
    opt_vals = [optimism[0], optimism[1]]

    bars = ax.barh(y_pos, opt_vals, height=bar_height, color=colors,
                   edgecolor='white', linewidth=0.6, alpha=0.85, zorder=3)

    # Value labels (inside bars for gene, outside for pathway)
    ax.text(opt_vals[0] - 0.005, y_pos[0], f'{opt_vals[0]:.4f}',
            va='center', ha='right', fontsize=8, color='white', fontweight='bold')
    ax.text(opt_vals[1] + 0.003, y_pos[1], f'{opt_vals[1]:.4f}',
            va='center', ha='left', fontsize=8, color=PATH_COLOR, fontweight='bold')

    # ── 88% reduction bracket annotation ──
    # Draw a curly brace between the two bars
    brace_x = 0.075
    ax.annotate('', xy=(brace_x, 0.25), xytext=(brace_x, 0.75),
                arrowprops=dict(arrowstyle='-', color='#555555', lw=0.8))
    ax.plot([brace_x - 0.003, brace_x + 0.003], [0.25, 0.25], color='#555555', lw=0.8)
    ax.plot([brace_x - 0.003, brace_x + 0.003], [0.75, 0.75], color='#555555', lw=0.8)
    ax.text(brace_x + 0.008, 0.50, '88%\nreduction',
            ha='left', va='center', fontsize=6.5, fontweight='bold',
            color='#2C3E50', linespacing=0.9)

    # ── Axis settings ──
    ax.set_yticks(y_pos)
    ax.set_yticklabels(['Gene-level', 'Pathway-level'], fontsize=8)
    ax.set_xlabel('Bootstrap Optimism', fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(0, 0.11)
    ax.set_title('Bootstrap optimism (lower = less overfitting)',
                 fontsize=9, fontweight='bold', color='#2C3E50', pad=8)
    add_panel_label(ax, 'B')

    # ── Bottom caption ──
    fig.text(0.5, 0.02,
             'GSE39582 XELOX cohort (n = 227, 79 events)  ·  λ.1se, 10-fold CV  ·  B = 200 bootstrap',
             ha='center', fontsize=5.5, color='#888888', style='italic')

    plt.savefig(OUT, dpi=DPI, bbox_inches='tight', pad_inches=0.08)
    print(f'Saved: {OUT}')
    return fig


if __name__ == '__main__':
    fig5_v3()
    print('Done.')
