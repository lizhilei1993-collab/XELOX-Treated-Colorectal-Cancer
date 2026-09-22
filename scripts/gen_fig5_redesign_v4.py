"""
Figure 5 Redesign v4 — Academic clarity: no overlapping labels
Story: pathway aggregation prevents overfitting collapse
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

# ── Style ──
GENE_COLOR = '#C94C4C'      # muted red
PATH_COLOR = '#5B9BD5'      # muted blue
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
    'savefig.pad_inches': 0.10,
})

# ── Data ──
apparent = [0.500, 0.500]
corrected = [0.406, 0.489]
optimism = [0.0945, 0.011]
epv = [1.1, 11.3]

BASE = r'/path/to/xelox_project'
OUT = os.path.join(BASE, 'results', 'figures_png', 'v5.11', 'Fig5_fair_comparison_v4.png')


def add_panel_label(ax, label, x=-0.12, y=1.04):
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=14, fontweight='bold', va='top', ha='left',
            color='#2C3E50')


def fig5_v4():
    fig = plt.figure(figsize=(7.4, 3.8))
    gs = fig.add_gridspec(1, 2, wspace=0.52, left=0.10, right=0.96,
                          top=0.82, bottom=0.16)

    # ═══════════════════════════════════════════
    # Panel A: Slope chart — with clear non-overlapping labels
    # ═══════════════════════════════════════════
    ax = fig.add_subplot(gs[0, 0])
    ax.set_facecolor(BG_COLOR)

    x_app, x_cor = 0, 1
    colors = [GENE_COLOR, PATH_COLOR]

    # Draw lines first (behind dots)
    for i, c in enumerate(colors):
        ax.plot([x_app, x_cor], [apparent[i], corrected[i]],
                color=c, linewidth=2.5, alpha=0.85, zorder=2,
                solid_capstyle='round')

    # Draw dots on top
    for i, c in enumerate(colors):
        ax.scatter([x_app], [apparent[i]], s=90, color=c,
                   edgecolors='white', linewidths=1.2, zorder=4)
        ax.scatter([x_cor], [corrected[i]], s=90, color=c,
                   edgecolors='white', linewidths=1.2, zorder=4)

    # ── Left side: Apparent value (one label, centered between dots) ──
    ax.text(x_app - 0.10, apparent[0] + 0.006, '0.500',
            ha='right', va='bottom', fontsize=8, color='#444444', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                      edgecolor='none', alpha=0.9))

    # ── Right side: Corrected values + EPV (staggered vertically) ──
    # Gene: below the dot
    ax.text(x_cor + 0.10, corrected[0] - 0.008, f'{corrected[0]:.3f}',
            ha='left', va='top', fontsize=9, color=GENE_COLOR, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                      edgecolor='none', alpha=0.9))
    # Gene EPV: further below
    ax.text(x_cor + 0.10, corrected[0] - 0.030, f'EPV = {epv[0]:.1f}',
            ha='left', va='top', fontsize=6.5, color=GENE_COLOR, alpha=0.85, style='italic')

    # Pathway: above the dot
    ax.text(x_cor + 0.10, corrected[1] + 0.008, f'{corrected[1]:.3f}',
            ha='left', va='bottom', fontsize=9, color=PATH_COLOR, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                      edgecolor='none', alpha=0.9))
    # Pathway EPV: further above
    ax.text(x_cor + 0.10, corrected[1] + 0.030, f'EPV = {epv[1]:.1f}',
            ha='left', va='bottom', fontsize=6.5, color=PATH_COLOR, alpha=0.85, style='italic')

    # ── Random reference line ──
    ax.axhline(0.500, color='#BBBBBB', linestyle='--', linewidth=0.9, alpha=0.7, zorder=1)
    # Random label: far right, below the line to avoid EPV=11.3
    ax.text(1.38, 0.485, 'Random\n(C = 0.500)',
            fontsize=6, color='#999999', va='top', ha='left',
            linespacing=0.85,
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                      edgecolor='#CCCCCC', alpha=0.95))

    # ── Optimism drop annotations (clearly separated from lines) ──
    # Gene: annotation on the LEFT side, well away from the line
    mid_gene = (apparent[0] + corrected[0]) / 2
    ax.annotate('', xy=(0.28, corrected[0] + 0.015), xytext=(0.28, apparent[0] - 0.015),
                arrowprops=dict(arrowstyle='<->', color=GENE_COLOR, lw=1.0))
    ax.text(0.10, mid_gene, f'$\\Delta$ = {optimism[0]:.3f}',
            ha='left', va='center', fontsize=7.5, color=GENE_COLOR, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                      edgecolor='none', alpha=0.95))

    # Pathway: annotation on the RIGHT side, well away from the line
    mid_path = (apparent[1] + corrected[1]) / 2
    ax.annotate('', xy=(0.72, corrected[1] + 0.005), xytext=(0.72, apparent[1] - 0.005),
                arrowprops=dict(arrowstyle='<->', color=PATH_COLOR, lw=1.0))
    ax.text(0.90, mid_path, f'$\\Delta$ = {optimism[1]:.3f}',
            ha='right', va='center', fontsize=7.5, color=PATH_COLOR, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                      edgecolor='none', alpha=0.95))

    # ── Axis ──
    ax.set_xlim(-0.22, 1.55)
    ax.set_ylim(0.36, 0.54)
    ax.set_xticks([x_app, x_cor])
    ax.set_xticklabels(['Apparent\nC-index', 'Corrected\nC-index'], fontsize=8.5)
    ax.set_ylabel('C-index', fontsize=9.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('Overfitting collapse under identical LASSO-Cox',
                 fontsize=9.5, fontweight='bold', color='#2C3E50', pad=14)
    add_panel_label(ax, 'A')

    # Legend
    handles = [mpatches.Patch(color=GENE_COLOR, label='Gene-level (74 features)'),
               mpatches.Patch(color=PATH_COLOR, label='Pathway-level (7 features)')]
    ax.legend(handles=handles, fontsize=7, loc='lower left',
              framealpha=0.95, edgecolor='#CCCCCC', handlelength=1.3)

    # ═══════════════════════════════════════════
    # Panel B: Optimism comparison — clear labels
    # ═══════════════════════════════════════════
    ax = fig.add_subplot(gs[0, 1])
    ax.set_facecolor(BG_COLOR)

    y_pos = [1, 0]
    bar_height = 0.48
    opt_vals = [optimism[0], optimism[1]]

    bars = ax.barh(y_pos, opt_vals, height=bar_height, color=colors,
                   edgecolor='white', linewidth=0.6, alpha=0.85, zorder=3)

    # Value labels: gene inside (white), pathway outside (color) with background
    ax.text(opt_vals[0] - 0.004, y_pos[0], f'{opt_vals[0]:.4f}',
            va='center', ha='right', fontsize=9, color='white', fontweight='bold')
    ax.text(opt_vals[1] + 0.004, y_pos[1], f'{opt_vals[1]:.4f}',
            va='center', ha='left', fontsize=9, color=PATH_COLOR, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                      edgecolor='none', alpha=0.9))

    # ── 88% reduction bracket + label ──
    # Horizontal bracket between bar ends
    bracket_y_gene = 1.0
    bracket_y_path = 0.0
    bracket_x = 0.078
    ax.plot([bracket_x, bracket_x], [bracket_y_path + 0.18, bracket_y_gene - 0.18],
            color='#555555', lw=0.8)
    # Tick marks
    ax.plot([bracket_x - 0.002, bracket_x + 0.002], [bracket_y_path + 0.18, bracket_y_path + 0.18],
            color='#555555', lw=0.8)
    ax.plot([bracket_x - 0.002, bracket_x + 0.002], [bracket_y_gene - 0.18, bracket_y_gene - 0.18],
            color='#555555', lw=0.8)
    # Label
    ax.text(bracket_x + 0.006, 0.50, '88%\nreduction',
            ha='left', va='center', fontsize=7, fontweight='bold',
            color='#2C3E50', linespacing=0.9,
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                      edgecolor='#CCCCCC', alpha=0.95))

    # ── Axis ──
    ax.set_yticks(y_pos)
    ax.set_yticklabels(['Gene-level', 'Pathway-level'], fontsize=8.5)
    ax.set_xlabel('Bootstrap Optimism', fontsize=9.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(0, 0.11)
    ax.set_title('Bootstrap optimism (lower = less overfitting)',
                 fontsize=9.5, fontweight='bold', color='#2C3E50', pad=14)
    add_panel_label(ax, 'B')

    # ── Bottom caption ──
    fig.text(0.5, 0.01,
             'GSE39582 XELOX cohort (n = 227, 79 events)  ·  λ.1se, 10-fold CV  ·  B = 200 bootstrap',
             ha='center', fontsize=6, color='#888888', style='italic')

    plt.savefig(OUT, dpi=DPI, bbox_inches='tight', pad_inches=0.10)
    print(f'Saved: {OUT}')
    return fig


if __name__ == '__main__':
    fig5_v4()
    print('Done.')
