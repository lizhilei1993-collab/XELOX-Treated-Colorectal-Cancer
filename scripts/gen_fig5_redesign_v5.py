"""
Figure 5 Redesign v5 — Final clean version
Changes:
  1. Panel A: Δ labels on line segments, no absolute corrected values
  2. Panel B: two-row layout — top: Optimism, bottom: EPV
  3. Remove "Random" dashed line, rename to subtle "Baseline"
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

# ── Style ──
GENE_COLOR = '#E64B35'
PATH_COLOR = '#4DBBD5'
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
OUT = os.path.join(BASE, 'results', 'figures_png', 'v5.11', 'Fig5_fair_comparison_v5.png')


def add_panel_label(ax, label, x=-0.14, y=1.03):
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=13, fontweight='bold', va='top', ha='left',
            color='#2C3E50')


def fig5_v5():
    fig = plt.figure(figsize=(7.2, 4.0))
    gs = fig.add_gridspec(1, 2, wspace=0.48, left=0.10, right=0.96,
                          top=0.88, bottom=0.12)

    # ═══════════════════════════════════════════
    # Panel A: Slope chart — Δ labels on segments
    # ═══════════════════════════════════════════
    ax = fig.add_subplot(gs[0, 0])
    ax.set_facecolor(BG_COLOR)

    x_app, x_cor = 0, 1
    colors = [GENE_COLOR, PATH_COLOR]

    # ── Connecting lines ──
    for i, c in enumerate(colors):
        ax.plot([x_app, x_cor], [apparent[i], corrected[i]],
                color=c, linewidth=2.4, alpha=0.85, zorder=2,
                solid_capstyle='round')

    # ── Dots (triangles) ──
    for i, c in enumerate(colors):
        ax.scatter([x_app], [apparent[i]], s=90, color=c, marker='^',
                   edgecolors='white', linewidths=1.0, zorder=4)
        ax.scatter([x_cor], [corrected[i]], s=90, color=c, marker='v',
                   edgecolors='white', linewidths=1.0, zorder=4)

    # ── Left: Apparent label (both = 0.500, single label) ──
    ax.text(x_app - 0.06, 0.500, '0.500',
            ha='right', va='center', fontsize=8.5, color='#444444', fontweight='bold')

    # ── Δ labels along the line segments ──
    # Gene: big drop, label on lower portion of segment, slightly left-down
    gene_line_y_at_07 = apparent[0] + (corrected[0] - apparent[0]) * 0.72
    ax.text(0.66, gene_line_y_at_07 + 0.006,
            r'$\Delta$ = $-0.095$',
            ha='center', va='bottom', fontsize=8, color=GENE_COLOR, fontweight='bold',
            rotation=-55,  # follow the line angle
            )

    # Pathway: small drop, label on right portion of segment
    path_line_y_at_07 = apparent[1] + (corrected[1] - apparent[1]) * 0.7
    ax.text(0.72, path_line_y_at_07 - 0.008,
            r'$\Delta$ = $-0.011$',
            ha='center', va='top', fontsize=8, color=PATH_COLOR, fontweight='bold',
            rotation=-5,  # nearly horizontal
            )

    # ── Axis ──
    ax.set_xlim(-0.22, 1.18)
    ax.set_ylim(0.375, 0.53)
    ax.set_xticks([x_app, x_cor])
    ax.set_xticklabels(['Apparent\nC-index', 'Corrected\nC-index'], fontsize=8)
    ax.set_ylabel('C-index', fontsize=9)
    # Explicit y-ticks: skip 0.500 (we label it manually), keep clean spacing
    ax.set_yticks([0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 0.52])
    ax.set_yticklabels(['0.38', '0.40', '0.42', '0.44', '0.46', '0.48', '0.52'])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('Overfitting collapse under identical LASSO-Cox',
                 fontsize=9, fontweight='bold', color='#2C3E50', pad=8)
    add_panel_label(ax, 'A')

    # Legend
    handles = [mpatches.Patch(color=GENE_COLOR, label='Gene-level (74 features)'),
               mpatches.Patch(color=PATH_COLOR, label='Pathway-level (7 features)')]
    ax.legend(handles=handles, fontsize=6.5, loc='lower left',
              framealpha=0.95, edgecolor='#CCCCCC', handlelength=1.2)

    # ═══════════════════════════════════════════
    # Panel B: Two-row comparison
    # ═══════════════════════════════════════════
    gs_b = gs[0, 1].subgridspec(2, 1, hspace=0.45)

    # ── Row 1: Optimism (Pathway wins) ──
    ax1 = fig.add_subplot(gs_b[0])
    ax1.set_facecolor(BG_COLOR)

    y_pos = [1, 0]
    bar_height = 0.50
    opt_vals = [optimism[0], optimism[1]]

    ax1.barh(y_pos, opt_vals, height=bar_height, color=colors,
             edgecolor='white', linewidth=0.6, alpha=0.85, zorder=3)

    # Value labels
    ax1.text(opt_vals[0] + 0.003, y_pos[0], f'{opt_vals[0]:.4f}',
             va='center', ha='left', fontsize=8, color=GENE_COLOR, fontweight='bold')
    ax1.text(opt_vals[1] + 0.003, y_pos[1], f'{opt_vals[1]:.4f}',
             va='center', ha='left', fontsize=8, color=PATH_COLOR, fontweight='bold')

    # 88% reduction bracket
    bracket_x = 0.005
    ax1.annotate('', xy=(bracket_x, 0.25), xytext=(bracket_x, 0.75),
                 arrowprops=dict(arrowstyle='<->', color='#555555', lw=0.9,
                                 mutation_scale=10))
    ax1.text(bracket_x + 0.005, 0.50, '88%\nreduction',
             ha='left', va='center', fontsize=7, fontweight='bold',
             color='#2C3E50', linespacing=0.9)

    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(['Gene', 'Pathway'], fontsize=8)
    ax1.set_xlabel('Bootstrap Optimism', fontsize=8)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.set_xlim(0, 0.10)
    ax1.set_title('Overfitting (lower = better)',
                  fontsize=8.5, fontweight='bold', color='#2C3E50', pad=6, loc='center')

    # ── Row 2: EPV (Pathway big win) ──
    ax2 = fig.add_subplot(gs_b[1])
    ax2.set_facecolor(BG_COLOR)

    epv_vals = [epv[0], epv[1]]

    ax2.barh(y_pos, epv_vals, height=bar_height, color=colors,
             edgecolor='white', linewidth=0.6, alpha=0.85, zorder=3)

    # Value labels
    ax2.text(epv_vals[0] + 0.3, y_pos[0], f'{epv_vals[0]:.1f}',
             va='center', ha='left', fontsize=8, color=GENE_COLOR, fontweight='bold')
    ax2.text(epv_vals[1] + 0.3, y_pos[1], f'{epv_vals[1]:.1f}',
             va='center', ha='left', fontsize=8, color=PATH_COLOR, fontweight='bold')

    # EPV threshold line at 10
    ax2.axvline(10, color='#CC8844', linestyle='--', linewidth=0.8, alpha=0.7, zorder=1)
    ax2.text(10.2, 0.55, 'EPV ≥ 10\n(reliable)',
             fontsize=6, color='#CC8844', va='bottom', ha='left', style='italic')

    # 10x annotation (arrow: Gene → Pathway, with tick marks)
    arrow_y = 0.25
    ax2.annotate('', xy=(epv_vals[1], arrow_y), xytext=(epv_vals[0], arrow_y),
                 arrowprops=dict(arrowstyle='->', color='#555555', lw=1.0,
                                 mutation_scale=12))
    # Tick marks at both ends
    tick_half = 0.08
    ax2.plot([epv_vals[0], epv_vals[0]], [arrow_y - tick_half, arrow_y + tick_half],
             color='#555555', lw=0.9, zorder=5)
    ax2.plot([epv_vals[1], epv_vals[1]], [arrow_y - tick_half, arrow_y + tick_half],
             color='#555555', lw=0.9, zorder=5)
    ax2.text((epv_vals[0] + epv_vals[1]) / 2, 0.38, '10× higher',
             ha='center', va='center', fontsize=7, fontweight='bold',
             color='#2C3E50')

    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(['Gene', 'Pathway'], fontsize=8)
    ax2.set_xlabel('Events Per Variable (EPV)', fontsize=8)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.set_xlim(0, 15)
    ax2.set_title('Statistical power (higher = better)',
                  fontsize=8.5, fontweight='bold', color='#2C3E50', pad=6, loc='center')

    # Panel B label on top
    add_panel_label(ax1, 'B')

    # ── Bottom caption ──
    fig.text(0.5, 0.01,
             'GSE39582 XELOX cohort (n = 227, 79 events)  ·  $\\lambda$.1se, 10-fold CV  ·  B = 200 bootstrap',
             ha='center', fontsize=6.5, color='#555555', style='italic')

    plt.savefig(OUT, dpi=DPI, bbox_inches='tight', pad_inches=0.08)
    print(f'Saved: {OUT}')
    return fig


if __name__ == '__main__':
    fig5_v5()
    print('Done.')
