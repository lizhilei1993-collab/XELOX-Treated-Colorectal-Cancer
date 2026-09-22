#!/usr/bin/env python3
"""
Figure 1 v5: Publication-quality study pipeline.
Outputs both SVG (vector) and high-res PNG.
Based on redraw_fig1_v4.py with improved layout and verified connections.
Journal of Translational Medicine: 170mm width, 300 DPI.
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patheffects as pe

DPI = 300
BASE = r'/path/to/xelox_project'
OUT_PNG = os.path.join(BASE, 'results', 'figures_png', 'v5.11')
OUT_SVG = os.path.join(BASE, 'results', 'figures_png', 'v5.11')
os.makedirs(OUT_PNG, exist_ok=True)

# ── Color Palette (academic, muted) ──
C = {
    'txt': '#212121',
    'mut': '#666666',
    'light': '#999999',
    # Box fills & edges
    'pink_fc': '#FDECEA',    'pink_ec': '#C62828',
    'amber_fc': '#FFF8E1',   'amber_ec': '#F57F17',
    'blue_fc': '#E3F2FD',    'blue_ec': '#1565C0',
    'teal_fc': '#E0F7FA',    'teal_ec': '#00838F',
    'green_fc': '#E8F5E9',   'green_ec': '#2E7D32',
    'purple_fc': '#F3E5F5',  'purple_ec': '#6A1B9A',
    'lime_fc': '#F1F8E9',    'lime_ec': '#558B2F',
    'orange_fc': '#FFF3E0',  'orange_ec': '#E65100',
}

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 9,
    'figure.dpi': DPI,
    'savefig.dpi': DPI,
    'savefig.facecolor': 'white',
    'figure.facecolor': 'white',
})


def rbox(ax, x, y, w, h, title, subtitle='', fc='white', ec='#333', lw=1.2, ts=10, ss=7.5):
    """Rounded box with title and optional subtitle."""
    box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.03",
                         facecolor=fc, edgecolor=ec, linewidth=lw, zorder=2)
    ax.add_patch(box)
    if subtitle:
        ax.text(x + w/2, y + h/2 + 0.08, title, ha='center', va='center',
                fontsize=ts, fontweight='bold', color=C['txt'], zorder=3)
        ax.text(x + w/2, y + h/2 - 0.18, subtitle, ha='center', va='center',
                fontsize=ss, color=C['mut'], zorder=3)
    else:
        ax.text(x + w/2, y + h/2, title, ha='center', va='center',
                fontsize=ts, fontweight='bold', color=C['txt'], zorder=3)


def arrow_down(ax, x, y1, y2, lw=1.5, color='#555'):
    """Vertical arrow from (x, y1) to (x, y2)."""
    ax.annotate('', xy=(x, y2), xytext=(x, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw,
                                shrinkA=0, shrinkB=0),
                zorder=1)


def arrow_right(ax, x1, y1, x2, y2, lw=1.2, color='#999', ls='dashed'):
    """Horizontal or angled arrow."""
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw,
                                linestyle=ls, shrinkA=2, shrinkB=2),
                zorder=1)


def sidebar(ax, x, y, text, color='#555', fs=8, ha='left'):
    """Side label."""
    ax.text(x, y, text, fontsize=fs, color=color, ha=ha, va='center')


def main():
    fig, ax = plt.subplots(figsize=(7, 9.5))
    ax.set_xlim(0, 7)
    ax.set_ylim(0, 9.5)
    ax.axis('off')

    # ── Layout constants ──
    # Main column: x = 0.5 to 4.0 (width 3.5)
    # Right panel: x = 4.8 to 6.8 (width 2.0)
    CX = 0.5       # left edge of main boxes
    BW = 3.5       # box width
    BH = 0.55      # box height
    M = 0.5        # vertical margin between boxes

    # ── ROW 1: GEO Data Retrieval ──
    y1 = 8.5
    rbox(ax, CX, y1, BW, BH, 'GEO Data Retrieval',
         '5 Discovery Cohorts (GSE28702, GSE69665, GSE39582, GSE72970, GSE83129)',
         fc=C['pink_fc'], ec=C['pink_ec'], lw=1.5, ts=11, ss=7.5)
    sidebar(ax, CX + BW + 0.15, y1 + BH/2, 'n = 1,015', C['pink_ec'], 9)

    # ── ROW 2: Quality Control ──
    y2 = y1 - BH - M
    arrow_down(ax, CX + BW/2, y1, y2 + BH)
    rbox(ax, CX, y2, BW, BH, 'Quality Control & Preprocessing',
         'Platform normalization, probe filtering (IQR-based), batch correction (ComBat)',
         fc=C['amber_fc'], ec=C['amber_ec'], lw=1.5, ts=10, ss=7.5)
    sidebar(ax, CX + BW + 0.15, y2 + BH/2, 'Excluded: n = 539', C['light'], 8)
    sidebar(ax, CX + BW + 0.15, y2 + BH/2 - 0.22, '(non-CRC, no chemo info)', '#bbb', 6.5)

    # ── ROW 3: XELOX mCRC ──
    y3 = y2 - BH - M
    arrow_down(ax, CX + BW/2, y2, y3 + BH)
    rbox(ax, CX, y3, BW, BH, 'XELOX / FOLFOX-treated mCRC',
         '5 Cohorts Combined (post-QC)',
         fc=C['orange_fc'], ec=C['orange_ec'], lw=1.5, ts=10, ss=7.5)
    sidebar(ax, CX + BW + 0.15, y3 + BH/2, 'n = 476', C['orange_ec'], 9)

    # ── ROW 4: GSE39582 Discovery ──
    y4 = y3 - BH - M
    arrow_down(ax, CX + BW/2, y3, y4 + BH)
    rbox(ax, CX, y4, BW, BH, 'GSE39582 Discovery Cohort',
         'Primary analysis dataset (Affymetrix HG-U133 Plus 2.0)',
         fc=C['blue_fc'], ec=C['blue_ec'], lw=2, ts=11, ss=7.5)
    sidebar(ax, CX + BW + 0.15, y4 + BH/2, 'n = 239  (79 events)', C['blue_ec'], 9)

    # ── ROW 5: Two analysis branches ──
    y5 = y4 - BH - M
    arrow_down(ax, CX + BW/2, y4, y5 + BH)

    # Left: Binary Classification (gene-level)
    w_split = 1.6
    rbox(ax, CX, y5, w_split, BH, 'Binary Classification',
         'Resistant (n=45) vs Sensitive (n=119)',
         fc=C['teal_fc'], ec=C['teal_ec'], lw=1.5, ts=9, ss=6.5)
    sidebar(ax, CX + w_split + 0.1, y5 + BH/2 + 0.12, 'Excluded: n = 75', '#999', 7)
    sidebar(ax, CX + w_split + 0.1, y5 + BH/2 - 0.08, '(intermediate RFS)', '#ccc', 6)

    # Right: Survival Analysis (pathway-level)
    gap = 0.3
    rx = CX + w_split + gap
    rbox(ax, rx, y5, w_split + 0.2, BH, 'Survival Analysis',
         'Relapse-Free Survival, 79 events',
         fc=C['green_fc'], ec=C['green_ec'], lw=1.5, ts=9, ss=6.5)

    # ── ROW 6: Two method branches ──
    y6 = y5 - BH - M
    arrow_down(ax, CX + w_split/2, y5, y6 + BH)
    arrow_down(ax, rx + (w_split + 0.2)/2, y5, y6 + BH)

    # Left: Gene-level ML
    rbox(ax, CX, y6, w_split, BH, 'Gene-Level ML',
         'LASSO-CV → XGBoost + SHAP',
         fc=C['teal_fc'], ec=C['teal_ec'], lw=1.5, ts=9, ss=6.5)

    # Right: Pathway-level Model
    rbox(ax, rx, y6, w_split + 0.2, BH, 'Pathway-Level Model',
         'ssGSEA → LASSO-Cox + AIC',
         fc=C['green_fc'], ec=C['green_ec'], lw=1.5, ts=9, ss=6.5)

    # ── ROW 7: Outputs ──
    y7 = y6 - BH - M
    arrow_down(ax, CX + w_split/2, y6, y7 + BH)
    arrow_down(ax, rx + (w_split + 0.2)/2, y6, y7 + BH)

    # Left output: 10-gene fingerprint
    rbox(ax, CX, y7, w_split, BH, '10-Gene Fingerprint',
         'Apparent AUC = 0.899',
         fc=C['teal_fc'], ec=C['teal_ec'], lw=1, ts=9, ss=6.5)

    # Right output: Nomogram
    rbox(ax, rx, y7, w_split + 0.2, BH, '7-Pathway Nomogram',
         'C-index = 0.676',
         fc=C['green_fc'], ec=C['green_ec'], lw=1, ts=9, ss=6.5)

    # ── ROW 8: External Validation (spans full width) ──
    y8 = y7 - BH - M
    arrow_down(ax, CX + w_split/2, y7, y8 + BH)
    arrow_down(ax, rx + (w_split + 0.2)/2, y7, y8 + BH)

    rbox(ax, CX, y8, rx + w_split + 0.2 - CX, BH, 'External Validation',
         'TCGA-COAD/READ (n=140) | GSE104645 (n=113) | GSE72970 (n=32) | GSE83129 (n=26)',
         fc=C['purple_fc'], ec=C['purple_ec'], lw=1.5, ts=10, ss=7.5)
    sidebar(ax, CX, y8 - 0.2, 'Total: n = 311', C['purple_ec'], 8)

    # ── ROW 9: Results ──
    y9 = y8 - BH - M
    arrow_down(ax, CX + (rx + w_split + 0.2 - CX)/2, y8, y9 + BH)
    rbox(ax, CX, y9, rx + w_split + 0.2 - CX, BH, 'Key Findings',
         'Gene-level: overfitting (AUC 0.54–0.64) | Pathway-level: generalizable (C-index 0.659)',
         fc=C['lime_fc'], ec=C['lime_ec'], lw=1.5, ts=10, ss=7.5)

    # ── RIGHT PANEL: Study Inclusion Summary ──
    panel_x = 4.8
    panel_w = 2.0
    panel_h = 3.5
    panel_y = y3 - 0.1
    box = FancyBboxPatch((panel_x, panel_y), panel_w, panel_h,
                         boxstyle="round,pad=0.03",
                         facecolor='#F5F5F5', edgecolor='#BDBDBD',
                         linewidth=1.2, zorder=1)
    ax.add_patch(box)

    # Dashed arrow from main pipeline to panel
    arrow_right(ax, CX + BW, (y2 + BH/2 + y3 + BH/2)/2,
                panel_x, panel_y + panel_h/2, lw=1, color='#BDBDBD')

    # Panel title
    ax.text(panel_x + panel_w/2, panel_y + panel_h - 0.2, 'Sample Flow Summary',
            ha='center', va='top', fontsize=9, fontweight='bold', color=C['txt'])

    # Panel content
    lines = [
        ('Initial pool', 'n = 1,015'),
        ('After QC', 'n = 476 (47%)'),
        ('GSE39582', 'n = 239'),
        ('  Classification', 'n = 164'),
        ('  Survival', 'n = 227'),
        ('External valid.', 'n = 311'),
    ]
    for i, (label, value) in enumerate(lines):
        ly = panel_y + panel_h - 0.55 - i * 0.42
        ax.text(panel_x + 0.15, ly, label, fontsize=7.5, color=C['mut'], va='center')
        ax.text(panel_x + panel_w - 0.15, ly, value, fontsize=7.5, fontweight='bold',
                color=C['txt'], ha='right', va='center')

    # Separator line
    sep_y = panel_y + panel_h - 0.42
    ax.plot([panel_x + 0.1, panel_x + panel_w - 0.1], [sep_y, sep_y],
            color='#E0E0E0', lw=0.8, zorder=1)

    plt.tight_layout(pad=0.3)
    return fig


if __name__ == '__main__':
    print('Generating Figure 1 (final)...')
    fig = main()

    # Save as PNG
    png_path = os.path.join(OUT_PNG, 'Fig1_pipeline.png')
    fig.savefig(png_path, dpi=DPI, bbox_inches='tight', facecolor='white')
    print(f'  PNG: {png_path} ({os.path.getsize(png_path)//1024} KB)')

    # Save as SVG
    svg_path = os.path.join(OUT_SVG, 'Fig1_pipeline.svg')
    fig.savefig(svg_path, format='svg', bbox_inches='tight', facecolor='white')
    print(f'  SVG: {svg_path} ({os.path.getsize(svg_path)//1024} KB)')

    plt.close(fig)
    print('Done!')
