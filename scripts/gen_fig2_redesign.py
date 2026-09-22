"""
Figure 2 — SHAP Analysis & Ensemble Voting Framework (Publication Quality)
Panel A: SHAP summary (top XGBoost features by mean |SHAP|)
Panel B: Algorithm voting matrix for the final 10-gene fingerprint
Exports: 600 DPI TIFF, SVG, PNG
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import numpy as np
import os
import pandas as pd

# ── Style (Nature Communications / Scientific Reports) ──
DPI = 600
FAMILY = 'Arial'
plt.rcParams.update({
    'font.family': FAMILY,
    'font.size': 8,
    'axes.linewidth': 0.6,
    'figure.dpi': DPI,
    'savefig.dpi': DPI,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})

# ── Colors ──
C_DARK = '#2C3E50'
C_GRAY = '#7F8C8D'
C_LIGHT = '#BDC3C7'
C_FINGERPRINT = '#C0392B'     # red — in final fingerprint
C_SHAP_ONLY = '#3498DB'       # blue — SHAP top but not in fingerprint
C_BG = '#FAFAFA'

BASE = r'/path/to/xelox_project'
OUT_DIR = os.path.join(BASE, 'results', 'figures_png', 'v5.11')

# ── Algorithm display names and colors ──
ALGO_MAP = {
    'LASSO':           ('LASSO',          '#E74C3C'),
    'Ridge':           ('Ridge',          '#E67E22'),
    'ElasticNet':      ('ElasticNet',     '#F1C40F'),
    'RandomForest':    ('Random Forest',  '#2ECC71'),
    'XGBoost_Gain':    ('XGB Gain',       '#1ABC9C'),
    'XGBoost_SHAP':    ('XGB SHAP',       '#3498DB'),
    'Boruta_Confirmed':('Boruta',         '#9B59B6'),
    'Boruta_Top15':    ('Boruta Top15',   '#8E44AD'),
    'LogReg_AIC':      ('LogReg AIC',     '#2980B9'),
    'PLS_DA':          ('PLS-DA',         '#16A085'),
    'CART':            ('CART',           '#D35400'),
}

ALGO_ORDER = ['LASSO', 'Ridge', 'ElasticNet', 'RandomForest',
              'XGBoost_Gain', 'XGBoost_SHAP', 'Boruta_Confirmed',
              'Boruta_Top15', 'LogReg_AIC', 'PLS_DA', 'CART']


def fig2():
    # ── Load data ──
    shap = pd.read_csv(os.path.join(BASE, 'results', 'tables', 'ml_phase2', 'SHAP_summary_phase2.csv'))
    fp = pd.read_csv(os.path.join(BASE, 'results', 'tables', 'ml_phase2', 'XELOX_resistance_fingerprint.csv'))

    # Fingerprint gene set
    fp_genes = set(fp['Gene'].tolist())

    # Clean gene labels
    def clean_gene(g):
        return g.replace(' /// ', '/')

    # ── Figure layout ──
    fig = plt.figure(figsize=(7.2, 6.5))
    gs = GridSpec(2, 1, figure=fig, hspace=0.40, height_ratios=[1, 1.3],
                  left=0.12, right=0.95, top=0.92, bottom=0.06)

    # ═══════════════════════════════════════════
    # Panel A: SHAP summary bar chart
    # ═══════════════════════════════════════════
    ax_a = fig.add_subplot(gs[0])
    ax_a.set_facecolor(C_BG)

    shap_sorted = shap.sort_values('MeanAbsSHAP', ascending=True).reset_index(drop=True)
    n_shap = len(shap_sorted)

    bar_colors = []
    for g in shap_sorted['Gene']:
        if g in fp_genes:
            bar_colors.append(C_FINGERPRINT)
        else:
            bar_colors.append(C_SHAP_ONLY)

    y_pos = np.arange(n_shap)
    bars = ax_a.barh(y_pos, shap_sorted['MeanAbsSHAP'], height=0.65,
                     color=bar_colors, edgecolor='white', linewidth=0.5, zorder=3)

    # Gene labels
    for i, row in shap_sorted.iterrows():
        label = clean_gene(row['Gene'])
        in_fp = row['Gene'] in fp_genes
        weight = 'bold' if in_fp else 'normal'
        color = C_DARK if in_fp else C_GRAY
        ax_a.text(-0.02, i, label, ha='right', va='center', fontsize=7.5,
                  fontweight=weight, color=color)

    # Value labels
    for i, (v, g) in enumerate(zip(shap_sorted['MeanAbsSHAP'], shap_sorted['Gene'])):
        ax_a.text(v + 0.01, i, f'{v:.3f}', va='center', ha='left',
                  fontsize=6.5, color=C_GRAY)

    # Fingerprint badge
    for i, g in enumerate(shap_sorted['Gene']):
        if g in fp_genes:
            ax_a.text(shap_sorted['MeanAbsSHAP'].max() * 0.95, i, '\u2605',
                      ha='right', va='center', fontsize=9, color=C_FINGERPRINT)

    ax_a.set_yticks([])
    ax_a.set_xlabel('Mean |SHAP| Value (XGBoost)', fontsize=9)
    ax_a.set_xlim(0, shap_sorted['MeanAbsSHAP'].max() * 1.25)
    ax_a.invert_yaxis()
    ax_a.spines['top'].set_visible(False)
    ax_a.spines['right'].set_visible(False)
    ax_a.spines['left'].set_visible(False)
    ax_a.set_title('XGBoost Feature Importance (SHAP)', fontsize=10, fontweight='bold',
                   color=C_DARK, pad=8)

    # Legend for Panel A
    legend_a = [mpatches.Patch(color=C_FINGERPRINT, label='In 10-gene fingerprint'),
                mpatches.Patch(color=C_SHAP_ONLY, label='Top SHAP (not in fingerprint)')]
    ax_a.legend(handles=legend_a, fontsize=6.5, loc='lower right',
                framealpha=0.95, edgecolor='#CCCCCC')
    ax_a.text(-0.15, 1.05, 'A', transform=ax_a.transAxes,
              fontsize=14, fontweight='bold', va='top', ha='left', color=C_DARK)

    # ═══════════════════════════════════════════
    # Panel B: Algorithm voting matrix (dot plot)
    # ═══════════════════════════════════════════
    ax_b = fig.add_subplot(gs[1])
    ax_b.set_facecolor(C_BG)

    # Build voting matrix
    fp_data = fp.sort_values('AlgorithmVotes', ascending=False).reset_index(drop=True)
    n_genes = len(fp_data)
    n_algos = len(ALGO_ORDER)

    # Gene labels (y-axis)
    gene_labels = [clean_gene(g) for g in fp_data['Gene']]
    vote_counts = fp_data['AlgorithmVotes'].values

    # Draw dot matrix
    dot_size = 120
    for gi, row in fp_data.iterrows():
        algos_voted = [a.strip() for a in row['SupportingAlgorithms'].split(';')]
        for ai, algo in enumerate(ALGO_ORDER):
            if algo in algos_voted:
                algo_label, algo_color = ALGO_MAP[algo]
                ax_b.scatter(ai, gi, s=dot_size, color=algo_color,
                             edgecolors='white', linewidths=0.5, zorder=4, marker='o')
            else:
                # Subtle background dot
                ax_b.scatter(ai, gi, s=20, color='#EAECEE',
                             edgecolors='none', zorder=1, marker='o')

    # Gene labels (left) with vote count
    for gi, (gl, vc) in enumerate(zip(gene_labels, vote_counts)):
        ax_b.text(-0.5, gi, gl, ha='right', va='center', fontsize=7.5,
                  fontweight='bold', color=C_DARK)
        ax_b.text(n_algos + 0.3, gi, f'{vc}/11', ha='left', va='center',
                  fontsize=7, color=C_FINGERPRINT, fontweight='bold')

    # Algorithm labels (top)
    for ai, algo in enumerate(ALGO_ORDER):
        algo_label, algo_color = ALGO_MAP[algo]
        ax_b.text(ai, -0.8, algo_label, ha='center', va='bottom',
                  fontsize=6, color=algo_color, fontweight='bold', rotation=45)

    # Workflow arrow annotation
    ax_b.annotate('Top SHAP \u2192 Ensemble Voting \u2192 10-Gene Fingerprint',
                  xy=(n_algos / 2, n_genes + 0.3), fontsize=8, fontweight='bold',
                  ha='center', va='bottom', color=C_DARK,
                  bbox=dict(boxstyle='round,pad=0.4', facecolor='#EBF5FB',
                            edgecolor='#3498DB', linewidth=0.8))

    ax_b.set_xlim(-6, n_algos + 2.5)
    ax_b.set_ylim(-1.5, n_genes + 1.0)
    ax_b.set_xticks([])
    ax_b.set_yticks([])
    ax_b.spines['top'].set_visible(False)
    ax_b.spines['right'].set_visible(False)
    ax_b.spines['bottom'].set_visible(False)
    ax_b.spines['left'].set_visible(False)
    ax_b.set_title('Ensemble Voting Framework: Algorithm Consensus',
                   fontsize=10, fontweight='bold', color=C_DARK, pad=8)

    # Algorithm color legend (bottom)
    algo_handles = []
    for algo in ALGO_ORDER:
        _, ac = ALGO_MAP[algo]
        al, _ = ALGO_MAP[algo]
        algo_handles.append(mpatches.Patch(color=ac, label=al))
    ax_b.legend(handles=algo_handles, fontsize=5.5, loc='lower center',
                ncol=4, framealpha=0.95, edgecolor='#CCCCCC',
                bbox_to_anchor=(0.5, -0.18))

    ax_b.text(-0.15, 1.05, 'B', transform=ax_b.transAxes,
              fontsize=14, fontweight='bold', va='top', ha='left', color=C_DARK)

    # ── Cohort info ──
    fig.text(0.5, 0.005,
             'GSE39582 XELOX cohort (n = 227)  ·  XGBoost AUC = 0.934  ·  '
             '10 algorithms, consensus threshold \u2265 5 votes',
             ha='center', fontsize=6, color=C_LIGHT, style='italic')

    # ── Export ──
    # PNG
    png_path = os.path.join(OUT_DIR, 'Fig2_SHAP_voting.png')
    plt.savefig(png_path, dpi=DPI, bbox_inches='tight', pad_inches=0.1)
    print(f'PNG saved: {png_path}')

    # TIFF (600 DPI)
    tiff_path = os.path.join(OUT_DIR, 'Fig2_SHAP_voting.tif')
    plt.savefig(tiff_path, dpi=600, format='tiff', bbox_inches='tight', pad_inches=0.1)
    print(f'TIFF saved: {tiff_path}')

    # SVG (vector)
    svg_path = os.path.join(OUT_DIR, 'Fig2_SHAP_voting.svg')
    plt.savefig(svg_path, format='svg', bbox_inches='tight', pad_inches=0.1)
    print(f'SVG saved: {svg_path}')

    return fig


if __name__ == '__main__':
    fig2()
    print('Done.')
