"""
Figure 4 — Pathway-based Prognostic Risk Model: Construction & Validation
3-panel: A = Forest plot, B = KM curves, C = Risk score density
(Panel D removed — coefficient bar chart was redundant with forest plot)
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import numpy as np
import os
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

# ── Settings ──
DPI = 300
FAMILY = 'Arial'
plt.rcParams.update({
    'font.family': FAMILY,
    'font.size': 8,
    'axes.linewidth': 0.6,
    'figure.dpi': DPI,
    'savefig.dpi': DPI,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.12,
})

# ── Colors ──
C_POS = '#E64B35'
C_NEG = '#4DBBD5'
C_DARK = '#2C3E50'
C_GRAY = '#888888'

BASE = r'/path/to/xelox_project'
OUT = os.path.join(BASE, 'results', 'figures_png', 'v5.11', 'Fig4_model_validation.png')

# ── Variable name mapping (consistent across all figures) ──
VAR_LABELS = {
    'HALLMARK_TGF_BETA_SIGNALING':       'TGF-\u03b2 signaling (HALLMARK)',
    'HALLMARK_WNT_BETA_CATENIN_SIGNALING': 'Wnt/\u03b2-catenin (HALLMARK)',
    'KEGG_ECM_RECEPTOR_INTERACTION':     'ECM receptor interaction (KEGG)',
    'KEGG_TGF_BETA_SIGNALING_PATHWAY':   'TGF-\u03b2 pathway (KEGG)',
    'KEGG_PATHWAYS_IN_CANCER':           'Pathways in cancer (KEGG)',
    'HALLMARK_MYC_TARGETS_V2':           'MYC targets V2 (HALLMARK)',
    'KEGG_COLORECTAL_CANCER':            'Colorectal cancer (KEGG)',
    'LOCATION_DISTAL':                   'Distal tumor location',
}


def fig4():
    cox = pd.read_csv(os.path.join(BASE, 'results', 'tables', 'nomogram', 'cox_final_results.csv'))
    risk = pd.read_csv(os.path.join(BASE, 'results', 'tables', 'nomogram', 'risk_scores.csv'))

    fig = plt.figure(figsize=(6.7, 8.5))
    gs = GridSpec(2, 2, figure=fig, hspace=0.38, wspace=0.30,
                  height_ratios=[1, 1])

    # ═══════════════════════════════════════════
    # A: Forest plot (top-left)
    # ═══════════════════════════════════════════
    ax_a = fig.add_subplot(gs[0, 0])
    cox_sorted = cox.sort_values('HR')
    ypos = list(range(len(cox_sorted)))

    for idx, (i, row) in enumerate(cox_sorted.iterrows()):
        c = C_POS if row['HR'] > 1 else C_NEG
        sig = row['p_value'] < 0.05
        alpha = 1.0 if sig else 0.5
        lw = 2.0 if sig else 1.5
        ax_a.plot([row['HR_lower'], row['HR_upper']], [idx, idx],
                  color=c, linewidth=lw, alpha=alpha)
        ax_a.plot(row['HR'], idx, 'o', color=c, markersize=6,
                  markeredgecolor='white', markeredgewidth=0.5, alpha=alpha)
        # p-value annotation
        pval = row['p_value']
        ptxt = f'p < 0.001' if pval < 0.001 else f'p = {pval:.3f}'
        ax_a.text(3.8, idx, ptxt, ha='right', va='center', fontsize=5.5,
                  color=c if sig else '#AAAAAA')

    ax_a.axvline(1, color=C_DARK, linewidth=0.8, linestyle='--', alpha=0.4)
    ax_a.set_yticks(ypos)
    ax_a.set_yticklabels([VAR_LABELS.get(v, v.replace('_', ' ')) for v in cox_sorted['Variable']],
                         fontsize=7)
    ax_a.set_xlabel('Hazard Ratio (95% CI)', fontsize=8)
    ax_a.set_xlim(0.2, 4.2)
    ax_a.spines['top'].set_visible(False)
    ax_a.spines['right'].set_visible(False)
    ax_a.set_title('Multivariable Cox Regression', fontsize=9, fontweight='bold', color=C_DARK, pad=6)
    ax_a.text(-0.12, 1.03, 'A', transform=ax_a.transAxes,
              fontsize=13, fontweight='bold', va='top', ha='left', color=C_DARK)

    # ═══════════════════════════════════════════
    # B: KM curves (top-right)
    # ═══════════════════════════════════════════
    ax_b = fig.add_subplot(gs[0, 1])
    df = risk.copy()
    df['rfs_delay'] = pd.to_numeric(df['rfs_delay'], errors='coerce')
    df['rfs_event'] = pd.to_numeric(df['rfs_event'], errors='coerce')
    df = df.dropna(subset=['rfs_delay', 'rfs_event'])

    high = df[df['risk_group'] == 'High risk']
    low = df[df['risk_group'] == 'Low risk']

    kmf_h = KaplanMeierFitter()
    kmf_l = KaplanMeierFitter()
    kmf_h.fit(high['rfs_delay'].values, high['rfs_event'].values,
              label=f'High Risk (n={len(high)})')
    kmf_l.fit(low['rfs_delay'].values, low['rfs_event'].values,
              label=f'Low Risk (n={len(low)})')

    kmf_h.plot_survival_function(ax=ax_b, color=C_POS, linewidth=2, ci_show=True, ci_alpha=0.12)
    kmf_l.plot_survival_function(ax=ax_b, color=C_NEG, linewidth=2, ci_show=True, ci_alpha=0.12)

    lr = logrank_test(high['rfs_delay'], low['rfs_delay'],
                      high['rfs_event'], low['rfs_event'])
    p_text = 'log-rank p < 0.001' if lr.p_value < 0.001 else f'log-rank p = {lr.p_value:.4f}'
    ax_b.text(0.55, 0.92, p_text, transform=ax_b.transAxes, fontsize=7.5,
              fontweight='bold',
              bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.9, edgecolor='#CCCCCC'))
    ax_b.set_xlabel('Time (months)', fontsize=8)
    ax_b.set_ylabel('RFS Probability', fontsize=8)
    ax_b.set_ylim(0, 1.05)
    ax_b.legend(fontsize=7, framealpha=0.95, edgecolor='#CCCCCC', loc='lower left')
    ax_b.spines['top'].set_visible(False)
    ax_b.spines['right'].set_visible(False)
    ax_b.set_title('Kaplan-Meier Survival Curves', fontsize=9, fontweight='bold', color=C_DARK, pad=6)
    ax_b.text(-0.12, 1.03, 'B', transform=ax_b.transAxes,
              fontsize=13, fontweight='bold', va='top', ha='left', color=C_DARK)

    # ═══════════════════════════════════════════
    # C: Risk score density (bottom, spanning full width)
    # ═══════════════════════════════════════════
    ax_c = fig.add_subplot(gs[1, :])
    ax_c.hist(low['risk_score'], bins=25, alpha=0.55, color=C_NEG, label='Low Risk', density=True)
    ax_c.hist(high['risk_score'], bins=25, alpha=0.55, color=C_POS, label='High Risk', density=True)
    ax_c.set_xlabel('Risk Score (Linear Predictor)', fontsize=8)
    ax_c.set_ylabel('Density', fontsize=8)
    ax_c.legend(fontsize=7, framealpha=0.95, edgecolor='#CCCCCC')
    ax_c.spines['top'].set_visible(False)
    ax_c.spines['right'].set_visible(False)
    ax_c.set_title('Risk Score Distribution', fontsize=9, fontweight='bold', color=C_DARK, pad=6)
    ax_c.text(-0.06, 1.03, 'C', transform=ax_c.transAxes,
              fontsize=13, fontweight='bold', va='top', ha='left', color=C_DARK)

    # ── Figure caption (bottom) ──
    fig.text(0.5, 0.01,
             'Figure 4. Construction and validation of the 7-pathway prognostic risk model. '
             '(A) Multivariable Cox regression with LASSO-selected pathway features. '
             '(B) Kaplan-Meier RFS curves stratified by risk group. '
             '(C) Risk score density distributions.',
             ha='center', fontsize=6.5, color=C_GRAY, style='italic',
             wrap=True)

    plt.savefig(OUT, dpi=DPI, bbox_inches='tight', pad_inches=0.12)
    print(f'Saved: {OUT}')
    return fig


if __name__ == '__main__':
    fig4()
    print('Done.')
