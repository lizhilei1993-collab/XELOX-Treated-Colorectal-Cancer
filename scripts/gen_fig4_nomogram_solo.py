"""
Figure 4D — Standalone Nomogram (Publication Quality)
7-Pathway + Location Cox model for predicting RFS.
Designed for BMC Cancer / Briefings in Bioinformatics.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os
import pandas as pd
from lifelines import KaplanMeierFitter

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
    'savefig.pad_inches': 0.15,
})

# ── Colors ──
C_POS = '#E64B35'   # risk-increasing (red)
C_NEG = '#4DBBD5'   # protective (teal)
C_DARK = '#2C3E50'
C_GRAY = '#888888'
C_BG = '#FAFAFA'

BASE = r'/path/to/xelox_project'
OUT = os.path.join(BASE, 'results', 'figures_png', 'v5.11', 'Fig4D_nomogram_standalone.png')

# ── Data (from cox_final_results.csv, ln(HR) = coefficient) ──
variables = [
    ('LOCATION_DISTAL',        'Distal tumor\nlocation',       np.log(2.143),  'p = 0.007'),
    ('HALLMARK_TGF_BETA',      'TGF-\u03b2 signaling\n(HALLMARK)',   np.log(1.872),  'p = 0.041'),
    ('KEGG_ECM_RECEPTOR',      'ECM receptor\ninteraction',    np.log(1.652),  'p = 0.014'),
    ('KEGG_COLORECTAL',        'Colorectal cancer\n(KEGG)',    np.log(1.520),  'p = 0.013'),
    ('HALLMARK_WNT',           'Wnt/\u03b2-catenin\n(HALLMARK)', np.log(1.457),  'p = 0.001'),
    ('HALLMARK_MYC',           'MYC targets V2\n(HALLMARK)',   np.log(0.681),  'p = 0.013'),
    ('KEGG_PATHWAYS_CANCER',   'Pathways in\ncancer',          np.log(0.614),  'p = 0.070'),
    ('KEGG_TGF_BETA',          'TGF-\u03b2 pathway\n(KEGG)',   np.log(0.492),  'p = 0.007'),
]

n_vars = len(variables)
coeffs = np.array([v[2] for v in variables])
coef_abs = np.abs(coeffs)
max_abs = coef_abs.max()

# ── Point system: max coefficient = 100 points ──
POINTS_PER_UNIT = 100.0 / max_abs
var_points = np.round(coef_abs * POINTS_PER_UNIT).astype(int)
total_max = var_points.sum()  # ≈ 430

# ── Representative variable values for tick marks ──
# (continuous variables: quartiles; binary: 0/1)
var_ranges = [
    (0, 1),         # Distal (binary)
    (0.3, 0.5, 0.7, 0.9),  # TGF-β HALLMARK (continuous)
    (0.3, 0.5, 0.7, 0.9),  # ECM receptor (continuous)
    (0.3, 0.5, 0.7, 0.9),  # Colorectal cancer (continuous)
    (0.3, 0.5, 0.7, 0.9),  # Wnt/β-catenin (continuous)
    (0.3, 0.5, 0.7, 0.9),  # MYC targets (continuous)
    (0.3, 0.5, 0.7, 0.9),  # Pathways in cancer (continuous)
    (0, 1),         # TGF-β KEGG (binary in model context)
]


def fig4_nomogram_solo():
    fig, ax = plt.subplots(figsize=(7.2, 10.5))
    ax.set_facecolor(C_BG)
    ax.set_xlim(0, 160)
    ax.set_ylim(-2.5, n_vars + 3.2)
    ax.axis('off')

    # ═══════════════════════════════════════════
    # Title block
    # ═══════════════════════════════════════════
    ax.text(80, n_vars + 2.8, 'Nomogram for Predicting Recurrence-Free Survival',
            ha='center', va='center', fontsize=13, fontweight='bold', color=C_DARK)
    ax.text(80, n_vars + 2.2, '7-Pathway Signature + Tumor Location  |  LASSO-Cox (λ.1se)',
            ha='center', va='center', fontsize=8.5, color=C_GRAY)
    ax.text(80, n_vars + 1.7, 'GSE39582 XELOX cohort (n = 227, 79 events)  ·  C-index = 0.676',
            ha='center', va='center', fontsize=7.5, color=C_GRAY, style='italic')

    # ═══════════════════════════════════════════
    # "Points" scale at top
    # ═══════════════════════════════════════════
    pts_y = n_vars + 0.9
    pts_left = 40
    pts_right = pts_left + 100
    ax.text(pts_left - 2, pts_y + 0.35, 'Points', ha='right', va='center',
            fontsize=8, fontweight='bold', color=C_DARK)
    ax.plot([pts_left, pts_right], [pts_y, pts_y], color=C_DARK, lw=1.2)
    for p in range(0, 101, 10):
        x = pts_left + p
        tick_h = 0.18 if p % 20 == 0 else 0.1
        ax.plot([x, x], [pts_y - tick_h, pts_y + tick_h], color=C_DARK, lw=0.8)
        if p % 20 == 0:
            ax.text(x, pts_y - 0.28, str(p), ha='center', va='top',
                    fontsize=6.5, color=C_DARK)

    # ═══════════════════════════════════════════
    # Variable scales
    # ═══════════════════════════════════════════
    bar_left = pts_left
    bar_height = 0.5
    y_spacing = 1.0

    for i, (var_id, label, coef, pval) in enumerate(variables):
        y = n_vars - 0.2 - i * y_spacing
        pt_max = var_points[i]
        coef_val = coeffs[i]
        c = C_POS if coef_val > 0 else C_NEG

        # Variable name (left)
        ax.text(2, y, label, ha='left', va='center', fontsize=7,
                color=C_DARK, fontweight='bold')

        # Scale bar background
        bar_bottom = y - bar_height / 2
        rect = plt.Rectangle((bar_left, bar_bottom), pt_max, bar_height,
                              fc=c, ec='white', lw=0.5, alpha=0.25)
        ax.add_patch(rect)

        # Scale bar border
        ax.plot([bar_left, bar_left + pt_max], [bar_bottom, bar_bottom], color=c, lw=1.0)
        ax.plot([bar_left, bar_left + pt_max], [bar_bottom + bar_height, bar_bottom + bar_height],
                color=c, lw=1.0)

        # Tick marks on the scale
        rng = var_ranges[i]
        n_ticks = len(rng)
        for j, val in enumerate(rng):
            # Map value to points: 0 value = 0 points, max value = pt_max points
            tx = bar_left + (val / max(rng)) * pt_max
            ax.plot([tx, tx], [bar_bottom - 0.08, bar_bottom], color=c, lw=0.6)

        # For binary: label 0 and 1
        if len(rng) == 2:
            ax.text(bar_left, bar_bottom - 0.18, 'No', ha='center', va='top',
                    fontsize=5.5, color=C_GRAY)
            ax.text(bar_left + pt_max, bar_bottom - 0.18, 'Yes', ha='center', va='top',
                    fontsize=5.5, color=C_GRAY)
        else:
            # Continuous: show quartile labels
            for j, val in enumerate(rng):
                tx = bar_left + (val / max(rng)) * pt_max
                ax.text(tx, bar_bottom - 0.18, f'{val:.1f}', ha='center', va='top',
                        fontsize=5, color=C_GRAY)

        # HR label at end of bar
        hr_val = np.exp(coef_val)
        ax.text(bar_left + pt_max + 3, y, f'HR = {hr_val:.3f}',
                ha='left', va='center', fontsize=5.5, color=C_GRAY)

        # Direction arrow (use matplotlib arrow to avoid font glyph issues)
        if coef_val > 0:
            ax.annotate('', xy=(bar_left + pt_max + 2, y),
                        xytext=(bar_left + pt_max - 0.5, y),
                        arrowprops=dict(arrowstyle='->', color=c, lw=1.2, mutation_scale=8))
        else:
            ax.annotate('', xy=(bar_left - 2, y),
                        xytext=(bar_left + 0.5, y),
                        arrowprops=dict(arrowstyle='->', color=c, lw=1.2, mutation_scale=8))

    # ═══════════════════════════════════════════
    # Total Points line
    # ═══════════════════════════════════════════
    total_y = -0.3
    ax.text(2, total_y + 0.25, 'Total Points', ha='left', va='center',
            fontsize=8.5, fontweight='bold', color=C_DARK)

    ax.plot([bar_left, bar_left + total_max], [total_y, total_y], color=C_DARK, lw=1.5)
    n_total = total_max // 50
    for p in range(0, total_max + 1, 50):
        x = bar_left + p
        tick_h = 0.12 if p % 100 == 0 else 0.08
        ax.plot([x, x], [total_y - tick_h, total_y + tick_h], color=C_DARK, lw=0.8)
        ax.text(x, total_y - 0.25, str(p), ha='center', va='top',
                fontsize=6, color=C_DARK)

    # ═══════════════════════════════════════════
    # Linear Predictor
    # ═══════════════════════════════════════════
    lp_y = -0.85
    ax.text(2, lp_y + 0.25, 'Linear\nPredictor', ha='left', va='center',
            fontsize=7, fontweight='bold', color=C_DARK)

    # Map total points to LP: LP = (total / 100) * max_abs - mean_lp
    # For display: range from min to max of possible LPs
    lp_min = -(total_max / 100) * max_abs / 2
    lp_max = (total_max / 100) * max_abs / 2

    ax.plot([bar_left, bar_left + total_max], [lp_y, lp_y], color='#95A5A6', lw=1.0)
    for p in range(0, total_max + 1, 50):
        x = bar_left + p
        lp_val = (p / 100.0) * max_abs - max_abs / 2
        ax.plot([x, x], [lp_y - 0.06, lp_y + 0.06], color='#95A5A6', lw=0.6)
        if p % 100 == 0:
            ax.text(x, lp_y - 0.22, f'{lp_val:.2f}', ha='center', va='top',
                    fontsize=5.5, color=C_GRAY)

    # ═══════════════════════════════════════════
    # Survival Probability Axes (12/36/60 months)
    # ═══════════════════════════════════════════
    surv_configs = [
        (12, '#3498DB', '12-month RFS'),
        (36, '#E74C3C', '36-month RFS'),
        (60, '#27AE60', '60-month RFS'),
    ]

    # Calculate baseline survival from KM
    risk = pd.read_csv(os.path.join(BASE, 'results', 'tables', 'nomogram', 'risk_scores.csv'))
    risk['rfs_delay'] = pd.to_numeric(risk['rfs_delay'], errors='coerce')
    risk['rfs_event'] = pd.to_numeric(risk['rfs_event'], errors='coerce')
    risk = risk.dropna(subset=['rfs_delay', 'rfs_event'])

    kmf = KaplanMeierFitter()
    kmf.fit(risk['rfs_delay'].values, risk['rfs_event'].values)

    for ti, (t_months, sc, label) in enumerate(surv_configs):
        sy = -1.55 - ti * 0.5

        # Get baseline survival
        try:
            s0 = kmf.survival_function_at_times(t_months).values[0]
        except:
            s0 = 0.5

        ax.text(2, sy, label, ha='left', va='center', fontsize=7,
                fontweight='bold', color=sc)

        # Survival axis: for each total point value, compute survival probability
        # S(t) = S0(t)^exp(LP), where LP is linear predictor
        n_pts = 100
        total_vals = np.linspace(0, total_max, n_pts)
        surv_probs = []
        for tp in total_vals:
            lp = (tp / 100.0) * max_abs - max_abs / 2
            s_prob = s0 ** np.exp(lp)
            s_prob = np.clip(s_prob, 0.01, 0.99)
            surv_probs.append(s_prob)

        surv_probs = np.array(surv_probs)

        # Draw the survival axis
        ax.plot([bar_left, bar_left + total_max], [sy, sy], color=sc, lw=1.0, alpha=0.8)

        # Label ends
        ax.text(bar_left - 1, sy, f'{surv_probs[0]:.0%}', ha='right', va='center',
                fontsize=5.5, color=sc, fontweight='bold')
        ax.text(bar_left + total_max + 1, sy, f'{surv_probs[-1]:.0%}', ha='left', va='center',
                fontsize=5.5, color=sc, fontweight='bold')

        # Tick marks at 10% intervals
        for prob_level in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
            # Find x position where survival = prob_level
            idx = np.argmin(np.abs(surv_probs - prob_level))
            if 5 < idx < n_pts - 5:
                x = bar_left + total_vals[idx]
                ax.plot([x, x], [sy - 0.06, sy + 0.06], color=sc, lw=0.5)
                if prob_level in [0.2, 0.4, 0.6, 0.8]:
                    ax.text(x, sy - 0.17, f'{prob_level:.0%}', ha='center', va='top',
                            fontsize=4.5, color=sc)

        # "High Risk" and "Low Risk" labels
        ax.text(bar_left + total_max / 2, sy + 0.18, f'\u2190 Lower risk   |   Higher risk \u2192',
                ha='center', va='center', fontsize=5, color=C_GRAY, style='italic')

    # ═══════════════════════════════════════════
    # Usage instructions
    # ═══════════════════════════════════════════
    inst_y = -3.2
    ax.text(80, inst_y,
            'Instructions: Locate each variable\'s value \u2192 draw a vertical line to the Points scale \u2192 '
            'Sum all points \u2192 Find the total on the Total Points axis \u2192 '
            'Read the corresponding survival probabilities below.',
            ha='center', va='center', fontsize=6.5, color=C_GRAY, style='italic',
            wrap=True)

    # ═══════════════════════════════════════════
    # Panel label
    # ═══════════════════════════════════════════
    ax.text(-5, n_vars + 2.8, 'D', ha='right', va='center',
            fontsize=14, fontweight='bold', color=C_DARK)

    plt.savefig(OUT, dpi=DPI, bbox_inches='tight', pad_inches=0.15)
    print(f'Saved: {OUT}')
    return fig


if __name__ == '__main__':
    fig4_nomogram_solo()
    print('Done.')
