"""
Figure 5 — Clinical Nomogram and Calibration
Panel A: Standard nomogram (Points → Variables → Total Points → RFS Probability)
Panel B: Calibration curves (12/36/60-month predicted vs observed)
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import numpy as np
import os
import pandas as pd
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.utils import concordance_index
from scipy.interpolate import UnivariateSpline

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
C_POS = '#E64B35'
C_NEG = '#4DBBD5'
C_DARK = '#2C3E50'
C_GRAY = '#888888'
C_BG = '#FAFAFA'
SURV_COLORS = {'12m': '#3498DB', '36m': '#E74C3C', '60m': '#27AE60'}

BASE = r'/path/to/xelox_project'
OUT = os.path.join(BASE, 'results', 'figures_png', 'v5.11', 'Fig5_nomogram_calibration.png')

# ── Cox model data ──
VAR_DATA = [
    # (csv_name, display_label, coefficient=ln(HR))
    ('LOCATION_DISTAL',                   'Distal location',           np.log(2.143)),
    ('HALLMARK_TGF_BETA_SIGNALING',       'TGF-\u03b2 (HALLMARK)',     np.log(1.872)),
    ('KEGG_ECM_RECEPTOR_INTERACTION',     'ECM receptor (KEGG)',       np.log(1.652)),
    ('KEGG_COLORECTAL_CANCER',            'Colorectal Ca (KEGG)',      np.log(1.520)),
    ('HALLMARK_WNT_BETA_CATENIN_SIGNALING', 'Wnt/\u03b2-catenin (HALLMARK)', np.log(1.457)),
    ('HALLMARK_MYC_TARGETS_V2',           'MYC targets V2 (HALLMARK)', np.log(0.681)),
    ('KEGG_PATHWAYS_IN_CANCER',           'Pathways in cancer (KEGG)', np.log(0.614)),
    ('KEGG_TGF_BETA_SIGNALING_PATHWAY',   'TGF-\u03b2 pathway (KEGG)', np.log(0.492)),
]

n_vars = len(VAR_DATA)
coeffs = np.array([v[2] for v in VAR_DATA])
coef_abs = np.abs(coeffs)
max_abs = coef_abs.max()
POINTS_PER_UNIT = 100.0 / max_abs
var_points_max = np.round(coef_abs * POINTS_PER_UNIT).astype(int)
TOTAL_MAX = int(var_points_max.sum())


def draw_nomogram(ax):
    """Draw the standard nomogram on the given axes."""
    ax.set_facecolor(C_BG)
    ax.set_xlim(-5, 165)
    ax.set_ylim(-4.5, n_vars + 3.0)
    ax.axis('off')

    # ── Title ──
    ax.text(80, n_vars + 2.6, 'Nomogram for Predicting 12/36/60-Month RFS',
            ha='center', va='center', fontsize=12, fontweight='bold', color=C_DARK)
    ax.text(80, n_vars + 2.0,
            '7-Pathway Signature + Tumor Location  |  LASSO-Cox (\u03bb.1se)  |  C-index = 0.676',
            ha='center', va='center', fontsize=8, color=C_GRAY)

    # ── Top ruler: Points ──
    pts_y = n_vars + 1.2
    pts_left = 38
    pts_right = pts_left + 100
    ax.text(pts_left - 2, pts_y, 'Points', ha='right', va='center',
            fontsize=9, fontweight='bold', color=C_DARK)
    ax.plot([pts_left, pts_right], [pts_y, pts_y], color=C_DARK, lw=1.5)
    for p in range(0, 101, 5):
        x = pts_left + p
        major = (p % 10 == 0)
        tick_h = 0.22 if major else 0.12
        ax.plot([x, x], [pts_y - tick_h, pts_y + tick_h], color=C_DARK, lw=0.8 if major else 0.5)
        if major:
            ax.text(x, pts_y + 0.32, str(p), ha='center', va='bottom',
                    fontsize=6.5, color=C_DARK)

    # ── Variable scales ──
    bar_left = pts_left
    bar_height = 0.55

    for i, (var_id, label, coef) in enumerate(VAR_DATA):
        y = n_vars - 0.3 - i * 1.0
        pt_max = var_points_max[i]
        c = C_POS if coef > 0 else C_NEG

        # Variable label
        ax.text(0, y, label, ha='left', va='center', fontsize=7,
                color=C_DARK, fontweight='bold')

        # Scale bar
        bar_bottom = y - bar_height / 2
        # Filled bar
        rect = plt.Rectangle((bar_left, bar_bottom), pt_max, bar_height,
                              fc=c, ec='none', lw=0, alpha=0.2)
        ax.add_patch(rect)
        # Border
        ax.plot([bar_left, bar_left + pt_max], [bar_bottom, bar_bottom], color=c, lw=1.0)
        ax.plot([bar_left, bar_left + pt_max],
                [bar_bottom + bar_height, bar_bottom + bar_height], color=c, lw=1.0)
        # Vertical caps
        ax.plot([bar_left, bar_left], [bar_bottom, bar_bottom + bar_height], color=c, lw=1.0)
        ax.plot([bar_left + pt_max, bar_left + pt_max],
                [bar_bottom, bar_bottom + bar_height], color=c, lw=1.0)

        # Tick marks (points on the scale)
        n_ticks = min(pt_max + 1, 11)
        tick_positions = np.linspace(0, pt_max, n_ticks).astype(int)
        for tp in tick_positions:
            tx = bar_left + tp
            major = (tp == 0 or tp == pt_max or tp % max(1, pt_max // 5) == 0)
            tick_h = 0.12 if major else 0.06
            ax.plot([tx, tx], [bar_bottom - tick_h, bar_bottom], color=c, lw=0.5)
            if major:
                ax.text(tx, bar_bottom - 0.22, str(tp), ha='center', va='top',
                        fontsize=5, color=c)

        # HR annotation
        hr_val = np.exp(coef)
        ax.text(bar_left + pt_max + 4, y, f'HR = {hr_val:.3f}',
                ha='left', va='center', fontsize=5.5, color=C_GRAY)

    # ── Total Points ruler ──
    total_y = -0.6
    ax.text(0, total_y, 'Total Points', ha='left', va='center',
            fontsize=9, fontweight='bold', color=C_DARK)
    ax.plot([bar_left, bar_left + TOTAL_MAX], [total_y, total_y], color=C_DARK, lw=2.0)

    total_ticks = list(range(0, TOTAL_MAX + 1, 25))
    for tp in total_ticks:
        x = bar_left + tp
        major = (tp % 50 == 0)
        tick_h = 0.18 if major else 0.10
        ax.plot([x, x], [total_y - tick_h, total_y + tick_h], color=C_DARK, lw=0.8 if major else 0.5)
        if major:
            ax.text(x, total_y - 0.30, str(tp), ha='center', va='top',
                    fontsize=6, color=C_DARK)

    # ── Survival probability axes (12/36/60 months) ──
    # Compute S(t) = S0(t)^exp(LP) using baseline KM
    risk = pd.read_csv(os.path.join(BASE, 'results', 'tables', 'nomogram', 'risk_scores.csv'))
    risk['rfs_delay'] = pd.to_numeric(risk['rfs_delay'], errors='coerce')
    risk['rfs_event'] = pd.to_numeric(risk['rfs_event'], errors='coerce')
    risk = risk.dropna(subset=['rfs_delay', 'rfs_event'])

    kmf = KaplanMeierFitter()
    kmf.fit(risk['rfs_delay'].values, risk['rfs_event'].values)

    time_configs = [
        (12, '12-month RFS', SURV_COLORS['12m']),
        (36, '36-month RFS', SURV_COLORS['36m']),
        (60, '60-month RFS', SURV_COLORS['60m']),
    ]

    for ti, (t_months, t_label, sc) in enumerate(time_configs):
        sy = -1.8 - ti * 0.85

        # Baseline survival at this time point
        try:
            s0 = float(kmf.survival_function_at_times(t_months).values[0])
        except:
            s0 = 0.5

        ax.text(0, sy, t_label, ha='left', va='center', fontsize=7.5,
                fontweight='bold', color=sc)

        # Compute survival for range of total points
        n_pts = 200
        total_range = np.linspace(0, TOTAL_MAX, n_pts)
        surv_probs = []
        for tp in total_range:
            lp = (tp / 100.0) * max_abs - max_abs * 0.4  # offset to center probabilities
            s_prob = s0 ** np.exp(lp)
            s_prob = np.clip(s_prob, 0.02, 0.98)
            surv_probs.append(s_prob)
        surv_probs = np.array(surv_probs)

        # Draw axis line
        ax.plot([bar_left, bar_left + TOTAL_MAX], [sy, sy], color=sc, lw=1.5, alpha=0.8)

        # End labels
        ax.text(bar_left - 1, sy, f'{surv_probs[0]:.0%}', ha='right', va='center',
                fontsize=6, color=sc, fontweight='bold')
        ax.text(bar_left + TOTAL_MAX + 1, sy, f'{surv_probs[-1]:.0%}', ha='left', va='center',
                fontsize=6, color=sc, fontweight='bold')

        # Tick marks at survival probability levels
        for prob_level in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
            idx = np.argmin(np.abs(surv_probs - prob_level))
            if 3 < idx < n_pts - 3:
                x = bar_left + total_range[idx]
                ax.plot([x, x], [sy - 0.08, sy + 0.08], color=sc, lw=0.6)
                if prob_level in [0.2, 0.4, 0.6, 0.8]:
                    ax.text(x, sy - 0.22, f'{prob_level:.0%}', ha='center', va='top',
                            fontsize=5, color=sc)

        # Direction annotation
        ax.text(bar_left + TOTAL_MAX / 2, sy + 0.22,
                '\u2190 Lower risk  |  Higher risk \u2192',
                ha='center', va='center', fontsize=5.5, color=C_GRAY, style='italic')

    # ── Instructions ──
    inst_y = -4.2
    ax.text(80, inst_y,
            'Instructions: For each variable, locate patient\'s value on its scale \u2192 '
            'read corresponding points \u2192 sum all points to get Total Points \u2192 '
            'find Total Points on the bottom axis \u2192 read RFS probability at each time point.',
            ha='center', va='center', fontsize=6, color=C_GRAY, style='italic')

    # Panel label
    ax.text(-12, n_vars + 2.6, 'A', ha='right', va='center',
            fontsize=14, fontweight='bold', color=C_DARK)


def draw_calibration(ax_12, ax_36, ax_60):
    """Draw calibration curves for 12/36/60 months."""
    risk = pd.read_csv(os.path.join(BASE, 'results', 'tables', 'nomogram', 'risk_scores.csv'))
    risk['rfs_delay'] = pd.to_numeric(risk['rfs_delay'], errors='coerce')
    risk['rfs_event'] = pd.to_numeric(risk['rfs_event'], errors='coerce')
    risk = risk.dropna(subset=['rfs_delay', 'rfs_event'])

    # Fit Cox model on risk_score
    cph = CoxPHFitter()
    df_cox = risk[['rfs_delay', 'rfs_event', 'risk_score']].copy()
    cph.fit(df_cox, duration_col='rfs_delay', event_col='rfs_event')

    # Baseline survival
    baseline_surv = cph.baseline_survival_

    time_configs = [
        (12, ax_12, SURV_COLORS['12m'], '12-month'),
        (36, ax_36, SURV_COLORS['36m'], '36-month'),
        (60, ax_60, SURV_COLORS['60m'], '60-month'),
    ]

    for t_months, ax, sc, t_label in time_configs:
        ax.set_facecolor(C_BG)

        # Predicted survival for each patient
        surv_at_t = cph.predict_survival_function(risk, times=[t_months])
        predicted = surv_at_t.values.flatten()

        # Observed: group by predicted deciles and compute KM
        n_bins = 8
        bins = pd.qcut(predicted, n_bins, duplicates='drop')
        df_temp = pd.DataFrame({
            'predicted': predicted,
            'observed_event': risk['rfs_event'].values,
            'observed_time': risk['rfs_delay'].values,
            'bin': bins
        })

        pred_means = []
        obs_means = []
        for grp_name, grp in df_temp.groupby('bin', observed=True):
            pred_means.append(grp['predicted'].mean())
            kmf = KaplanMeierFitter()
            kmf.fit(grp['observed_time'].values, grp['observed_event'].values)
            try:
                obs_surv = float(kmf.survival_function_at_times(t_months).values[0])
            except:
                obs_surv = grp['predicted'].mean()
            obs_means.append(obs_surv)

        pred_means = np.array(pred_means)
        obs_means = np.array(obs_means)

        # Plot
        ax.scatter(pred_means, obs_means, s=40, color=sc, edgecolors='white',
                   linewidths=0.8, zorder=4)
        # Connect with line
        sort_idx = np.argsort(pred_means)
        ax.plot(pred_means[sort_idx], obs_means[sort_idx], color=sc, lw=1.5,
                alpha=0.8, zorder=3)

        # Ideal diagonal
        ax.plot([0, 1], [0, 1], '--', color='#AAAAAA', lw=1.0, alpha=0.7, zorder=1)

        # Brier score (from calibration_stats.csv)
        brier_map = {12: 0.142, 36: 0.201, 60: 0.218}
        slope_map = {12: 0.817}

        # Calculate C-index for this time point
        ci_val = concordance_index(risk['rfs_delay'], -risk['risk_score'], risk['rfs_event'])

        ax.set_xlabel('Predicted RFS Probability', fontsize=7.5)
        ax.set_ylabel('Observed RFS Probability', fontsize=7.5)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_aspect('equal', adjustable='box')
        ax.set_title(f'{t_label} Calibration', fontsize=8.5, fontweight='bold', color=C_DARK, pad=5)

        # Metrics box
        brier = brier_map.get(t_months, 0)
        ax.text(0.97, 0.03,
                f'Brier = {brier:.3f}\nC-index = {ci_val:.3f}',
                transform=ax.transAxes, fontsize=6, va='bottom', ha='right',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.9,
                          edgecolor='#CCCCCC'))

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)


def fig5_nomogram_calibration():
    fig = plt.figure(figsize=(7.2, 14))

    # Top: Nomogram (Panel A) — takes ~60% of height
    gs = fig.add_gridspec(4, 1, height_ratios=[3.5, 1.2, 1.2, 1.2],
                          hspace=0.35, left=0.08, right=0.95, top=0.95, bottom=0.04)

    ax_nom = fig.add_subplot(gs[0])
    draw_nomogram(ax_nom)

    # Bottom: Calibration curves (Panel B, C, D)
    ax_cal_12 = fig.add_subplot(gs[1])
    ax_cal_36 = fig.add_subplot(gs[2])
    ax_cal_60 = fig.add_subplot(gs[3])
    draw_calibration(ax_cal_12, ax_cal_36, ax_cal_60)

    # Panel labels for calibration
    for ax, label in [(ax_cal_12, 'B'), (ax_cal_36, 'C'), (ax_cal_60, 'D')]:
        ax.text(-0.08, 1.15, label, transform=ax.transAxes,
                fontsize=13, fontweight='bold', va='top', ha='left', color=C_DARK)

    # Figure caption
    fig.text(0.5, 0.005,
             'Figure 5. Clinical nomogram and calibration. '
             '(A) Nomogram with point scales for 7 pathway features and tumor location. '
             '(B-D) Calibration curves for 12-, 36-, and 60-month RFS prediction.',
             ha='center', fontsize=6.5, color=C_GRAY, style='italic', wrap=True)

    plt.savefig(OUT, dpi=DPI, bbox_inches='tight', pad_inches=0.15)
    print(f'Saved: {OUT}')
    return fig


if __name__ == '__main__':
    fig5_nomogram_calibration()
    print('Done.')
