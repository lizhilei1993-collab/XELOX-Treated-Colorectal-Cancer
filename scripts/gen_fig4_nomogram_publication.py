"""
Figure 4 — Publication-Quality Clinical Nomogram (Python)
Mimics rms::nomogram() output from R.
Standard layout: Points → Variables → Total Points → LP → RFS Probabilities
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np
import os
import pandas as pd
from lifelines import KaplanMeierFitter

# ── Style ──
DPI = 300
FAMILY = 'Arial'
plt.rcParams.update({
    'font.family': FAMILY,
    'font.size': 8,
    'axes.linewidth': 0.6,
    'figure.dpi': DPI,
    'savefig.dpi': DPI,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.2,
})

C_RISK = '#C0392B'      # risk-increasing (dark red)
C_PROT = '#2980B9'       # protective (dark blue)
C_DARK = '#2C3E50'
C_GRAY = '#7F8C8D'
C_LIGHTGRAY = '#BDC3C7'
C_BG = '#FAFAFA'

BASE = r'/path/to/xelox_project'
OUT  = os.path.join(BASE, 'results', 'figures_png', 'v5.11', 'Fig4_nomogram_publication.png')

# ── Cox model data (from cox_final_results.csv) ──
# Ordered by |coefficient| descending for visual hierarchy
VARIABLES = [
    # (label, ln(HR), HR, p-value, range_type, tick_values)
    ('Tumor location\n(distal)',          np.log(2.143),  2.143,  0.007,  'binary', [0, 1]),
    ('TGF-β signaling\n(HALLMARK)',       np.log(1.872),  1.872,  0.041,  'cont',   [0.0, 0.3, 0.5, 0.7, 0.9, 1.0]),
    ('ECM receptor interaction\n(KEGG)',   np.log(1.652),  1.652,  0.014,  'cont',   [0.0, 0.3, 0.5, 0.7, 0.9, 1.0]),
    ('Colorectal cancer\n(KEGG)',          np.log(1.520),  1.520,  0.013,  'cont',   [0.0, 0.3, 0.5, 0.7, 0.9, 1.0]),
    ('Wnt/β-catenin\n(HALLMARK)',         np.log(1.457),  1.457,  0.001,  'cont',   [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]),
    ('MYC targets V2\n(HALLMARK)',        np.log(0.681),  0.681,  0.013,  'cont',   [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]),
    ('Pathways in cancer\n(KEGG)',        np.log(0.614),  0.614,  0.070,  'cont',   [0.0, 0.3, 0.5, 0.7, 0.9, 1.0]),
    ('TGF-β pathway\n(KEGG)',             np.log(0.492),  0.492,  0.007,  'binary', [0, 1]),
]

N_VARS = len(VARIABLES)
COEFS = np.array([v[1] for v in VARIABLES])
COEF_ABS = np.abs(COEFS)
MAX_ABS = COEF_ABS.max()
POINTS_PER_UNIT = 100.0 / MAX_ABS
VAR_POINTS_MAX = np.round(COEF_ABS * POINTS_PER_UNIT).astype(int)
TOTAL_MAX = int(VAR_POINTS_MAX.sum())


def draw_ruler(ax, y, x_start, x_end, n_major, n_minor_per_major,
               labels, color=C_DARK, lw=1.2, tick_up=0.25, tick_down=0.25,
               label_above=True, fontsize=7, label_fmt='{}'):
    """Draw a horizontal ruler with tick marks and labels."""
    ax.plot([x_start, x_end], [y, y], color=color, lw=lw, zorder=3)
    total_range = x_end - x_start
    major_step = total_range / n_major
    minor_step = major_step / n_minor_per_major

    for i in range(n_major * n_minor_per_major + 1):
        x = x_start + i * minor_step
        is_major = (i % n_minor_per_major == 0)
        if is_major:
            tu, td = tick_up, tick_down
            mlw = lw * 0.9
        else:
            tu, td = tick_up * 0.5, tick_down * 0.5
            mlw = lw * 0.5
        ax.plot([x, x], [y - td, y + tu], color=color, lw=mlw, zorder=3)

    # Major labels
    if labels is not None:
        for i, lbl in enumerate(labels):
            x = x_start + i * major_step
            if isinstance(lbl, float):
                txt = f'{lbl:.1f}'
            else:
                txt = str(lbl)
            if label_above:
                ax.text(x, y + tick_up + 0.08, txt, ha='center', va='bottom',
                        fontsize=fontsize, color=color)
            else:
                ax.text(x, y - tick_down - 0.08, txt, ha='center', va='top',
                        fontsize=fontsize, color=color)


def fig4_nomogram():
    fig = plt.figure(figsize=(7.5, 11.5))
    ax = fig.add_axes([0.05, 0.04, 0.90, 0.92])  # [left, bottom, width, height]
    ax.set_facecolor(C_BG)

    # Layout constants
    LABEL_WIDTH = 22    # left margin for variable names
    SCALE_LEFT = 24     # where scales start
    SCALE_WIDTH = 55    # width of points scale (0–100)
    SCALE_RIGHT = SCALE_LEFT + SCALE_WIDTH

    ax.set_xlim(-2, 95)
    ax.set_ylim(-7.5, N_VARS + 5.5)
    ax.axis('off')

    # ═══════════════════════════════════════════
    # Title
    # ═══════════════════════════════════════════
    title_y = N_VARS + 4.8
    ax.text(45, title_y, 'Nomogram for Predicting 12/36/60-Month Recurrence-Free Survival',
            ha='center', va='center', fontsize=13, fontweight='bold', color=C_DARK)
    ax.text(45, title_y - 0.8,
            'LASSO-Cox (λ.1se)  ·  7-Pathway Signature + Tumor Location  ·  C-index = 0.676',
            ha='center', va='center', fontsize=8.5, color=C_GRAY)
    ax.text(45, title_y - 1.4,
            'GSE39582 XELOX cohort (n = 227, 79 events)',
            ha='center', va='center', fontsize=7.5, color=C_LIGHTGRAY, style='italic')

    # ═══════════════════════════════════════════
    # Top ruler: Points (0–100)
    # ═══════════════════════════════════════════
    pts_y = N_VARS + 2.8
    ax.text(SCALE_LEFT - 1.5, pts_y, 'Points', ha='right', va='center',
            fontsize=9.5, fontweight='bold', color=C_DARK)

    # Draw the points ruler
    pts_labels = list(range(0, 101, 10))
    draw_ruler(ax, pts_y, SCALE_LEFT, SCALE_RIGHT, 10, 2,
               labels=pts_labels, color=C_DARK, lw=1.5,
               tick_up=0.30, tick_down=0.30, label_above=False, fontsize=6.5)

    # ═══════════════════════════════════════════
    # Variable scales
    # ═══════════════════════════════════════════
    for i, (label, coef, hr, pval, rtype, ticks) in enumerate(VARIABLES):
        y = N_VARS - 0.5 - i * 1.0
        pt_max = VAR_POINTS_MAX[i]
        color = C_RISK if coef > 0 else C_PROT

        # Variable name (left side)
        ax.text(LABEL_WIDTH, y, label, ha='right', va='center', fontsize=7.5,
                color=C_DARK, fontweight='bold', linespacing=0.9)

        # Scale bar: maps [0, pt_max] points
        # The bar starts at SCALE_LEFT and extends pt_max units
        bar_left = SCALE_LEFT
        bar_right = SCALE_LEFT + pt_max
        bar_h = 0.32

        # Filled background
        rect = plt.Rectangle((bar_left, y - bar_h), pt_max, bar_h * 2,
                              fc=color, ec='none', alpha=0.12, zorder=1)
        ax.add_patch(rect)

        # Top and bottom borders
        ax.plot([bar_left, bar_right], [y + bar_h, y + bar_h], color=color, lw=0.8, zorder=2)
        ax.plot([bar_left, bar_right], [y - bar_h, y - bar_h], color=color, lw=0.8, zorder=2)
        # End caps
        ax.plot([bar_left, bar_left], [y - bar_h, y + bar_h], color=color, lw=0.8, zorder=2)
        ax.plot([bar_right, bar_right], [y - bar_h, y + bar_h], color=color, lw=0.8, zorder=2)

        # Tick marks with values
        n_ticks = len(ticks)
        for j, val in enumerate(ticks):
            if rtype == 'binary':
                tx = bar_left + val * pt_max
            else:
                tx = bar_left + (val / max(ticks)) * pt_max

            # Vertical tick
            ax.plot([tx, tx], [y - bar_h, y - bar_h - 0.15], color=color, lw=0.5, zorder=2)

            # Label
            if rtype == 'binary':
                lbl_txt = 'Yes' if val == 1 else 'No'
            else:
                lbl_txt = f'{val:.1f}' if val != int(val) else f'{int(val)}'
            ax.text(tx, y - bar_h - 0.25, lbl_txt, ha='center', va='top',
                    fontsize=5.5, color=color)

        # Point scale above the bar (showing 0, pt_max/2, pt_max)
        pt_ticks = np.linspace(0, pt_max, min(pt_max + 1, 6)).astype(int)
        for pt in pt_ticks:
            tx = bar_left + pt
            ax.plot([tx, tx], [y + bar_h, y + bar_h + 0.10], color=color, lw=0.4, zorder=2)
            if pt == 0 or pt == pt_max or pt == pt_max // 2:
                ax.text(tx, y + bar_h + 0.15, str(pt), ha='center', va='bottom',
                        fontsize=4.5, color=color, alpha=0.8)

    # ═══════════════════════════════════════════
    # Total Points ruler
    # ═══════════════════════════════════════════
    total_y = -1.2
    ax.text(LABEL_WIDTH, total_y, 'Total\nPoints', ha='right', va='center',
            fontsize=9.5, fontweight='bold', color=C_DARK, linespacing=0.85)

    total_labels = list(range(0, TOTAL_MAX + 1, 50))
    n_major = TOTAL_MAX // 50
    draw_ruler(ax, total_y, SCALE_LEFT, SCALE_LEFT + TOTAL_MAX, n_major, 2,
               labels=total_labels, color=C_DARK, lw=1.8,
               tick_up=0.28, tick_down=0.28, label_above=False, fontsize=6)

    # ═══════════════════════════════════════════
    # Linear Predictor ruler
    # ═══════════════════════════════════════════
    lp_y = -2.5
    ax.text(LABEL_WIDTH, lp_y, 'Linear\nPredictor', ha='right', va='center',
            fontsize=8, fontweight='bold', color=C_GRAY, linespacing=0.85)

    lp_left = SCALE_LEFT
    lp_right = SCALE_LEFT + TOTAL_MAX

    # LP values: map total points to LP
    # LP = sum(coef_i * x_i), for display we show range
    lp_labels = []
    for tp in range(0, TOTAL_MAX + 1, 50):
        lp_val = (tp / 100.0) * MAX_ABS - MAX_ABS * 0.5
        lp_labels.append(f'{lp_val:.1f}')

    draw_ruler(ax, lp_y, lp_left, lp_right, n_major, 2,
               labels=lp_labels, color=C_GRAY, lw=1.0,
               tick_up=0.18, tick_down=0.18, label_above=False, fontsize=5.5)

    # ═══════════════════════════════════════════
    # Survival Probability axes (12/36/60 months)
    # ═══════════════════════════════════════════
    risk = pd.read_csv(os.path.join(BASE, 'results', 'tables', 'nomogram', 'risk_scores.csv'))
    risk['rfs_delay'] = pd.to_numeric(risk['rfs_delay'], errors='coerce')
    risk['rfs_event'] = pd.to_numeric(risk['rfs_event'], errors='coerce')
    risk = risk.dropna(subset=['rfs_delay', 'rfs_event'])

    kmf = KaplanMeierFitter()
    kmf.fit(risk['rfs_delay'].values, risk['rfs_event'].values)

    surv_configs = [
        (12, '12-Month\nRFS Probability', '#3498DB'),
        (36, '36-Month\nRFS Probability', '#E74C3C'),
        (60, '60-Month\nRFS Probability', '#27AE60'),
    ]

    for ti, (t_months, t_label, sc) in enumerate(surv_configs):
        sy = -3.8 - ti * 1.2

        # Label
        ax.text(LABEL_WIDTH, sy, t_label, ha='right', va='center',
                fontsize=7.5, fontweight='bold', color=sc, linespacing=0.85)

        # Baseline survival
        try:
            s0 = float(kmf.survival_function_at_times(t_months).values[0])
        except:
            s0 = 0.5

        # Compute survival probability across total points range
        n_pts = 300
        total_range = np.linspace(0, TOTAL_MAX, n_pts)
        surv_probs = []
        for tp in total_range:
            lp = (tp / 100.0) * MAX_ABS - MAX_ABS * 0.5
            s_prob = s0 ** np.exp(lp)
            s_prob = np.clip(s_prob, 0.02, 0.98)
            surv_probs.append(s_prob)
        surv_probs = np.array(surv_probs)

        # Draw the axis line (non-linear, use piecewise linear)
        axis_left = SCALE_LEFT
        axis_right = SCALE_LEFT + TOTAL_MAX
        ax.plot([axis_left, axis_right], [sy, sy], color=sc, lw=1.5, alpha=0.6, zorder=1)

        # End labels
        ax.text(axis_left - 0.8, sy, f'{surv_probs[0]:.0%}', ha='right', va='center',
                fontsize=7, color=sc, fontweight='bold')
        ax.text(axis_right + 0.8, sy, f'{surv_probs[-1]:.0%}', ha='left', va='center',
                fontsize=7, color=sc, fontweight='bold')

        # Tick marks at specific survival probabilities
        target_probs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        for prob in target_probs:
            idx = np.argmin(np.abs(surv_probs - prob))
            if 5 < idx < n_pts - 5:
                x = axis_left + total_range[idx]
                major = prob in [0.2, 0.4, 0.6, 0.8]
                tu = 0.18 if major else 0.10
                td = 0.18 if major else 0.10
                ax.plot([x, x], [sy - td, sy + tu], color=sc, lw=0.6 if major else 0.4)
                if major:
                    ax.text(x, sy - 0.30, f'{prob:.0%}', ha='center', va='top',
                            fontsize=5.5, color=sc)

        # Direction arrow annotation
        ax.text(axis_left + TOTAL_MAX / 2, sy + 0.35,
                '← Lower risk  ·  Higher risk →',
                ha='center', va='center', fontsize=5.5, color=C_LIGHTGRAY, style='italic')

    # ═══════════════════════════════════════════
    # Divider lines between sections
    # ═══════════════════════════════════════════
    for divider_y in [N_VARS + 1.8, -0.5, -1.8, -3.1]:
        ax.plot([SCALE_LEFT - 2, SCALE_LEFT + TOTAL_MAX + 2], [divider_y, divider_y],
                color=C_LIGHTGRAY, lw=0.3, alpha=0.5, zorder=0)

    # ═══════════════════════════════════════════
    # Instructions (bottom)
    # ═══════════════════════════════════════════
    inst_y = -6.8
    ax.text(45, inst_y,
            'Instructions: For each variable, locate the patient\'s value on its scale → '
            'draw a vertical line to read points → sum all points on the Total Points axis → '
            'project down to read 12/36/60-month RFS probability.',
            ha='center', va='center', fontsize=6.5, color=C_GRAY, style='italic',
            wrap=True)

    ax.text(45, inst_y - 0.8,
            'Significant predictors (p < 0.05) carry more points. '
            'Red = risk-increasing, Blue = protective.  |  B = 200 bootstrap for calibration.',
            ha='center', va='center', fontsize=6, color=C_LIGHTGRAY, style='italic')

    # ═══════════════════════════════════════════
    # p-value annotations (right side)
    # ═══════════════════════════════════════════
    for i, (label, coef, hr, pval, rtype, ticks) in enumerate(VARIABLES):
        y = N_VARS - 0.5 - i * 1.0
        bar_right = SCALE_LEFT + VAR_POINTS_MAX[i]
        color = C_RISK if coef > 0 else C_PROT
        sig = '*' if pval < 0.05 else ''
        ptxt = f'p={pval:.3f}{sig}'
        ax.text(bar_right + 1.5, y, ptxt, ha='left', va='center',
                fontsize=5, color=color, alpha=0.7)

    plt.savefig(OUT, dpi=DPI, bbox_inches='tight', pad_inches=0.2)
    print(f'Saved: {OUT}')
    return fig


if __name__ == '__main__':
    fig4_nomogram()
    print('Done.')
