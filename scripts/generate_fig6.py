"""
Generate Figure 6B: Feature Type x Modeling Method Comparison (C-index / AUC)

Output: results/figures_png/Fig6B_method_comparison.png
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# ==========================================================
# Fig6A: Forest Plot of 23-gene AIC-selected Cox Model (unchanged)
# ==========================================================

df = pd.read_csv(
    r"/path/to/xelox_project\results\tables\gene_level\gene_cox_aic_model.csv"
)
df = df.sort_values('HR', ascending=True).reset_index(drop=True)
log_hr = np.log(df['HR'].values)
log_lower = np.log(df['HR_lower'].values)
log_upper = np.log(df['HR_upper'].values)
p_colors = []
for p in df['p_value']:
    if p < 0.001: p_colors.append('#2166ac')
    elif p < 0.01: p_colors.append('#4393c3')
    elif p < 0.05: p_colors.append('#92c5de')
    else: p_colors.append('#d1e5f0')

fig, ax = plt.subplots(figsize=(8, 7))
y_pos = np.arange(len(df))
ax.axvline(x=0, color='gray', linestyle='--', linewidth=0.8, alpha=0.6)
for i in range(len(df)):
    ax.errorbar(log_hr[i], y_pos[i],
                xerr=[[log_hr[i] - log_lower[i]], [log_upper[i] - log_hr[i]]],
                fmt='o', color=p_colors[i], ecolor=p_colors[i], capsize=3, capthick=1.5,
                markersize=7, markeredgecolor='white', markeredgewidth=0.5, zorder=5)
ax.set_yticks(y_pos)
ax.set_yticklabels(df['Gene'].str.replace('`', '').values, fontsize=9)
xticks_log = np.linspace(-1.2, 0.9, 8)
ax.set_xticks(xticks_log)
ax.set_xticklabels([f'{x:.2f}' for x in np.exp(xticks_log)], fontsize=9)
ax.set_xlabel('Hazard Ratio (95% CI)', fontsize=11, fontweight='bold')
ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim())
ax2.set_xticks([])
ax2.text(0.5, 1.1, 'Protective -> Risk', transform=ax.transAxes,
         ha='center', fontsize=10, style='italic', color='gray')
ax.axvspan(0, ax.get_xlim()[0], alpha=0.04, color='#2166ac', zorder=-1)
ax.axvspan(0, ax.get_xlim()[1], alpha=0.04, color='#b2182b', zorder=-1)
for i in range(len(df)):
    hr_text = f'{df["HR"].iloc[i]:.2f} ({df["HR_lower"].iloc[i]:.2f}-{df["HR_upper"].iloc[i]:.2f})'
    ax.text(log_upper[i] + 0.02, y_pos[i], hr_text, fontsize=7.5, va='center')
ax.set_ylabel('Gene', fontsize=11, fontweight='bold')
ax.set_title('A  Gene-Level Ablation Study: 23-Gene AIC-Selected Cox Model',
             fontsize=12, fontweight='bold', loc='left', pad=15)
ax.tick_params(axis='y', length=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.text(0.01, -0.12,
        'EPV = 3.4 (79 events / 23 variables) | C-index = 0.834 | Severe overfitting risk',
        transform=ax.transAxes, fontsize=9, style='italic', color='#b2182b')
plt.tight_layout()
fig.savefig(
    r"/path/to/xelox_project\results\figures_png\Fig6A_ablation_forest.png",
    dpi=300, bbox_inches='tight'
)
plt.close()
print("[OK] Fig6A saved")

# ==========================================================
# Fig6B: Method Comparison Bar Chart - NO overlapping text
# Key changes:
#   1) Insight text moved OUTSIDE the plot area (below the figure, as a caption block)
#   2) Larger canvas with more whitespace
#   3) All annotations carefully positioned to never touch bars
# ==========================================================

fig, ax = plt.subplots(figsize=(12, 9))

models = [
    ('Gene-level\nXGBoost', 'XGBoost', 0.888, '#d1e5f0', '-', 'Cross-platform\nfailure'),
    ('Pathway-level\nLASSO-PRS', 'LASSO', 0.100, '#fddbc7', '-', 'Jaccard=0.10\n49.2% lambda.min fallback'),
    ('Gene-level\nAIC Cox', 'Backward AIC', 0.834, '#f4a582', '3.4', 'Severe\noverfitting'),
    ('Pathway-level\nAIC Nomogram', 'Backward AIC', 0.676, '#2166ac', '9.9', 'Bootstrap-corrected\nC-index=0.692'),
]

n = len(models)
x_pos = np.arange(n)
bar_width = 0.50
colors = [m[3] for m in models]

bars = ax.bar(x_pos, [m[2] for m in models], bar_width,
              color=colors, edgecolor='gray', linewidth=1.3, alpha=0.88)

# --- Value labels on top of each bar ---
for i, bar in enumerate(bars):
    val = models[i][2]
    epv = models[i][4]
    # Main value above bar
    ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.02,
            f'{val:.3f}', ha='center', va='bottom', fontsize=13, fontweight='bold')
    # EPV label inside bar for tall bars, below short bar
    if epv != '-':
        if val > 0.5:   # tall bar: put EPV text inside the bar body
            y_epv = bar.get_height() * 0.55
            ax.text(bar.get_x() + bar.get_width()/2., y_epv,
                    f'EPV={epv}', ha='center', va='center',
                    fontsize=10, color='white', fontweight='bold')
        else:           # short bar: put EPV text just above the value
            ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.06,
                    f'EPV={epv}', ha='center', va='bottom',
                    fontsize=9, color='#333333', fontweight='bold')

# --- X-axis labels ---
ax.set_xticks(x_pos)
ax.set_xticklabels([m[0] for m in models], fontsize=11, fontweight='bold')

# --- Note under each x-tick ---
for i, m in enumerate(models):
    ax.text(x_pos[i], -0.14, m[5],
            ha='center', va='top', fontsize=9, color='#555555',
            style='italic', transform=ax.get_xaxis_transform())

# --- Y-axis ---
ax.set_ylabel('Discrimination Metric', fontsize=12, fontweight='bold')
ax.set_ylim(0, 1.15)

# Reference lines
ax.axhline(y=0.7, color='gray', linestyle=':', linewidth=0.8, alpha=0.6)
ax.axhline(y=0.5, color='gray', linestyle=':', linewidth=0.8, alpha=0.6)
ax.text(n - 0.35, 0.703, '0.7 (clinically useful)', fontsize=9, color='gray', alpha=0.7, va='bottom')
ax.text(n - 0.35, 0.497, '0.5 (random)', fontsize=9, color='gray', alpha=0.7, va='top')

# --- Top grouping brackets (method names) ---
brackets = [
    (-0.25, 0.25, 'XGBoost', '#4393c3'),
    (0.75, 1.25, 'LASSO', '#d6604d'),
    (1.75, 2.25, 'Backward AIC', '#2166ac'),
    (2.75, 3.25, 'Backward AIC', '#2166ac'),
]
for start, end, label, clr in brackets:
    mid = (start + end) / 2.
    ax.annotate('', xy=(start, 1.04), xytext=(end, 1.04),
                arrowprops=dict(arrowstyle='-', color='gray', lw=0.8))
    ax.text(mid, 1.08, label, ha='center', va='bottom', fontsize=9, color=clr, fontweight='bold')

# Feature-type grouping
ax.text(0.5, 1.19, '<- Gene-level ->', ha='center', va='bottom', fontsize=10, fontweight='bold', color='#333')
ax.text(2.5, 1.19, '<- Pathway-level ->', ha='center', va='bottom', fontsize=10, fontweight='bold', color='#333')

# Style
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_title('B  Feature Type x Modeling Method Comparison',
             fontsize=13, fontweight='bold', loc='left', pad=28)

# --- INSIGHT BOX: placed at bottom-left, well clear of all bars ---
props = dict(boxstyle='round,pad=0.45', facecolor='#f5f5f5', edgecolor='#bbbbbb', alpha=0.95)
insight_lines = [
    'Methodological Insight:',
    '* Gene-level + backward AIC: highest apparent C-index (0.834)',
    '  but critically low EPV (3.4) => severe overfitting',
    '* Pathway-level + backward AIC: lower but honest C-index (0.676)',
    '  with adequate EPV (9.9), clinically useful',
    '* Pathway aggregation = implicit regularization that improves',
    '  reproducibility and constrains model space',
]
insight_text = '\n'.join(insight_lines)
ax.text(0.01, -0.38, insight_text, transform=ax.transAxes, fontsize=8.8,
        verticalalignment='top', horizontalalignment='left', bbox=props,
        family='monospace')

plt.tight_layout()
plt.subplots_adjust(bottom=0.30)  # extra room for the insight box
fig.savefig(
    r"/path/to/xelox_project\results\figures_png\Fig6B_method_comparison.png",
    dpi=300, bbox_inches='tight'
)
plt.close()
print("[OK] Fig6B saved")
print("Done! Both Figure 6 panels generated.")
