#!/usr/bin/env python3
"""
Generate Figure 4B: Kaplan-Meier risk stratification curves.
Creates high-quality, standardized KM plot with risk table.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

# Paths
BASE = r'/path/to/xelox_project'
DATA = os.path.join(BASE, 'results', 'tables', 'nomogram', 'risk_scores.csv')
OUT_HD = os.path.join(BASE, 'results', 'figures_png', 'hd', 'Fig4B_prs_km_curve.png')

# Publication settings
DPI = 300
FIG_WIDTH = 6.5
FIG_HEIGHT = 5.5

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': DPI,
    'savefig.dpi': DPI,
    'savefig.bbox': 'tight',
})

print('Generating Figure 4B: KM Risk Stratification...')

# Load risk scores
df = pd.read_csv(DATA)
print(f'  Loaded {len(df)} samples')

# Ensure survival columns exist
if 'rfs_delay' not in df.columns or 'rfs_event' not in df.columns:
    # Try loading clinical data and merging
    clin = pd.read_csv(os.path.join(BASE, 'results', 'tables', 'GSE39582_xelox_groups.csv'))
    df = df.merge(clin[['sample_id', 'rfs_delay', 'rfs_event']], on='sample_id', how='left')
    df = df.dropna(subset=['rfs_delay', 'rfs_event'])
    print(f'  After merge: {len(df)} samples')

# Split by risk group
high_risk = df[df['risk_group'] == 'High risk']
low_risk = df[df['risk_group'] == 'Low risk']

print(f'  High risk: n={len(high_risk)}')
print(f'  Low risk: n={len(low_risk)}')

# Create figure with 2 panels (KM curve + risk table)
fig = plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT))
gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.08)

# Main KM plot
ax_km = fig.add_subplot(gs[0])

# Fit KM curves
kmf_high = KaplanMeierFitter()
kmf_low = KaplanMeierFitter()

kmf_high.fit(high_risk['rfs_delay'], high_risk['rfs_event'], label='High Risk')
kmf_low.fit(low_risk['rfs_delay'], low_risk['rfs_event'], label='Low Risk')

# Plot
kmf_low.plot_survival_function(ax=ax_km, color='#4DBBD5', linewidth=2.5, ci_show=True, ci_alpha=0.15)
kmf_high.plot_survival_function(ax=ax_km, color='#E64B35', linewidth=2.5, ci_show=True, ci_alpha=0.15)

# Log-rank test
results = logrank_test(high_risk['rfs_delay'], low_risk['rfs_delay'],
                       event_observed_A=high_risk['rfs_event'],
                       event_observed_B=low_risk['rfs_event'])

# Add p-value annotation
p_val = results.p_value
if p_val < 0.001:
    p_str = 'p < 0.001'
else:
    p_str = f'p = {p_val:.4f}'

ax_km.text(0.62, 0.95, f'Log-rank {p_str}',
           transform=ax_km.transAxes, fontsize=11, fontweight='bold',
           verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))

ax_km.set_xlabel('')
ax_km.set_ylabel('RFS Probability', fontsize=11)
ax_km.set_title('Kaplan-Meier Risk Stratification\nGSE39582 XELOX Cohort (n = 227)',
                fontsize=12, fontweight='bold', pad=10)
ax_km.set_ylim(0, 1.05)
ax_km.set_xlim(0, 170)
ax_km.legend(loc='lower left', framealpha=0.9)
ax_km.grid(True, alpha=0.3, linestyle='--')
ax_km.set_xticklabels([])  # Hide x-labels for risk table alignment

# Risk table
ax_risk = fig.add_subplot(gs[1], sharex=ax_km)

time_points = [0, 12, 24, 36, 48, 60, 72, 84, 96, 120, 144, 168]

# Count at risk
n_low = []
n_high = []
for t in time_points:
    n_low.append(sum(low_risk['rfs_delay'] >= t))
    n_high.append(sum(high_risk['rfs_delay'] >= t))

ax_risk.bar([t - 2 for t in time_points], n_low, width=4, color='#4DBBD5', alpha=0.7, label='Low Risk')
ax_risk.bar([t + 2 for t in time_points], n_high, width=4, color='#E64B35', alpha=0.7, label='High Risk')

ax_risk.set_xlabel('Time (months)', fontsize=11)
ax_risk.set_ylabel('No. at Risk', fontsize=10)
ax_risk.set_xlim(0, 170)
ax_risk.set_ylim(0, max(max(n_low), max(n_high)) * 1.15)
ax_risk.set_xticks(time_points)
ax_risk.legend(loc='upper right', fontsize=8, ncol=2)
ax_risk.grid(True, axis='y', alpha=0.3, linestyle='--')

plt.tight_layout()

# Save
os.makedirs(os.path.dirname(OUT_HD), exist_ok=True)
fig.savefig(OUT_HD, dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close(fig)

print(f'  Saved: {OUT_HD}')
print('Done!')
