#!/usr/bin/env python3
"""
Generate Figure 4C: Forest plot of the 8-variable AIC-selected Cox regression model.
Based on cox_final_results.csv (Script 29 output).
Variables: 7 pathways + location_distal
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Paths
BASE = r'/path/to/xelox_project'
CSV_PATH = os.path.join(BASE, 'results', 'tables', 'nomogram', 'cox_final_results.csv')
OUT_PNG = os.path.join(BASE, 'results', 'figures_png', 'Fig4C_cox_forest.png')
OUT_HD = os.path.join(BASE, 'results', 'figures_png', 'hd', 'Fig4C_cox_forest.png')

# Load data
df = pd.read_csv(CSV_PATH)
print("Loaded Cox final results:")
print(df)

# Variable display names
name_map = {
    'HALLMARK_TGF_BETA_SIGNALING': 'TGF-β signaling (Hallmark)',
    'HALLMARK_WNT_BETA_CATENIN_SIGNALING': 'Wnt/β-catenin signaling (Hallmark)',
    'KEGG_ECM_RECEPTOR_INTERACTION': 'ECM-receptor interaction (KEGG)',
    'KEGG_TGF_BETA_SIGNALING_PATHWAY': 'TGF-β signaling pathway (KEGG)',
    'KEGG_PATHWAYS_IN_CANCER': 'Pathways in cancer (KEGG)',
    'HALLMARK_MYC_TARGETS_V2': 'MYC targets v2 (Hallmark)',
    'KEGG_COLORECTAL_CANCER': 'Colorectal cancer (KEGG)',
    'location_distal': 'Tumor location (distal vs. proximal)',
}

df['Display'] = df['Variable'].map(name_map)
df['Significant'] = df['p_value'] < 0.05

# Sort by HR descending
df = df.sort_values('HR', ascending=True).reset_index(drop=True)

# Create figure
fig, ax = plt.subplots(figsize=(8, 5), dpi=300)

y_positions = range(len(df))
colors = ['#C0392B' if sig else '#7F8C8D' for sig in df['Significant']]

# Plot CIs
for i, row in df.iterrows():
    y = i
    hr = row['HR']
    lower = row['HR_lower']
    upper = row['HR_upper']
    col = colors[i]
    # Horizontal line for CI
    ax.plot([lower, upper], [y, y], color=col, linewidth=2, solid_capstyle='round')
    # Point for HR
    ax.plot(hr, y, 'o', color=col, markersize=7, markeredgecolor='white', markeredgewidth=0.8, zorder=3)

# Reference line at HR=1
ax.axvline(x=1, color='#2C3E50', linewidth=0.8, linestyle='--', alpha=0.7)

# Y-axis labels
ax.set_yticks(list(y_positions))
ax.set_yticklabels(df['Display'], fontsize=9)

# X-axis
ax.set_xlabel('Hazard Ratio (95% CI)', fontsize=10)
ax.set_xlim(0.3, 3.0)
ax.set_xticks([0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0])
ax.set_xticklabels(['0.4', '0.6', '0.8', '1.0', '1.2', '1.5', '2.0', '2.5', '3.0'], fontsize=8)

# Add HR text annotations
for i, row in df.iterrows():
    y = i
    hr = row['HR']
    lower = row['HR_lower']
    upper = row['HR_upper']
    p = row['p_value']
    sig_mark = '*' if p < 0.05 else ''
    text = f"{hr:.2f} ({lower:.2f}–{upper:.2f}){sig_mark}"
    ax.text(3.05, y, text, va='center', ha='left', fontsize=8, color='#2C3E50')

# Title
ax.set_title(
    'Multivariable Cox Regression: 8-Variable Nomogram Model\n'
    'GSE39582 XELOX cohort (n = 227, 79 RFS events, EPV = 9.9)',
    fontsize=11, fontweight='bold', pad=10
)

# Legend
red_patch = mpatches.Patch(color='#C0392B', label='p < 0.05')
gray_patch = mpatches.Patch(color='#7F8C8D', label='p ≥ 0.05')
ax.legend(handles=[red_patch, gray_patch], loc='lower right', fontsize=8, framealpha=0.9)

# Spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(0.6)
ax.spines['bottom'].set_linewidth(0.6)
ax.tick_params(axis='both', length=3)

# Tight layout
plt.tight_layout()

# Save
os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
os.makedirs(os.path.dirname(OUT_HD), exist_ok=True)
fig.savefig(OUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(OUT_HD, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Saved: {OUT_PNG}")
print(f"Saved HD: {OUT_HD}")

plt.close(fig)
