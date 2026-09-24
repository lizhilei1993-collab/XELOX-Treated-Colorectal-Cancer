#!/usr/bin/env python3
"""
Figure 2: SHAP Analysis from XGBoost model + training data CSV.
Uses XGBoost's built-in predcontrib for SHAP values.
No external SHAP library needed.
"""
import os, sys, warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xgboost as xgb

# Paths
BASE = r'/path/to/xelox_project'
MODEL_PATH = os.path.join(BASE, 'results/tables/ml_phase2/XGBoost_model_phase2.xgb')
TRAIN_CSV = r'/tmp/scratch\training_shap_data.csv'
OUT_DIR = os.path.join(BASE, 'results/figures_png/v5.11')
os.makedirs(OUT_DIR, exist_ok=True)

# Settings
DPI = 300
FAMILY = 'Arial'
W1 = 3.35   # single column 85mm

plt.rcParams.update({
    'font.family': FAMILY, 'font.size': 9,
    'axes.labelsize': 9, 'axes.titlesize': 0,
    'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'legend.fontsize': 7, 'figure.dpi': DPI,
    'savefig.dpi': DPI, 'savefig.facecolor': 'white',
    'figure.facecolor': 'white', 'savefig.bbox': 'tight',
})

print('=== Figure 2: SHAP Analysis ===')

# ── Load model & data ──
model = xgb.Booster()
model.load_model(MODEL_PATH)
fnames = model.feature_names
print(f'Model features: {fnames}')

train = pd.read_csv(TRAIN_CSV, index_col=0)
print(f'Training data: {train.shape}')

# Extract features
x = train[fnames].values
y = train['label'].values
group = train['group'].values

# Create DMatrix
dmat = xgb.DMatrix(x, label=y, feature_names=fnames)

# Compute SHAP values via pred_contribs
shap_contrib = model.predict(dmat, pred_contribs=True)
bias = shap_contrib[:, -1]
shap_vals = shap_contrib[:, :-1]
print(f'SHAP matrix: {shap_vals.shape}')

# Sort features by mean |SHAP|
mean_abs = np.mean(np.abs(shap_vals), axis=0)
sorted_idx = np.argsort(mean_abs)[::-1]
sorted_names = [fnames[i] for i in sorted_idx]

# ── Panel A: Beeswarm ──
print('Creating beeswarm plot...')
fig_a, ax = plt.subplots(figsize=(W1, 3.8))

np.random.seed(42)
for feat_idx, fi in enumerate(sorted_idx):
    vals = shap_vals[:, fi]
    jitter = np.random.normal(0, 0.2, size=len(vals))
    for v, j in zip(vals, jitter):
        # Color by group: resistant=red, sensitive=blue
        color = '#C0392B' if v > 0 else '#3498DB'
        ax.plot(v, feat_idx + j, 'o', color=color, markersize=2.5, alpha=0.4, markeredgecolor='none')

ax.set_yticks(range(len(sorted_names)))
ax.set_yticklabels(sorted_names, fontsize=7.5)
ax.set_xlabel('SHAP value (impact on model output)', fontsize=9)
ax.axvline(0, color='gray', linestyle='--', linewidth=0.5)
ax.invert_yaxis()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(0.5)
ax.spines['bottom'].set_linewidth(0.5)

# Legend
ax.plot([], [], 'o', color='#C0392B', markersize=6, label='Risk (high output)')
ax.plot([], [], 'o', color='#3498DB', markersize=6, label='Protective (low output)')
ax.legend(fontsize=6, loc='lower left', bbox_to_anchor=(-0.5, -0.08), framealpha=0.8, edgecolor='#ccc')

out_a = os.path.join(OUT_DIR, 'Fig2A_SHAP_beeswarm.png')
fig_a.savefig(out_a, dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close(fig_a)
print(f'  Saved: Fig2A ({os.path.getsize(out_a)//1024} KB)')

# ── Panel B: Mean |SHAP| bar ──
print('Creating bar plot...')
fig_b, ax = plt.subplots(figsize=(W1, 3.8))

std_abs = np.std(np.abs(shap_vals), axis=0)
vals_sorted = mean_abs[sorted_idx]
stds_sorted = std_abs[sorted_idx]

colors = ['#C0392B' if i < 3 else '#4DBBD5' for i in range(len(sorted_idx))]
bars = ax.barh(range(len(sorted_idx)), vals_sorted, xerr=stds_sorted, color=colors,
               edgecolor='white', linewidth=0.5, capsize=3,
               error_kw={'linewidth': 0.8, 'color': '#666'})

ax.set_yticks(range(len(sorted_idx)))
ax.set_yticklabels(sorted_names, fontsize=7.5)
ax.set_xlabel('Mean |SHAP value| ± SD', fontsize=9)
ax.invert_yaxis()

for i, v in enumerate(vals_sorted):
    ax.text(v + 0.01, i - 0.18, f'{v:.3f}', va='bottom', fontsize=6.5, color='#333')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(0.5)
ax.spines['bottom'].set_linewidth(0.5)

out_b = os.path.join(OUT_DIR, 'Fig2B_SHAP_bar.png')
fig_b.savefig(out_b, dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close(fig_b)
print(f'  Saved: Fig2B ({os.path.getsize(out_b)//1024} KB)')

# ── Combined Figure (for HTML gallery) ──
print('Creating combined figure...')
fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.8))
plt.subplots_adjust(wspace=0.35, left=0.18, right=0.97)

for ax_idx, ax in enumerate(axes):
    if ax_idx == 0:  # Beeswarm
        np.random.seed(42)
        for feat_idx, fi in enumerate(sorted_idx):
            vals = shap_vals[:, fi]
            jitter = np.random.normal(0, 0.2, size=len(vals))
            for v, j in zip(vals, jitter):
                color = '#C0392B' if v > 0 else '#3498DB'
                ax.plot(v, feat_idx + j, 'o', color=color, markersize=2, alpha=0.35, markeredgecolor='none')
        ax.set_yticks(range(len(sorted_names)))
        ax.set_yticklabels(sorted_names, fontsize=7.5)
        ax.set_xlabel('SHAP value', fontsize=9)
        ax.axvline(0, color='gray', linestyle='--', linewidth=0.5)
        ax.invert_yaxis()
        ax.plot([], [], 'o', color='#C0392B', label='Resistant', markersize=5)
        ax.plot([], [], 'o', color='#3498DB', label='Sensitive', markersize=5)
        ax.legend(fontsize=6, loc='lower left', bbox_to_anchor=(-0.35, -0.08), framealpha=0.8, edgecolor='#ccc')
    else:  # Bar
        ax.barh(range(len(sorted_idx)), vals_sorted, xerr=stds_sorted, color=colors,
               edgecolor='white', linewidth=0.5, capsize=3,
               error_kw={'linewidth': 0.8, 'color': '#666'})
        ax.set_yticks(range(len(sorted_idx)))
        ax.set_yticklabels(sorted_names, fontsize=7.5)
        ax.set_xlabel('Mean |SHAP| ± SD', fontsize=9)
        ax.invert_yaxis()
        for i, v in enumerate(vals_sorted):
            ax.text(v + 0.01, i - 0.18, f'{v:.3f}', va='bottom', fontsize=6.5, color='#333')
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    fig.text(0.03 if ax_idx == 0 else 0.52, 1.06, chr(65 + ax_idx),
             fontsize=13, fontweight='bold', transform=ax.transAxes, va='top')

out_combined = os.path.join(OUT_DIR, 'Fig2_SHAP.png')
fig.savefig(out_combined, dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'  Saved: Combined ({os.path.getsize(out_combined)//1024} KB)')

print('\nDone!')
print(f'Output: {OUT_DIR}')
