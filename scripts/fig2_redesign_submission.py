#!/usr/bin/env python3
"""
Figure 2 Redesign — FINAL v3
BMC Cancer / Scientific Reports submission.
Panel A: Ensemble Voting lollipop (descending by votes)
Panel B: XGBoost SHAP beeswarm
Two-column (170 mm), 600 DPI TIFF + editable SVG.
"""
import os, sys, warnings, io
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import xgboost as xgb
from PIL import Image

# ── Paths ──
BASE = r'/path/to/xelox_project'
MODEL_PATH   = os.path.join(BASE, 'results/tables/ml_phase2/XGBoost_model_phase2.xgb')
TRAIN_CSV    = r'E:\tmp\training_shap_data.csv'
VOTES_CSV    = os.path.join(BASE, 'results/tables/ml_phase2/MultiAlgo_feature_votes.csv')
FP_CSV       = os.path.join(BASE, 'results/tables/ml_phase2/XELOX_resistance_fingerprint.csv')
OUT_DIR      = os.path.join(BASE, 'results/figures_png/v5.11')
os.makedirs(OUT_DIR, exist_ok=True)

# ── Publication settings ──
DPI = 600
FAMILY = 'Arial'
W_INCH, H_INCH = 6.70, 6.80   # 170 mm wide, taller for bottom text

plt.rcParams.update({
    'font.family': FAMILY, 'font.size': 8,
    'axes.labelsize': 8.5, 'axes.titlesize': 9.5,
    'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
    'legend.fontsize': 6.5,
    'figure.dpi': DPI, 'savefig.dpi': DPI,
    'savefig.facecolor': 'white', 'figure.facecolor': 'white',
    'axes.linewidth': 0.6,
    'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    'xtick.major.size': 3,   'ytick.major.size': 3,
})

# ── Colors ──
C_FP     = '#D62728'
C_NON_FP = '#7F7F7F'
C_RISK   = '#E64B35'
C_PROT   = '#357ABD'
C_BG     = '#FCFCFC'

print('=== Figure 2 v3 ===')

# ══════════════════════════════════════════════════
# 1. Load data
# ══════════════════════════════════════════════════
model = xgb.Booster()
model.load_model(MODEL_PATH)
fnames = model.feature_names
train  = pd.read_csv(TRAIN_CSV, index_col=0)
x_arr  = train[fnames].values
y_arr  = train['label'].values
dmat   = xgb.DMatrix(x_arr, label=y_arr, feature_names=fnames)
shap_v = model.predict(dmat, pred_contribs=True)[:, :-1]
ma     = np.mean(np.abs(shap_v), axis=0)
so     = np.argsort(ma)[::-1]
sn     = [fnames[i] for i in so]

n_r = int(np.sum(y_arr == 1))
n_s = int(np.sum(y_arr == 0))
n_t = len(y_arr)

votes_df = pd.read_csv(VOTES_CSV)
fp_set   = set(pd.read_csv(FP_CSV)['Gene'].tolist())

# ── Top 20 DESCENDING ──
top20 = votes_df.nlargest(20, 'Votes').sort_values('Votes', ascending=False).reset_index(drop=True)
# Add fingerprint flag
top20['is_fp'] = [g in fp_set for g in top20['Gene']]

N_TOP = len(top20)
fp_count = sum(top20['is_fp'])
print(f'Top 20: {fp_count} FP genes, {N_TOP-fp_count} non-FP')
for i, row in top20.iterrows():
    g = row['Gene']
    short = (g.split(' /// ')[0] + ' (*)') if ' /// ' in g else g
    mark = ' ***' if row['is_fp'] else ''
    print(f'  y={i:2d}  v={row["Votes"]:2d}  {short}{mark}')

# Gene labels
labels = []
for g in top20['Gene']:
    labels.append((g.split(' /// ')[0] + ' (*)') if ' /// ' in g else g)

# ══════════════════════════════════════════════════
# 2. Figure
# ══════════════════════════════════════════════════
fig = plt.figure(figsize=(W_INCH, H_INCH))

gs = fig.add_gridspec(
    2, 1,
    height_ratios=[1.0, 0.85],
    hspace=0.38,
    left=0.14, right=0.88, top=0.92, bottom=0.12
)

ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[1, 0])

# ──────────────────────────────────────────────
# PANEL A: Ensemble Voting (descending)
# ──────────────────────────────────────────────
ax_a.set_facecolor(C_BG)

ypos = np.arange(N_TOP)
votes_arr = top20['Votes'].values
is_fp     = top20['is_fp'].values
max_v     = int(votes_arr.max())

# Lollipop stems
STEM_X0 = 0.1
for i in range(N_TOP):
    v  = votes_arr[i]
    fp = is_fp[i]
    c  = C_FP if fp else C_NON_FP
    lw = 2.0 if fp else 1.2
    ax_a.hlines(i, STEM_X0, v, colors=c, linewidth=lw, alpha=0.75, zorder=1)
    ax_a.scatter(v, i, s=64 if fp else 40, c=c,
                 edgecolors='white', linewidth=0.6, zorder=3, alpha=0.95)

# Vote numbers next to circles
for i in range(N_TOP):
    v  = votes_arr[i]
    fp = is_fp[i]
    ax_a.text(v + 0.4, i, str(v),
              va='center', ha='left', fontsize=7,
              fontweight='bold' if fp else 'normal',
              color=C_FP if fp else '#444')

# Y labels
ax_a.set_yticks(ypos)
ax_a.set_yticklabels(labels, fontsize=7.5, ha='right', va='center')
ax_a.tick_params(axis='y', pad=4)
for i in range(N_TOP):
    if is_fp[i]:
        ax_a.get_yticklabels()[i].set_fontweight('bold')
        ax_a.get_yticklabels()[i].set_color(C_FP)

# Threshold
ax_a.axvline(x=5, color='#888', linestyle='--', linewidth=0.7, zorder=0, alpha=0.7)
# Lower half of the dashed line (y≈15 = bottom half of chart)
ax_a.text(5.05, 15, '≥5 votes', fontsize=6.5, color='#666',
          va='top', ha='left', fontstyle='italic')

# Axes
ax_a.set_xlabel('Algorithm votes (of 11)', fontsize=8.5, labelpad=8)
ax_a.set_xlim(0, max_v + 1.5)
ax_a.set_ylim(N_TOP - 0.5, -0.5)   # y=0 at top, no extra padding
ax_a.xaxis.set_major_locator(plt.MultipleLocator(2))
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)
ax_a.spines['left'].set_linewidth(0.5)
ax_a.spines['bottom'].set_linewidth(0.5)

# Info box — bottom-right, inside axes
info_txt = (
    'Gene pool: 2,537 genes\n'
    'GDSC IC50 pre-filtered: 2,528\n'
    '11 algorithms executed\n'
    '173 genes received votes\n'
    'Top 20 shown (≥3 votes)\n'
    f'GSE39582  n={n_t}  ({n_r} R / {n_s} S)'
)
ax_a.text(0.98, 0.02, info_txt, transform=ax_a.transAxes,
          fontsize=6.2, color='#333', ha='right', va='bottom',
          linespacing=1.25,
          bbox=dict(boxstyle='round,pad=0.35', facecolor='white',
                    edgecolor='#bbb', linewidth=0.5, alpha=0.95))

# Legend — above the info box, bottom-right
leg_a = ax_a.legend(
    handles=[
        Line2D([0],[0], marker='o', color='w', markerfacecolor=C_FP,
               markersize=7, label=f'10-gene fingerprint (n={fp_count})'),
        Line2D([0],[0], marker='o', color='w', markerfacecolor=C_NON_FP,
               markersize=5, label='Not selected'),
        Line2D([0],[0], linestyle='--', color='#888', linewidth=0.7,
               label='Consensus threshold (≥5)'),
    ],
    loc='lower right', bbox_to_anchor=(0.98, 0.35),
    fontsize=6.2, framealpha=0.92,
    edgecolor='#ccc', handletextpad=0.4, handlelength=1.8, borderpad=0.3)
leg_a.get_frame().set_linewidth(0.4)

# Panel label & title
ax_a.text(-0.16, 1.05, 'A', transform=ax_a.transAxes,
          fontsize=14, fontweight='bold', va='bottom', ha='left')
ax_a.text(0.5, 1.05, 'Multi-Algorithm Ensemble Voting for Gene Selection',
          transform=ax_a.transAxes, fontsize=9.5, fontweight='bold',
          ha='center', va='bottom')

# ──────────────────────────────────────────────
# PANEL B: SHAP Beeswarm
# ──────────────────────────────────────────────
ax_b.set_facecolor(C_BG)

N_F = len(so)
np.random.seed(42)

for rank, fi in enumerate(so):
    vals = shap_v[:, fi]
    jit  = np.random.normal(0, 0.16, len(vals))
    pos, neg = vals >= 0, vals < 0
    for v, j in zip(vals[pos], jit[pos]):
        ax_b.plot(v, rank + j, 'o', color=C_RISK, markersize=2.6,
                  alpha=0.5, markeredgecolor='none', zorder=2)
    for v, j in zip(vals[neg], jit[neg]):
        ax_b.plot(v, rank + j, 'o', color=C_PROT, markersize=2.6,
                  alpha=0.5, markeredgecolor='none', zorder=2)
    ax_b.axhspan(rank - 0.45, rank + 0.45, color=C_FP, alpha=0.04, zorder=0)

ax_b.axvline(0, color='#999', linewidth=0.5, zorder=1)

# Y labels — right-aligned for clean edge
disp = [(g.split(' /// ')[0] + ' (*)') if ' /// ' in g else g for g in sn]
ax_b.set_yticks(range(N_F))
ax_b.set_yticklabels(disp, fontsize=7.5, fontweight='bold', color=C_FP,
                     ha='right', va='center')
ax_b.tick_params(axis='y', pad=4)
ax_b.set_ylim(N_F - 0.5, -0.5)

ax_b.set_xlabel('SHAP value (impact on model output)', fontsize=8.5, labelpad=8)

# Right-side |SHAP| — right-aligned on = sign, pushed to 1.07
for rank, fi in enumerate(so):
    yf = 1.0 - rank / (N_F - 1) if N_F > 1 else 0.5
    ax_b.text(1.07, yf, f'|SHAP|={ma[fi]:.3f}',
              va='center', ha='left', fontsize=6, color='#777',
              transform=ax_b.transAxes)

ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)
ax_b.spines['left'].set_linewidth(0.5)
ax_b.spines['bottom'].set_linewidth(0.5)

# Legend — bottom-left, inside axes
leg_b = ax_b.legend(
    handles=[
        Line2D([0],[0], marker='o', color='w', markerfacecolor=C_RISK,
               markersize=5, label='High expr. → Risk ↑'),
        Line2D([0],[0], marker='o', color='w', markerfacecolor=C_PROT,
               markersize=5, label='High expr. → Protective ↓'),
    ],
    loc='lower left', bbox_to_anchor=(0.0, -0.02),
    fontsize=6.5, framealpha=0.92,
    edgecolor='#ccc', handletextpad=0.3, handlelength=1.5, borderpad=0.3)
leg_b.get_frame().set_linewidth(0.4)

# Performance — above the legend, far left
ax_b.text(0.02, 0.22, 'CV AUC = 0.899\nTest AUC = 0.888  Acc = 83.7%',
          transform=ax_b.transAxes, fontsize=6.2, color='#555',
          ha='left', va='bottom', fontstyle='italic',
          linespacing=1.2)

# Panel label & title
ax_b.text(-0.16, 1.05, 'B', transform=ax_b.transAxes,
          fontsize=14, fontweight='bold', va='bottom', ha='left')
ax_b.text(0.5, 1.05, 'XGBoost SHAP Analysis of the 10-Gene Fingerprint',
          transform=ax_b.transAxes, fontsize=9.5, fontweight='bold',
          ha='center', va='bottom')

# ──────────────────────────────────────────────
# Bottom area: workflow + footnote (below x-axis labels)
# ──────────────────────────────────────────────
# Workflow connector — below Panel B's x-axis label
fig.text(0.51, 0.04,
         '11-algorithm ensemble voting  →  ≥5 consensus  →  10-gene fingerprint  →  XGBoost SHAP validation',
         ha='center', va='top', fontsize=6.5,
         fontstyle='italic', color='#444',
         bbox=dict(boxstyle='round,pad=0.25', facecolor='#FFF8EE',
                   edgecolor='none', alpha=0.9))

# Footnote — below workflow text
fig.text(0.14, 0.01,
         '(*): Multi-gene probe group — first member shown (full probe set reported in Methods)',
         ha='left', va='bottom', fontsize=5.5, color='#888', fontstyle='italic')

# ══════════════════════════════════════════════════
# 3. Save
# ══════════════════════════════════════════════════
base = 'Fig2_EnsembleVoting_SHAP_redesign'
PAD = 0.15

out_png = os.path.join(OUT_DIR, f'{base}.png')
fig.savefig(out_png, dpi=300, format='png',
            bbox_inches='tight', pad_inches=PAD,
            facecolor='white', edgecolor='none')
print(f'  PNG: {os.path.getsize(out_png)//1024} KB')

out_svg = os.path.join(OUT_DIR, f'{base}.svg')
fig.savefig(out_svg, format='svg',
            bbox_inches='tight', pad_inches=PAD,
            facecolor='white', edgecolor='none')
print(f'  SVG: {os.path.getsize(out_svg)//1024} KB')

buf = io.BytesIO()
fig.savefig(buf, dpi=DPI, format='png',
            bbox_inches='tight', pad_inches=PAD,
            facecolor='white', edgecolor='none')
buf.seek(0)
Image.open(buf).convert('RGB').save(
    os.path.join(OUT_DIR, f'{base}.tiff'),
    format='TIFF', dpi=(DPI, DPI), compression='tiff_lzw')
print(f'  TIFF: {os.path.getsize(os.path.join(OUT_DIR, f"{base}.tiff"))//1024} KB')

plt.close(fig)
print('Done.')
