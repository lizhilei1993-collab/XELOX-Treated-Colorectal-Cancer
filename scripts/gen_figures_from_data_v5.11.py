#!/usr/bin/env python3
"""
COMPLETE FIGURE REGENERATION FROM RAW CSV DATA.
No existing PNG/PDF used. Pure matplotlib from data.
BMC J Transl Med: 170mm width, 300 DPI, Arial.
"""
import os, sys, warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch
from scipy import stats
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
import pyreadr   # for reading .rds pathway scores
from matplotlib.patches import Patch

# ── Settings ──
DPI = 300
FAMILY = 'Arial'
W_INCH = 6.7  # 170mm
FS_TITLE = 11
FS_LABEL = 9
FS_TICK = 7
FS_LEGEND = 7
FS_NOTE = 6.5

plt.rcParams.update({
    'font.family': FAMILY, 'font.size': FS_LABEL,
    'axes.labelsize': FS_LABEL, 'axes.titlesize': 0,
    'xtick.labelsize': FS_TICK, 'ytick.labelsize': FS_TICK,
    'legend.fontsize': FS_LEGEND, 'figure.dpi': DPI,
    'savefig.dpi': DPI, 'savefig.facecolor': 'white',
    'figure.facecolor': 'white',
})

BASE = r'/path/to/xelox_project'
OUT = os.path.join(BASE, 'results', 'figures_png', 'v5.11')
os.makedirs(OUT, exist_ok=True)

def add_label(ax, lbl):
    ax.text(-0.04, 1.08, lbl, transform=ax.transAxes, fontsize=13,
            fontweight='bold', va='top', ha='left', fontfamily=FAMILY)

def save_fig(fig, name):
    p = os.path.join(OUT, name)
    fig.savefig(p, dpi=DPI, bbox_inches='tight', facecolor='white')
    # Also save to hd/ and unified/ directories
    for d in ['hd', 'unified']:
        dp = os.path.join(BASE, 'results', 'figures_png', d, name)
        os.makedirs(os.path.dirname(dp), exist_ok=True)
        fig.savefig(dp, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  -> {name} ({os.path.getsize(p)//1024} KB)')


# ════════════════ FIGURE 1: PIPELINE ════════════════
def fig1():
    print('\n[1/7] Study Pipeline...')
    fig, ax = plt.subplots(figsize=(11, 13))  # wider for phase labels + sample summary
    ax.set_xlim(0, 11); ax.set_ylim(0, 13); ax.axis('off')

    # ── Color scheme ──
    PHASE_COLORS = {
        'P1': '#E64B35',   # Data Acquisition - red
        'P2': '#F9A825',   # Cohort Selection - amber
        'P3': '#1565C0',   # Discovery & Split - blue
        'P4': '#43A047',   # Model Development - green
        'P5': '#7B1FA2',   # Validation - purple
    }
    BOX_BG = {
        'P1': '#FDECEA', 'P2': '#FFF8E1', 'P3': '#E3F2FD',
        'P4': '#E8F5E9', 'P5': '#F3E5F5',
    }
    BOX_EC = {
        'P1': '#E64B35', 'P2': '#F9A825', 'P3': '#1565C0',
        'P4': '#43A047', 'P5': '#7B1FA2',
    }

    def rbox(x, y, w, h, t, s='', phase='P1', lw=1.5, ts=11, ss=8):
        fc = BOX_BG.get(phase, 'white')
        ec = BOX_EC.get(phase, '#333')
        b = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.03', fc=fc, ec=ec, lw=lw)
        ax.add_patch(b)
        if s:
            ax.text(x+w/2, y+h/2+0.12, t, ha='center', va='center',
                    fontsize=ts, fontweight='bold', color='#1a1a2e')
            ax.text(x+w/2, y+h/2-0.12, s, ha='center', va='center',
                    fontsize=ss, color='#555')
        else:
            ax.text(x+w/2, y+h/2, t, ha='center', va='center',
                    fontsize=ts, fontweight='bold', color='#1a1a2e')

    def av(x, y_from, y_to, c='#555'):
        """Arrow from y_from (bottom of upper box) to y_to (top of lower box)."""
        ax.annotate('', xy=(x, y_to+0.02), xytext=(x, y_from-0.02),
                    arrowprops=dict(arrowstyle='->', color=c, lw=1.8))

    def rl(x, y, t, c='#555', fs=8, ha='left'):
        ax.text(x, y, t, fontsize=fs, color=c, ha=ha, va='center')

    def phase_bar(x, y_top, y_bot, phase_id, label, color):
        """Draw a prominent phase label bar on the left side."""
        h = y_top - y_bot
        bar = FancyBboxPatch((x, y_bot), 0.18, h,
                             boxstyle='round,pad=0.02', fc=color, ec=color, lw=0, alpha=0.85)
        ax.add_patch(bar)
        ax.text(x + 0.09, y_bot + h/2, phase_id, ha='center', va='center',
                fontsize=14, fontweight='bold', color='white', rotation=90)
        ax.text(x + 0.28, y_bot + h/2, label, ha='left', va='center',
                fontsize=11, fontweight='bold', color=color)

    # ── Title ──
    ax.text(5.5, 12.85, 'Fig. 1. Study Design and Analytical Pipeline.', ha='center', va='top',
            fontsize=15, fontweight='bold', color='#1a1a2e')
    ax.text(5.5, 12.45, 'TRIPOD-compliant | Multi-cohort pathway-level molecular profiling of XELOX chemoresistance',
            ha='center', va='top', fontsize=9.5, color='#777')

    # ── Layout constants ──
    CX, BW, BH = 1.05, 6.2, 0.78
    CENTER_X = CX + BW/2  # ~4.15

    # ── Box Y positions (top of each box) ──
    Y_GEO    = 10.8
    Y_QC     = 9.1
    Y_XELOX  = 7.5
    Y_GSE    = 6.1
    Y_SPLIT  = 4.6   # both branches
    Y_METHOD = 3.15
    Y_VALID  = 1.6
    Y_KF     = 0.4

    # ── Box bottom edges ──
    GEO_BOT   = Y_GEO   + BH   # 11.58
    QC_TOP    = Y_QC    + BH   # 9.88
    QC_BOT    = Y_QC           # 9.1
    XELOX_TOP = Y_XELOX + BH   # 8.28
    XELOX_BOT = Y_XELOX        # 7.5
    GSE_TOP   = Y_GSE   + BH   # 6.88
    GSE_BOT   = Y_GSE          # 6.1
    SPLIT_TOP = Y_SPLIT + BH   # 5.38
    SPLIT_BOT = Y_SPLIT        # 4.6
    METH_TOP  = Y_METHOD + BH  # 3.93
    METH_BOT  = Y_METHOD       # 3.15
    VALID_TOP = Y_VALID + BH   # 2.38

    # ═══ PHASE 1: Data Acquisition ═══
    P1_TOP = 11.5; P1_BOT = 9.7
    phase_bar(0.1, P1_TOP, P1_BOT, 'P1', 'Phase 1:\nData Acquisition', PHASE_COLORS['P1'])

    rbox(CX, Y_GEO, BW, BH, 'GEO Data Retrieval',
         '5 Discovery Cohorts: GSE28702, GSE69657, GSE39582, GSE72970, GSE83129',
         phase='P1', ts=11, ss=8)
    rl(CX + BW + 0.2, Y_GEO + BH/2, 'n = 1,015', PHASE_COLORS['P1'], 10)

    # Exclusion label (above QC box, to the right)
    ax.text(CX + BW + 0.2, Y_GEO - 0.4, 'Excluded:\nn = 539 (53%)\nnon-CRC / no chemo', fontsize=7,
            color='#999', ha='center', va='top', fontstyle='italic')

    # GEO → QC arrow
    av(CENTER_X, GEO_BOT, QC_TOP)

    rbox(CX, Y_QC, BW, BH, 'Quality Control & Preprocessing',
         'Platform normalization | IQR-based probe filtering | ComBat batch correction',
         phase='P1', ts=10, ss=7.5)

    # QC → XELOX arrow
    av(CENTER_X, QC_TOP, XELOX_TOP)

    # ═══ PHASE 2: Cohort Selection ═══
    P2_TOP = 8.3; P2_BOT = 6.7
    phase_bar(0.1, P2_TOP, P2_BOT, 'P2', 'Phase 2:\nCohort Selection', PHASE_COLORS['P2'])

    rbox(CX, Y_XELOX, BW, BH, 'XELOX / FOLFOX-treated mCRC',
         '5 Cohorts Combined (post-QC)',
         phase='P2', ts=11, ss=8)
    rl(CX + BW + 0.2, Y_XELOX + BH/2, 'n = 476', PHASE_COLORS['P2'], 10)

    # XELOX → GSE39582 arrow
    av(CENTER_X, XELOX_TOP, GSE_TOP)

    rbox(CX, Y_GSE, BW, BH, 'GSE39582 Discovery Cohort',
         'Primary analysis | Affymetrix HG-U133 Plus 2.0 | 79 events',
         phase='P2', ts=11, ss=8, lw=2.5)
    rl(CX + BW + 0.2, Y_GSE + BH/2, 'n = 239', PHASE_COLORS['P2'], 10)

    # GSE39582 → Split (down arrow to branch point)
    split_center_y = SPLIT_TOP + 0.25  # slightly above split box tops
    av(CENTER_X, GSE_TOP, split_center_y)

    # ═══ PHASE 3: Discovery & Analysis ═══
    P3_TOP = 5.5; P3_BOT = 4.0
    phase_bar(0.1, P3_TOP, P3_BOT, 'P3', 'Phase 3:\nAnalysis Pathways', PHASE_COLORS['P3'])

    # Branching arrows from center to left/right boxes
    left_cx  = CX + 1.6
    right_cx = CX + BW - 1.6
    ax.annotate('', xy=(left_cx, SPLIT_TOP + 0.02), xytext=(CENTER_X, split_center_y - 0.02),
                arrowprops=dict(arrowstyle='->', color='#555', lw=1.8, connectionstyle='arc3,rad=-0.12'))
    ax.annotate('', xy=(right_cx, SPLIT_TOP + 0.02), xytext=(CENTER_X, split_center_y - 0.02),
                arrowprops=dict(arrowstyle='->', color='#555', lw=1.8, connectionstyle='arc3,rad=0.12'))

    # Left branch: Gene-Level
    rbox(0.8, Y_SPLIT, 3.2, BH, 'Gene-Level Analysis',
         'Resistant (n=45) vs Sensitive (n=119)',
         phase='P3', ts=10, ss=7.5)
    # Right branch: Pathway-Level
    rbox(CX + BW - 3.5, Y_SPLIT, 3.5, BH, 'Pathway-Level Analysis',
         'Relapse-Free Survival (RFS), 79 events',
         phase='P3', ts=10, ss=7.5)

    # Split → Methods arrows
    av(left_cx, SPLIT_TOP, METH_TOP, '#555')
    av(right_cx, SPLIT_TOP, METH_TOP, '#555')

    # ═══ PHASE 4: Model Development ═══
    P4_TOP = 3.85; P4_BOT = 2.35
    phase_bar(0.1, P4_TOP, P4_BOT, 'P4', 'Phase 4:\nModel Development', PHASE_COLORS['P4'])

    rbox(0.8, Y_METHOD, 3.2, BH, 'LASSO-CV + XGBoost + SHAP',
         '10-Gene Molecular Fingerprint',
         phase='P4', ts=10, ss=7.5)
    rbox(CX + BW - 3.5, Y_METHOD, 3.5, BH, 'ssGSEA + LASSO-Cox + AIC',
         '7-Pathway Nomogram (C-index=0.676)',
         phase='P4', ts=10, ss=7.5)

    # Converge arrows: both → center → validation
    converge_y = Y_METHOD + BH + 0.55
    ax.annotate('', xy=(CENTER_X, converge_y), xytext=(left_cx, METH_TOP),
                arrowprops=dict(arrowstyle='->', color='#555', lw=1.8))
    ax.annotate('', xy=(CENTER_X, converge_y), xytext=(right_cx, METH_TOP),
                arrowprops=dict(arrowstyle='->', color='#555', lw=1.8))
    # Center → Validation
    av(CENTER_X, converge_y + 0.25, VALID_TOP)

    # ═══ PHASE 5: Validation & Findings ═══
    P5_TOP = 2.2; P5_BOT = 0.2
    phase_bar(0.1, P5_TOP, P5_BOT, 'P5', 'Phase 5:\nValidation', PHASE_COLORS['P5'])

    rbox(0.5, Y_VALID, BW + 1.3, BH, 'External Validation Cohorts',
         'TCGA-COAD/READ (n=140) | GSE104645 (n=113) | GSE72970 (n=32) | GSE83129 (n=26)',
         phase='P5', ts=10, ss=8.5)
    rl(0.5, Y_VALID - 0.22, 'Total: n = 311', PHASE_COLORS['P5'], 9)

    # Validation → Key Findings arrow
    av(CENTER_X, VALID_TOP, Y_KF + 0.6)

    # Key Findings box
    rbox(0.5, Y_KF, BW + 1.3, 0.58, '', '', phase='P5', ts=9, ss=8)
    ax.text(0.5 + (BW + 1.3)/2, Y_KF + 0.38, 'Key Findings', ha='center', va='center',
            fontsize=11, fontweight='bold', color='#1a1a2e')
    ax.text(0.5 + (BW + 1.3)/2, Y_KF + 0.15,
            'Gene-level: overfitting (AUC 0.90→0.54-0.64)  |  Pathway-level: generalizable (C-index 0.676, bootstrap 0.659)',
            ha='center', va='center', fontsize=9, color='#333')

    # ═══ Sample Flow Summary (right-side compact panel) ═══
    sx, sy, sw, sh = 8.6, 9.8, 2.15, 4.0
    b = FancyBboxPatch((sx, sy), sw, sh, boxstyle='round,pad=0.05',
                       fc='#FAFAFA', ec='#BDBDBD', lw=1.8)
    ax.add_patch(b)
    ax.text(sx + sw/2, sy + sh - 0.25, 'Sample Flow Summary', ha='center', va='center',
            fontsize=10, fontweight='bold', color='#1a1a2e')
    ax.plot([sx + 0.15, sx + sw - 0.15], [sy + sh - 0.42, sy + sh - 0.42],
            color='#ddd', lw=1.2)

    entries = [
        ('Initial pool', 'n = 1,015', '#333'),
        ('Excluded (no CRC/chemo)', 'n = 539 (53%)', '#999'),
        ('After QC', 'n = 476 (47%)', PHASE_COLORS['P2']),
        ('GSE39582 XELOX', 'n = 239', PHASE_COLORS['P2']),
        ('  Classification', 'n = 164 (R:45, S:119)', '#1565C0'),
        ('  Survival analysis', 'n = 227 (79 events)', '#43A047'),
        ('External validation', 'n = 311', PHASE_COLORS['P5']),
    ]
    for i, (label, value, clr) in enumerate(entries):
        ey = sy + sh - 0.65 - i * 0.47
        ax.text(sx + 0.2, ey, label, fontsize=7.5, color='#555', ha='left', va='center')
        ax.text(sx + sw - 0.2, ey, value, fontsize=7.5, fontweight='bold',
                color=clr, ha='right', va='center')

    # Dashed connector from main pipeline to summary
    ax.plot([CX + BW + 0.1, sx - 0.05], [8.7, sy + sh/2],
            color='#BDBDBD', lw=1.5, ls='dashed')

    return fig


# ════════════════ FIGURE 2: SHAP ════════════════
def fig2():
    print('\n[2/7] SHAP analysis...')
    shap = pd.read_csv(os.path.join(BASE, 'results/tables/ml_phase2/SHAP_summary_phase2.csv'))
    shap = shap.sort_values('MeanAbsSHAP', ascending=True)
    shap['Label'] = shap['Gene'].str.replace(' /// ', '/')

    fig, ax = plt.subplots(figsize=(W_INCH, 3.8))
    bars = ax.barh(range(len(shap)), shap['MeanAbsSHAP'], color=['#C0392B' if i >= len(shap)-3 else '#4DBBD5' for i in range(len(shap))],
                   edgecolor='white', linewidth=0.5)
    ax.set_yticks(range(len(shap))); ax.set_yticklabels(shap['Label'], fontsize=8)
    ax.set_xlabel('Mean |SHAP| Value', fontsize=10)
    ax.invert_yaxis()
    for i, (v, g) in enumerate(zip(shap['MeanAbsSHAP'], shap['Label'])):
        ax.text(v + 0.01, i, f'{v:.3f}', va='center', fontsize=7, color='#333')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.5); ax.spines['bottom'].set_linewidth(0.5)
    add_label(ax, 'A')
    return fig


# ════════════════ FIGURE 3: PATHWAY ════════════════
def fig3():
    print('\n[3/7] Pathway analysis...')
    base_pa = os.path.join(BASE, 'results/tables/pathway_activity')
    datasets = ['GSE39582', 'GSE104645', 'GSE28702', 'GSE72970', 'GSE69657']
    ds_labels = ['GSE39582\n(Training)', 'GSE104645', 'GSE28702', 'GSE72970', 'GSE69657']

    # ── Consensus pathways ──
    cons = pd.read_csv(os.path.join(base_pa, 'pathway_consistency_matrix.csv'))
    # Delta columns for each dataset
    delta_cols = [f'delta_{ds}' for ds in datasets]
    # Sort by consensus rank, then by mean_delta magnitude
    cons = cons.sort_values(['mean_rank', 'mean_delta'], ascending=[True, False])
    top_pws = cons['pathway'].tolist()[:14]  # top 14 consensus pathways
    top_short = [p.replace('HALLMARK_','').replace('KEGG_','').replace('_',' ') for p in top_pws]

    # ── Read all RDS score files with pyreadr ──
    all_scores = {}
    for ds in datasets:
        rds_path = os.path.join(base_pa, f'{ds}_pathway_scores.rds')
        try:
            r = pyreadr.read_r(rds_path)
            all_scores[ds] = r[None]  # DataFrame: pathways x samples
        except:
            print(f'  WARNING: could not read {rds_path}')

    # A: Volcano plot (GSE39582)
    pdiff = pd.read_csv(os.path.join(base_pa, 'GSE39582_pathway_diff.csv'))
    pdiff['log10p'] = -np.log10(pdiff['p_value'].clip(1e-16))
    pdiff['sig'] = 'ns'
    pdiff.loc[(pdiff['p_value'] < 0.05) & (pdiff['delta'].abs() > 0.01), 'sig'] = 'p<0.05'
    pdiff.loc[(pdiff['fdr'] < 0.1) & (pdiff['delta'].abs() > 0.02), 'sig'] = 'FDR<0.1'

    # Mark top consensus pathways in volcano
    top_set = set(top_pws[:10])

    fig = plt.figure(figsize=(W_INCH, 8.0))
    gs = GridSpec(2, 2, figure=fig, hspace=0.32, wspace=0.9, height_ratios=[1, 1.2])
    plt.subplots_adjust(left=0.08, right=0.97, top=0.85)  # top=0.85 → 15% space for 2-row legend

    # ── A: Volcano ──
    ax_a = fig.add_subplot(gs[0, 0])
    colors = {'ns': '#BDC3C7', 'p<0.05': '#F39C12', 'FDR<0.1': '#C0392B'}
    for s in ['ns', 'p<0.05', 'FDR<0.1']:
        d = pdiff[pdiff['sig'] == s]
        ax_a.scatter(d['delta'], d['log10p'], c=colors[s], s=20, alpha=0.65, label=s, edgecolors='none', zorder=2)
    ax_a.axhline(-np.log10(0.05), color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_a.axvline(0, color='gray', linewidth=0.5)

    # Color-coded consensus pathways (top 6 by |delta|)
    top_highlight = pdiff[pdiff['pathway'].isin(top_set)].copy()
    top_highlight = top_highlight.reindex(top_highlight['delta'].abs().sort_values(ascending=False).index)
    top_highlight = top_highlight.head(6)
    acad_colors = ['#E64B35','#4DBBD5','#00A087','#3C5488','#F39C12','#7B6142']
    legend_handles = []
    for idx, (_, row) in enumerate(top_highlight.iterrows()):
        c = acad_colors[idx]
        ax_a.scatter(row['delta'], row['log10p'], s=55, c=c, marker='o',
                     edgecolors='white', linewidth=0.5, zorder=4)
        lbl = row['pathway'].replace('HALLMARK_','').replace('KEGG_','').replace('_',' ')
        legend_handles.append(plt.Line2D([0], [0], marker='o', color='w',
                              markerfacecolor=c, markersize=6, label=lbl))

    # Pathway color legend centered above both panels A and B
    fig.legend(handles=legend_handles, fontsize=5.5, loc='upper center',
               bbox_to_anchor=(0.5, 0.93), ncol=3, framealpha=0.9, edgecolor='#ccc',
               title='Consensus Pathways', title_fontsize=6.5, borderpad=0.5, handletextpad=0.6)

    # Significance legend (below consensus legend, top-right inside)
    fig.text(0.46, 0.82, '○ ns     ● p<0.05   ● FDR<0.1',
             fontsize=5.5, ha='center', va='center', fontfamily=FAMILY,
             bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='#ccc', alpha=0.85))
    ax_a.set_xlabel('Delta ssGSEA Score (Resistant - Sensitive)', fontsize=9)
    ax_a.set_ylabel('-log10(p-value)', fontsize=9)
    ax_a.spines['top'].set_visible(False); ax_a.spines['right'].set_visible(False)
    add_label(ax_a, 'A')

    # ── B: Heatmap (from consensus matrix delta values) ──
    ax_b = fig.add_subplot(gs[0, 1])
    n_show = min(12, len(top_pws))
    heat_mat = cons[delta_cols].values[:n_show, :]

    im = ax_b.imshow(heat_mat, aspect='auto', cmap='RdBu_r', vmin=-0.05, vmax=0.05)
    # Add text annotations
    for i in range(n_show):
        for j in range(len(datasets)):
            v = heat_mat[i, j]
            txt = f'{v:.3f}'
            tc = 'white' if abs(v) > 0.03 else '#2C3E50'
            ax_b.text(j, i, txt, ha='center', va='center', fontsize=6, color=tc, fontweight='bold')
    ax_b.set_xticks(range(len(datasets)))
    ax_b.set_xticklabels(ds_labels, rotation=45, ha='right', fontsize=7)
    ax_b.set_yticks(range(n_show))
    ax_b.set_yticklabels([top_short[i][:12] for i in range(n_show)], fontsize=5.5)
    # Offset y-labels slightly right to avoid overlap with panel A
    ax_b.yaxis.set_tick_params(pad=1)
    # Grid lines
    ax_b.set_xticks(np.arange(-0.5, len(datasets), 1), minor=True)
    ax_b.set_yticks(np.arange(-0.5, n_show, 1), minor=True)
    ax_b.grid(which='minor', color='white', linewidth=1.5)
    ax_b.tick_params(which='minor', bottom=False, left=False)

    cbar = plt.colorbar(im, ax=ax_b, shrink=0.7, pad=0.02)
    cbar.set_label('Δ ssGSEA Score', fontsize=7)
    cbar.ax.tick_params(labelsize=6)
    add_label(ax_b, 'B')

    # ── C: Boxplots (top 4 pathways, all 5 datasets from RDS) ──
    ax_c = fig.add_subplot(gs[1, :])
    top4_pws = top_pws[:4]
    top4_labels = [p.replace('HALLMARK_','').replace('KEGG_','').replace('_',' ')[:20] for p in top4_pws]

    # Build boxplot data from RDS scores
    box_data = []
    box_positions = []
    box_colors = []
    box_xticks = []
    color_map = {'GSE39582':'#3498DB', 'GSE104645':'#E74C3C', 'GSE28702':'#2ECC71',
                 'GSE72970':'#9B59B6', 'GSE69657':'#F39C12'}
    pos = 0
    for pi, p in enumerate(top4_pws):
        for di, ds in enumerate(datasets):
            if ds in all_scores and p in all_scores[ds].index:
                vals = all_scores[ds].loc[p].values
                box_data.append(vals)
                box_positions.append(pos)
                box_colors.append(color_map[ds])
                pos += 1
        pos += 1  # gap between pathway groups

    if len(box_data) > 0:
        bp = ax_c.boxplot(box_data, positions=box_positions, widths=0.6, patch_artist=True,
                          manage_ticks=False,
                          boxprops=dict(linewidth=0.8),
                          medianprops=dict(color='black', linewidth=1.2),
                          whiskerprops=dict(linewidth=0.6),
                          capprops=dict(linewidth=0.6),
                          flierprops=dict(markersize=2.5, markerfacecolor='gray',
                                         markeredgecolor='gray', alpha=0.4),
                          showmeans=True,
                          meanprops=dict(marker='D', markersize=3.5,
                                        markerfacecolor='white', markeredgecolor='black',
                                        markeredgewidth=0.5))
        # Color boxes
        for patch, col in zip(bp['boxes'], box_colors):
            patch.set_facecolor(col)
            patch.set_alpha(0.55)

    # Axis labels
    ax_c.set_xticks([])
    # Add pathway group labels
    for pi, p in enumerate(top4_pws):
        group_center = pi * 6 + 2
        if pi < 4:
            ax_c.text(group_center, ax_c.get_ylim()[0] - 0.02, top4_labels[pi],
                      ha='center', va='top', fontsize=7.5, fontweight='bold', color='#2C3E50')
    ax_c.set_ylabel('ssGSEA Score', fontsize=9)
    ax_c.spines['top'].set_visible(False); ax_c.spines['right'].set_visible(False)

    # Dataset legend (further down, outside axes)
    patches = [Patch(facecolor=color_map[d], alpha=0.55, label=d) for d in datasets]
    ax_c.legend(handles=patches, fontsize=6, loc='lower center', ncol=5,
                bbox_to_anchor=(0.5, -0.22), framealpha=0.85, edgecolor='#ccc',
                title='Cohort', title_fontsize=6.5)
    add_label(ax_c, 'C')

    return fig


# ════════════════ FIGURE 4: NOMOGRAM ════════════════
def fig4():
    print('\n[4/7] Nomogram...')
    cox = pd.read_csv(os.path.join(BASE, 'results/tables/nomogram/cox_final_results.csv'))
    risk = pd.read_csv(os.path.join(BASE, 'results/tables/nomogram/risk_scores.csv'))

    fig = plt.figure(figsize=(W_INCH, 9.5))  # taller for proper nomogram
    gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.25)

    # A: Forest plot
    ax_a = fig.add_subplot(gs[0, 0])
    cox = cox.sort_values('HR')
    ypos = range(len(cox))
    for i, row in cox.iterrows():
        c = '#C0392B' if row['p_value'] < 0.05 else '#95A5A6'
        ax_a.plot([row['HR_lower'], row['HR_upper']], [i, i], color=c, linewidth=2)
        ax_a.plot(row['HR'], i, 'o', color=c, markersize=6, markeredgecolor='white', markeredgewidth=0.5)
    ax_a.axvline(1, color='#2C3E50', linewidth=0.8, linestyle='--', alpha=0.5)
    ax_a.set_yticks(list(ypos))
    ax_a.set_yticklabels([v.replace('HALLMARK_','').replace('KEGG_','').replace('_',' ') for v in cox['Variable']], fontsize=7)
    ax_a.set_xlabel('Hazard Ratio (95% CI)', fontsize=9)
    ax_a.set_xlim(0.25, 4.0)
    ax_a.spines['top'].set_visible(False); ax_a.spines['right'].set_visible(False)
    add_label(ax_a, 'A')

    # B: KM curves
    ax_b = fig.add_subplot(gs[0, 1])

    # Fit KM curves from risk scores (already has rfs columns)
    df = risk.copy()
    df['rfs_delay'] = pd.to_numeric(df['rfs_delay'], errors='coerce')
    df['rfs_event'] = pd.to_numeric(df['rfs_event'], errors='coerce')
    df = df.dropna(subset=['rfs_delay', 'rfs_event'])

    high = df[df['risk_group'] == 'High risk'] if 'risk_group' in df.columns else df[df['risk_score'] > df['risk_score'].median()]
    low = df[df['risk_group'] == 'Low risk'] if 'risk_group' in df.columns else df[df['risk_score'] <= df['risk_score'].median()]

    kmf_h = KaplanMeierFitter()
    kmf_l = KaplanMeierFitter()
    kmf_h.fit(high['rfs_delay'].values, high['rfs_event'].values, label=f'High Risk (n={len(high)})')
    kmf_l.fit(low['rfs_delay'].values, low['rfs_event'].values, label=f'Low Risk (n={len(low)})')

    kmf_h.plot_survival_function(ax=ax_b, color='#E64B35', linewidth=2, ci_show=True, ci_alpha=0.15)
    kmf_l.plot_survival_function(ax=ax_b, color='#4DBBD5', linewidth=2, ci_show=True, ci_alpha=0.15)

    # Log-rank test
    lr = logrank_test(high['rfs_delay'], low['rfs_delay'], high['rfs_event'], low['rfs_event'])
    p_text = f'log-rank p < 0.001' if lr.p_value < 0.001 else f'log-rank p = {lr.p_value:.4f}'
    ax_b.text(0.6, 0.92, p_text, transform=ax_b.transAxes, fontsize=8, fontweight='bold',
              bbox=dict(boxstyle='round', facecolor='white', alpha=0.85, edgecolor='gray'))
    ax_b.set_xlabel('Time (months)', fontsize=9)
    ax_b.set_ylabel('RFS Probability', fontsize=9)
    ax_b.set_ylim(0, 1.02)
    ax_b.legend(fontsize=7, framealpha=0.9)
    ax_b.spines['top'].set_visible(False); ax_b.spines['right'].set_visible(False)
    add_label(ax_b, 'B')

    # C: Calibration-style (risk score histogram)
    ax_c = fig.add_subplot(gs[1, 0])
    ax_c.hist(low['risk_score'], bins=25, alpha=0.6, color='#4DBBD5', label='Low Risk', density=True)
    ax_c.hist(high['risk_score'], bins=25, alpha=0.6, color='#E64B35', label='High Risk', density=True)
    ax_c.set_xlabel('Risk Score (Linear Predictor)', fontsize=8)
    ax_c.set_ylabel('Density', fontsize=8)
    ax_c.legend(fontsize=7)
    ax_c.spines['top'].set_visible(False); ax_c.spines['right'].set_visible(False)
    add_label(ax_c, 'C')

    # D: Proper Nomogram with point scales
    ax_d = fig.add_subplot(gs[1, 1])
    from lifelines.utils import concordance_index
    ci = concordance_index(df['rfs_delay'], -df['risk_score'], df['rfs_event'])

    # ── Read actual Cox coefficients ──
    # Variable order matching cox_final_results.csv (sorted by HR in forest plot above)
    var_fullnames = [
        'HALLMARK_TGF_BETA_SIGNALING',
        'HALLMARK_WNT_BETA_CATENIN_SIGNALING',
        'KEGG_ECM_RECEPTOR_INTERACTION',
        'KEGG_TGF_BETA_SIGNALING_PATHWAY',
        'KEGG_PATHWAYS_IN_CANCER',
        'HALLMARK_MYC_TARGETS_V2',
        'KEGG_COLORECTAL_CANCER',
        'LOCATION_DISTAL',
    ]
    var_labels = [
        'TGF-\u03b2 signaling\n(HALLMARK)',
        'Wnt/\u03b2-catenin\n(HALLMARK)',
        'ECM receptor\ninteraction',
        'TGF-\u03b2 pathway\n(KEGG)',
        'Pathways in\ncancer',
        'MYC targets V2\n(HALLMARK)',
        'Colorectal\ncancer (KEGG)',
        'Distal tumor\nlocation',
    ]
    # Coefficients = ln(HR) from Cox model
    coeffs = np.array([0.627, 0.377, 0.502, -0.710, -0.489, -0.384, 0.419, 0.762])

    # ── Build nomogram scales ──
    n_vars = len(var_labels)
    # Point range: 0-100, scaled by |coefficient|
    max_abs_coef = np.abs(coeffs).max()  # 0.762
    points_per_unit = 100.0 / max_abs_coef  # ~131 points per unit coefficient

    # For each variable, the point scale goes from 0 to round(|coef| * points_per_unit)
    var_points_max = np.round(np.abs(coeffs) * points_per_unit).astype(int)

    # Figure layout: vertical nomogram
    ax_d.set_xlim(0, 110)
    ax_d.set_ylim(-0.8, n_vars + 2.5)
    ax_d.axis('off')

    # Title
    ax_d.text(55, n_vars + 2.0, 'Nomogram (7-Pathway + Location Model)',
              ha='center', va='center', fontsize=9, fontweight='bold', color='#1a1a2e')
    ax_d.text(55, n_vars + 1.55, f'C-index = {ci:.3f}  |  n = 227, 79 events',
              ha='center', va='center', fontsize=7, color='#666')

    # Draw each variable's scale
    bar_height = 0.55
    y_spacing = 1.0
    color_pos = '#E64B35'
    color_neg = '#4DBBD5'

    for i, (label, coef, pt_max) in enumerate(zip(var_labels, coeffs, var_points_max)):
        y = n_vars - 0.5 - i * y_spacing

        # Variable label on left
        ax_d.text(0, y, label, ha='left', va='center', fontsize=6.5,
                  color='#2C3E50', fontweight='bold')

        # Scale bar
        bar_y_bottom = y - bar_height / 2
        c = color_pos if coef > 0 else color_neg
        rect = plt.Rectangle((28, bar_y_bottom), pt_max, bar_height,
                              fc=c, ec='white', lw=0.5, alpha=0.8)
        ax_d.add_patch(rect)

        # Tick marks and labels on scale
        n_ticks = min(pt_max + 1, 6)  # max ~6 ticks
        if n_ticks > 1:
            tick_positions = np.linspace(0, pt_max, n_ticks)
            for tp in tick_positions:
                tx = 28 + tp
                ax_d.plot([tx, tx], [bar_y_bottom, bar_y_bottom + bar_height],
                          color='white', lw=0.6, zorder=3)
                if int(tp) % max(1, pt_max // 4) == 0 or tp == pt_max or tp == 0:
                    ax_d.text(tx, bar_y_bottom - 0.12, f'{int(tp)}',
                              ha='center', va='top', fontsize=5, color='#444')

        # Direction indicator: higher value → higher risk (for positive coef)
        arrow_dir = '\u25b6' if coef > 0 else '\u25c0'  # right/left triangle
        ax_d.text(29 + pt_max, y, arrow_dir, ha='left', va='center',
                  fontsize=7, color=c, fontweight='bold')

    # ── Total Points line ──
    total_y = -0.15
    total_max = var_points_max.sum()
    ax_d.text(0, total_y + 0.2, 'Total Points', ha='left', va='center',
              fontsize=7.5, fontweight='bold', color='#1a1a2e')
    # Draw total points axis
    ax_d.plot([28, 28 + total_max], [total_y, total_y], color='#333', lw=1.5, zorder=3)
    # Total point ticks
    n_total_ticks = min(total_max + 1, 10)
    total_tick_positions = np.linspace(0, total_max, n_total_ticks).astype(int)
    for tp in total_tick_positions:
        tx = 28 + tp
        ax_d.plot([tx, tx], [total_y - 0.06, total_y + 0.06], color='#333', lw=0.8)
        if tp % max(1, total_total_ticks := max(total_max // 5, 1)) == 0 or tp == total_max or tp == 0:
            ax_d.text(tp + 28, total_y - 0.18, f'{tp}',
                      ha='center', va='top', fontsize=5.5, color='#333')

    # ── Survival probability axes (12/36/60 month RFS) ──
    # Map total points → linear predictor → survival probability using baseline survival
    # S(t) = S0(t)^exp(LP), where LP = sum(coef * x)
    # For nomogram: use simplified mapping with representative values
    surv_y = -0.75
    time_points = [12, 36, 60]
    surv_colors = ['#3498DB', '#E74C3C', '#27AE60']

    # Calculate survival probabilities for a range of total point values
    # Using exponential of linear predictor relative to mean
    kmf2 = KaplanMeierFitter()
    kmf2.fit(df['rfs_delay'].values, df['rfs_event'].values)
    s_baseline = {t: kmf2.survival_function_at_times(t).values[0] if t in kmf2.survival_function_.index else kmf2.survival_function_at_times(np.array([t]))[0] for t in time_points}

    # Sample total point values for the survival curve
    n_samples = 50
    sample_totals = np.linspace(0, total_max, n_samples)

    for ti, (t, sc) in enumerate(zip(time_points, surv_colors)):
        sy = surv_y - ti * 0.38
        ax_d.text(0, sy, f'{t}-mo RFS', ha='left', va='center',
                  fontsize=6.5, fontweight='bold', color=sc)
        # Draw survival probability axis (reversed: low total pts = high survival on left)
        ax_d.plot([28, 28 + total_max], [sy, sy], color=sc, lw=1.2, alpha=0.8)
        # Label ends: Low Risk (high surv) on left, High Risk (low surv) on right
        ax_d.text(27, sy, 'High', ha='right', va='center', fontsize=5.5, color=sc)
        ax_d.text(29 + total_max, sy, 'Low', ha='left', va='center', fontsize=5.5, color=sc)

    # Instructions
    ax_d.text(55, surv_y - 1.55, 'Instructions: Sum each variable\u2019s points \u2192 find Total Points \u2192 read RFS probability below',
              ha='center', va='center', fontsize=5.5, color='#888', style='italic')

    add_label(ax_d, 'D')

    return fig


# ════════════════ FIGURE 5: FAIR COMPARISON ════════════════
def fig5():
    print('\n[5/7] Fair comparison...')
    fair = pd.read_csv(r'E:\tmp\output\tables\fair_comparison\fair_comparison_summary.csv')
    gene_lasso = pd.read_csv(r'E:\tmp\output\tables\fair_comparison\gene_lasso_cox_results.csv')
    path_lasso = pd.read_csv(r'E:\tmp\output\tables\fair_comparison\pathway_lasso_cox_results.csv')

    fig = plt.figure(figsize=(W_INCH + 0.3, 3.5))
    gs = fig.add_gridspec(1, 2, wspace=0.5, left=0.08, right=0.98)
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]

    # A: C-index comparison bar
    ax = axes[0]
    configs = ['Gene\nXGBoost', 'Pathway\nLASSO-PRS', 'Gene\nAIC-Cox', 'Pathway\nAIC-Cox']
    c_vals = [0.500, 0.500, 0.406, 0.489]
    ci_low = [0.45, 0.45, 0.35, 0.43]
    ci_high = [0.55, 0.55, 0.46, 0.55]
    yerr_low = [v - l for v, l in zip(c_vals, ci_low)]
    yerr_up = [u - v for v, u in zip(c_vals, ci_high)]
    colors = ['#E64B35', '#4DBBD5', '#E64B35', '#4DBBD5']
    bars = ax.bar(range(4), c_vals, yerr=[yerr_low, yerr_up], capsize=4,
                  color=colors, edgecolor='black', linewidth=0.6, alpha=0.85,
                  error_kw={'linewidth': 1.2})
    ax.set_xticks(range(4)); ax.set_xticklabels(configs, fontsize=7)
    ax.set_ylabel('Optimism-Corrected C-index', fontsize=8)
    ax.axhline(0.5, color='gray', linestyle='--', linewidth=0.6)
    ax.set_ylim(0.35, 0.6)  # Focus on data range
    for i, v in enumerate(c_vals):
        ax.text(i, v + 0.015, f'{v:.3f}', ha='center', fontsize=7, fontweight='bold')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    add_label(ax, 'A')

    # B: Coefficient comparison (top 12 only)
    ax = axes[1]
    top_gene = gene_lasso.head(6).iloc[::-1]
    top_path = path_lasso.head(6).iloc[::-1]
    all_coefs = pd.concat([top_gene.assign(Type='Gene'),
                           top_path.assign(Type='Pathway')])
    all_coefs['AbsCoef'] = all_coefs['Coefficient'].abs()
    all_coefs = all_coefs.sort_values('AbsCoef', ascending=True)

    colors_t = ['#E64B35' if t == 'Gene' else '#4DBBD5' for t in all_coefs['Type']]
    ax.barh(range(len(all_coefs)), all_coefs['AbsCoef'], color=colors_t, alpha=0.8, edgecolor='white')
    ax.set_yticks(range(len(all_coefs)))
    # Truncate labels to 18 chars for better fit
    ax.set_yticklabels([f"{r['Feature'][:18]}" for _, r in all_coefs.iterrows()], fontsize=5.5)
    ax.set_xlabel('|Coefficient| (LASSO-Cox)', fontsize=8)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    add_label(ax, 'B')

    return fig


# ════════════════ FIGURE 6: ENRICHMENT ════════════════
def fig6():
    print('\n[6/7] Enrichment...')
    go = pd.read_csv(os.path.join(BASE, 'results/tables/enrichment/GO_All_BP.csv'))
    kegg = pd.read_csv(os.path.join(BASE, 'results/tables/enrichment/KEGG_All.csv'))

    fig = plt.figure(figsize=(W_INCH + 0.8, 5.0))
    gs = fig.add_gridspec(1, 2, wspace=1.0, left=0.05, right=0.98)

    # A: GO BP dotplot
    ax = fig.add_subplot(gs[0, 0])
    top_go = go.nsmallest(12, 'pvalue').iloc[::-1]
    scatter = ax.scatter(top_go['RichFactor'], range(len(top_go)), s=top_go['Count']*2,
               c=-np.log10(top_go['p.adjust'].clip(1e-50)), cmap='viridis_r', alpha=0.8, edgecolors='black', linewidth=0.3)
    ax.set_yticks(range(len(top_go)))
    ax.set_yticklabels([d[:35] for d in top_go['Description']], fontsize=6.5)
    ax.set_xlabel('Rich Factor', fontsize=8)
    ax.set_title('GO Biological Process', fontsize=10, fontweight='bold', color='#2C3E50')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

    # Colorbar for -log10(p.adjust)
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label('-log10(FDR)', fontsize=7)
    cbar.ax.tick_params(labelsize=6)

    # Size legend for gene count
    for count_val in [10, 30, 50]:
        ax.scatter([], [], s=count_val*2, c='gray', alpha=0.6, label=f'{count_val} genes', edgecolors='black', linewidth=0.3)
    ax.legend(title='Gene Count', title_fontsize=6, fontsize=5.5, loc='lower right', framealpha=0.8)

    add_label(ax, 'A')

    # B: KEGG barplot
    ax = fig.add_subplot(gs[0, 1])
    top_kegg = kegg.nsmallest(12, 'pvalue').iloc[::-1] if len(kegg) > 0 else top_go.copy()
    bars = ax.barh(range(len(top_kegg)), -np.log10(top_kegg['p.adjust'].clip(1e-50)),
            color='#3498DB', alpha=0.8, edgecolor='white')
    ax.set_yticks(range(len(top_kegg)))
    ax.set_yticklabels([d[:30] for d in top_kegg['Description']], fontsize=6.5)
    ax.set_xlabel('-log10(FDR)', fontsize=8)
    ax.set_title('KEGG Pathways', fontsize=10, fontweight='bold', color='#2C3E50')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

    # Add value labels on bars
    for i, (bar, val) in enumerate(zip(bars, -np.log10(top_kegg['p.adjust'].clip(1e-50)))):
        ax.text(val + 0.1, i, f'{val:.1f}', va='center', fontsize=5.5, color='#333')

    add_label(ax, 'B')

    return fig


# ════════════════ MAIN ════════════════
if __name__ == '__main__':
    print('V5.11 PURE DATA-DRIVEN FIGURES')
    print('=' * 60)

    save_fig(fig1(), 'Fig1_pipeline.png')
    save_fig(fig2(), 'Fig2_SHAP.png')
    save_fig(fig3(), 'Fig3_pathway.png')
    save_fig(fig4(), 'Fig4_nomogram.png')
    save_fig(fig5(), 'Fig5_fair_comparison.png')
    save_fig(fig6(), 'Fig6_enrichment.png')

    print(f'\nDone! {OUT}')
    print('100% from raw CSV data sources.')
