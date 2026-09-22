"""Create Figure 1: Multi-phase analytical pipeline flowchart.
   Redesigned with strict grid alignment and uniform box sizes."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import os

OUT = r"/path/to/xelox_project\results\figures_png"

# --- Layout constants (grid system) ---
FIG_W = 24
FIG_H = 22
BOX_W = 3.8
BOX_H = 0.95
H_GAP = 0.6       # horizontal gap between boxes
V_PHASE = 2.6     # vertical gap between phases
ROW_H = 1.4       # height for two-row algorithm section
PHASE_LABEL_X = 0.35
TITLE_Y = 21.5
HEADER_FONT = 13
BOX_FONT = 10
SUB_FONT = 8

fig, ax = plt.subplots(1, 1, figsize=(FIG_W, FIG_H))
ax.set_xlim(0, FIG_W)
ax.set_ylim(0, FIG_H)
ax.axis('off')

COLORS = {
    'data': '#4472C4',
    'gene': '#548235',
    'ml': '#ED7D31',
    'validation': '#C55A11',
    'pathway': '#7030A0',
    'prs': '#2E75B6',
    'meta': '#C00000',
}

def draw_box(x, y, w, h, main_text, color, sub_text=None, fs=BOX_FONT):
    """Draw a rounded rectangle box with optional subtext below main text."""
    box = FancyBboxPatch((x, y), w, h,
                         boxstyle="round,pad=0.12", facecolor=color,
                         edgecolor='white', linewidth=2, alpha=0.9)
    ax.add_patch(box)
    y_main = y + h/2 + (0.12 if sub_text else 0)
    ax.text(x + w/2, y_main, main_text, ha='center', va='center',
            fontsize=fs, fontweight='bold', color='white')
    if sub_text:
        ax.text(x + w/2, y + h/2 - 0.28, sub_text, ha='center', va='top',
                fontsize=SUB_FONT, color='white', alpha=0.92)

def phase_header(y_center, text, color):
    """Draw phase label circle and header text."""
    ax.text(PHASE_LABEL_X, y_center, 'P' + text[0], ha='center', va='center',
            fontsize=9, fontweight='bold', color='#555',
            bbox=dict(boxstyle="circle,pad=0.22", facecolor='#EDEDED',
                      edgecolor='#BBBBBB'))
    ax.text(FIG_W/2, y_center + 0.45, text,
            ha='center', va='bottom', fontsize=HEADER_FONT,
            fontweight='bold', color=color)

def arrow_down(cx, y_from, y_to):
    """Vertical downward arrow from (cx, y_from) to (cx, y_to)."""
    ax.annotate('', xy=(cx, y_to), xytext=(cx, y_from),
                arrowprops=dict(arrowstyle='->', color='#777', lw=1.8))

def arrow_down_to_box(cx, y_from, x_target, y_target_top):
    """Arrow from center point to top of target box (L-shape: down then horizontal)."""
    mid_y = (y_from + y_target_top) / 2.
    # vertical segment
    ax.annotate('', xy=(cx, y_target_top + 0.15), xytext=(cx, y_from - 0.1),
                arrowprops=dict(arrowstyle='->', color='#888', lw=1.5))

# ================================================================
# TITLE
# ================================================================
ax.text(FIG_W/2, TITLE_Y,
        'Analytical Pipeline: XELOX Resistance Molecular Fingerprint Discovery',
        ha='center', va='center', fontsize=18, fontweight='bold', color='#333')

# ================================================================
# PHASE 1: Data Sources & Preprocessing
# ================================================================
P1_Y = 19.0
phase_header(P1_Y + 0.7, 'Phase 1: Data Sources & Preprocessing', COLORS['data'])

n_data = 5
data_total_w = n_data * BOX_W + (n_data - 1) * H_GAP
data_x_start = (FIG_W - data_total_w) / 2.
data_sources = [
    ('GSE39582\n(XELOX, n=164)'),
    ('GSE104645\n(Oxaliplatin, n=193)'),
    ('GSE28702\n(mFOLFOX6, n=83)'),
    ('GSE72970\n(FOLFOX, n=124)'),
    ('GSE69657\n(XELOX, n=30)')
]
data_centers = []
for i, label in enumerate(data_sources):
    bx = data_x_start + i * (BOX_W + H_GAP)
    draw_box(bx, P1_Y, BOX_W, BOX_H, label, COLORS['data'])
    data_centers.append(bx + BOX_W / 2.)

# TCGA box (wide, centered below data sources)
tcga_w = 14
tcga_x = (FIG_W - tcga_w) / 2.
TCGA_Y = P1_Y - 1.5
draw_box(tcga_x, TCGA_Y, tcga_w, BOX_H, 'TCGA-COAD/READ (RNA-seq, n=647)',
         '#5B9BD5', 'Cross-platform Validation')

# Arrows: data -> TCGA
for cx in data_centers:
    arrow_down(cx, P1_Y, TCGA_Y + BOX_H)

# Arrow: TCGA -> Phase 2
arrow_down(FIG_W/2, TCGA_Y, TCGA_Y - 1.2)

# ================================================================
# PHASE 2: Gene Pool Construction
# ================================================================
P2_Y = TCGA_Y - 2.2
phase_header(P2_Y + 0.7, 'Phase 2: Gene Pool Construction', COLORS['gene'])

p2_labels = [
    ('Differential Expression', 'limma |logFC|>0.5 adj.p<0.05'),
    ('WGCNA Co-expression Network', 'beta=6, scale-free R^2>0.85'),
    ('Capecitabine Metabolism Genes', 'CES1/2, TYMS, DPYD, DCK'),
]
n_p2 = len(p2_labels)
p2_total = n_p2 * BOX_W + (n_p2 - 1) * H_GAP
p2_x_start = (FIG_W - p2_total) / 2.
p2_boxes_right = []   # track for arrows
for i, (main, sub) in enumerate(p2_labels):
    bx = p2_x_start + i * (BOX_W + H_GAP)
    draw_box(bx, P2_Y, BOX_W, BOX_H, main, COLORS['gene'], sub)
    p2_boxes_right.append(bx)

# Candidate genes box (rightmost)
cand_w = 3.5
cand_x = p2_x_start + p2_total + 0.5
draw_box(cand_x, P2_Y, cand_w, BOX_H, '2,537\nCandidate Genes', '#70AD47',
         'Multi-tier integration')
p2_cand_cx = cand_x + cand_w / 2.

# Central arrow down
arrow_down(FIG_W/2, P2_Y, P2_Y - V_PHASE + 0.3)

# ================================================================
# PHASE 3: 10-Algorithm Ensemble Feature Selection
# ================================================================
P3_TOP = P2_Y - V_PHASE
phase_header(P3_TOP + 0.7, 'Phase 3: 10-Algorithm Ensemble + GDSC Dual-Drug Filter', COLORS['ml'])

algos_row1 = ['LASSO', 'Ridge', 'Elastic Net', 'Random Forest', 'XGBoost(Gain)']
algos_row2 = ['XGBoost(SHAP)', 'Boruta(Conf)', 'Boruta(Top15)', 'LogReg(AIC)', 'PLS-DA']

algo_w = 2.2
algo_h = 0.65
n_algos = len(algos_row1)
algo_total = n_algos * algo_w + (n_algos - 1) * H_GAP
algo_x_start = (FIG_W - algo_total) / 2.

P3_ROW1_Y = P3_TOP - 0.2
for i, name in enumerate(algos_row1):
    draw_box(algo_x_start + i * (algo_w + H_GAP), P3_ROW1_Y, algo_w, algo_h,
             name, COLORS['ml'], fs=9)

P3_ROW2_Y = P3_ROW1_Y - algo_h - 0.25
for i, name in enumerate(algos_row2):
    draw_box(algo_x_start + i * (algo_w + H_GAP), P3_ROW2_Y, algo_w, algo_h,
             name, COLORS['ml'], fs=9)

# GDSC filter bar
gdsc_y = P3_ROW2_Y - 0.9
gdsc_w = algo_total + 1.0
gdsc_x = (FIG_W - gdsc_w) / 2.
draw_box(gdsc_x, gdsc_y, gdsc_w, 0.7, 'GDSC Drug Sensitivity Pre-Filter',
         '#F4B183', 'Oxaliplatin & 5-FU IC50  Spearman|rho|>0.2')

# Fingerprint result
fp_y = gdsc_y - 1.25
draw_box((FIG_W - 11)/2., fp_y, 11, 0.9, '> or = 5/10 Algorithms => 10-Gene Fingerprint',
         COLORS['ml'], '+ SHAP Interpretability (Tree SHAP)')

arrow_down(FIG_W/2, fp_y, fp_y - V_PHASE + 0.3)

# ================================================================
# PHASES 4-5: Multi-Cohort Validation + ComBat Batch Correction
# ================================================================
P4_TOP = fp_y - V_PHASE
phase_header(P4_TOP + 0.7, 'Phases 4-5: Multi-Cohort Validation + ComBat Correction', COLORS['validation'])

p4_labels = [
    ('Same-Platform Validation', 'GPL570: GSE72970, GSE69657'),
    ('Cross-Platform Validation', 'TCGA RNA-seq (n=53)'),
    ('ComBat Batch Correction', 'GPL570 cohorts harmonized'),
    ('Result', 'AUC 0.54-0.64 Limited Gen.'),
]
n_p4 = len(p4_labels)
p4_total = n_p4 * BOX_W + (n_p4 - 1) * H_GAP
p4_x_start = (FIG_W - p4_total) / 2.
P4_Y = P4_TOP - 0.2
for i, (main, sub) in enumerate(p4_labels):
    bx = p4_x_start + i * (BOX_W + H_GAP)
    clr = COLORS['validation'] if i < 3 else '#C55A11'
    draw_box(bx, P4_Y, BOX_W, BOX_H, main, clr, sub if i < 3 else None)

arrow_down(FIG_W/2, P4_Y, P4_Y - V_PHASE + 0.3)

# ================================================================
# PHASES 6-8: Pathway-Level PRS -> Cox Regression -> Nomogram
# ================================================================
P5_TOP = P4_Y - V_PHASE
phase_header(P5_TOP + 0.7, 'Phases 6-8: Pathway-Level PRS -> Cox Regression -> Nomogram', COLORS['prs'])

p5_labels = [
    ('44 Oxaliplatin-Resistance Pathways', 'KEGG + Hallmark gene sets'),
    ('ssGSEA Pathway Activity Scoring', 'GSVA v2.6.1 across 5 cohorts'),
    ('3-Pathway LASSO-PRS', 'Chemokine / ECM / BER'),
    ('Cox + Nomogram', '8-var model C-index=0.676'),
]
p5_clrs = [COLORS['pathway'], COLORS['prs'], COLORS['prs'], '#2E75B6']
P5_Y = P5_TOP - 0.2
for i, ((main, sub), c) in enumerate(zip(p5_labels, p5_clrs)):
    bx = p4_x_start + i * (BOX_W + H_GAP)
    draw_box(bx, P5_Y, BOX_W, BOX_H, main, c, sub)

arrow_down(FIG_W/2, P5_Y, P5_Y - V_PHASE + 0.3)

# ================================================================
# PHASES 9-10: Meta-Analysis + Enrichment Characterization
# ================================================================
P6_TOP = P5_Y - V_PHASE
phase_header(P6_TOP + 0.7, 'Phases 9-10: Meta-Analysis + Enrichment Characterization', COLORS['meta'])

p6_labels = [
    ('Stouffer Meta-Analysis', '5 Cohorts => 253 Genes FDR<0.05'),
    ('ORA + GSEA: 5 Core Themes', 'Immune/EMT/Metabolic/Signal/DDR'),
    ('Systems Biology Landscape', 'of XELOX Resistance'),
]
n_p6 = len(p6_labels)
p6_total = n_p6 * BOX_W + (n_p6 - 1) * H_GAP
p6_x_start = (FIG_W - p6_total) / 2.
P6_Y = P6_TOP - 0.2
for i, (main, sub) in enumerate(p6_labels):
    bx = p6_x_start + i * (BOX_W + H_GAP)
    draw_box(bx, P6_Y, BOX_W, BOX_H, main, COLORS['meta'], sub)

plt.tight_layout()
out_path = os.path.join(OUT, 'Fig1_pipeline.png')
plt.savefig(out_path, dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print(f"[OK] Figure 1 created: {out_path} ({os.path.getsize(out_path)/1024:.1f} KB)")
