#!/usr/bin/env python3
"""
Figure 1: Clean SVG study design flowchart with precise coordinate math.
All arrow endpoints calculated from box centers. No manual coordinate guessing.
"""
import os

# ── SVG canvas ──
W, H = 800, 1100
OUT_DIR = r"/path/to/xelox_project\results\figures_png\v5.11"

# ── Color palette (soft academic) ──
COLORS = {
    'stage1':   {'fill': '#FFEBEE', 'stroke': '#C62828'},  # red - discovery
    'stage2':   {'fill': '#FFF8E1', 'stroke': '#F57F17'},  # amber - QC
    'stage3':   {'fill': '#FFF3E0', 'stroke': '#E65100'},  # orange - XELOX
    'stage4':   {'fill': '#E3F2FD', 'stroke': '#1565C0'},  # blue - GSE39582
    'gene':     {'fill': '#E0F7FA', 'stroke': '#00838F'},  # teal - gene-level
    'pathway':  {'fill': '#E8F5E9', 'stroke': '#2E7D32'},  # green - pathway
    'output':   {'fill': '#F3E5F5', 'stroke': '#6A1B9A'},  # purple - validation
    'result':   {'fill': '#F1F8E9', 'stroke': '#558B2F'},  # lime - results
    'sidebar':  {'fill': '#F5F5F5', 'stroke': '#BDBDBD'},  # grey - summary
}
TXT_DARK = '#212121'
TXT_MUTED = '#555555'
TXT_LIGHT = '#999999'
ARROW_COLOR = '#616161'

# ── Box helpers ──
class Box:
    def __init__(self, x, y, w, h):
        self.x = x; self.y = y; self.w = w; self.h = h
    @property
    def cx(self): return self.x + self.w/2
    @property
    def cy(self): return self.y + self.h/2
    @property
    def top(self): return self.y
    @property
    def bot(self): return self.y + self.h
    @property
    def left(self): return self.x
    @property
    def right(self): return self.x + self.w
    @property
    def top_cx(self): return (self.cx, self.top)
    @property
    def bot_cx(self): return (self.cx, self.bot)

# ── SVG helpers ──
import xml.sax.saxutils as saxutils

def esc(t):
    return saxutils.escape(str(t))

def svg_box(box, fill, stroke, sw=1.5, rx=6):
    return f'<rect x="{box.x}" y="{box.y}" width="{box.w}" height="{box.h}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>'

def svg_text(x, y, text, size=14, weight='bold', color=TXT_DARK, anchor='middle'):
    return f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-family="Arial,Helvetica,sans-serif" font-size="{size}" font-weight="{weight}" fill="{color}">{esc(text)}</text>'

def svg_subtitle(x, y, text, size=11, color=TXT_MUTED):
    return f'<text x="{x}" y="{y}" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="{size}" fill="{color}">{esc(text)}</text>'

def svg_arrow(x1, y1, x2, y2, color=ARROW_COLOR, sw=1.5):
    """Arrow from (x1,y1) to (x2,y2) with marker-end."""
    return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{sw}" marker-end="url(#arrowhead)"/>'

def svg_label(x, y, text, color=TXT_MUTED, size=11, anchor='start'):
    return f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-family="Arial,Helvetica,sans-serif" font-size="{size}" fill="{color}">{esc(text)}</text>'

# ── Layout: single column, boxes 500px wide, centered ──
BOX_W = 500
BOX_H = 60
LEFT = (W - BOX_W) // 2   # 150
GAP = 40                  # vertical gap between boxes
BRANCH_W = 230            # branch box width
BRANCH_GAP = 40           # gap between branches

# ── Build boxes from top to bottom ──
y = 100  # starting Y after title area

# Phase label positions
phase_labels = []

# Stage 1: GEO Data Retrieval
b1 = Box(LEFT, y, BOX_W, BOX_H)
y += BOX_H + GAP
phase_labels.append((b1, 'PHASE 1: Data Acquisition'))

# Stage 2: Quality Control
b2 = Box(LEFT, y, BOX_W, BOX_H)
y += BOX_H + GAP
phase_labels.append((b2, 'PHASE 2: Preprocessing'))

# Stage 3: XELOX mCRC
b3 = Box(LEFT, y, BOX_W, BOX_H)
y += BOX_H + GAP
phase_labels.append((b3, 'PHASE 3: Cohort Assembly'))

# Stage 4: GSE39582 Discovery
b4 = Box(LEFT, y, BOX_W, BOX_H + 10)
y += BOX_H + 10 + GAP
phase_labels.append((b4, 'PHASE 4: Discovery Cohort'))

# Stage 5: Two branches
bl = Box(LEFT, y, BRANCH_W, BOX_H)
br = Box(LEFT + BRANCH_W + BRANCH_GAP, y, BRANCH_W, BOX_H)
y += BOX_H + GAP

# Stage 6: Methods
b6l = Box(LEFT, y, BRANCH_W, BOX_H)
b6r = Box(LEFT + BRANCH_W + BRANCH_GAP, y, BRANCH_W, BOX_H)
y += BOX_H + GAP

# Stage 6b: Outputs
b7l = Box(LEFT, y, BRANCH_W, BOX_H)
b7r = Box(LEFT + BRANCH_W + BRANCH_GAP, y, BRANCH_W, BOX_H)
y += BOX_H + GAP

# Stage 7: Validation (full width)
b8 = Box(LEFT, y, BOX_W, BOX_H + 10)
y += BOX_H + 10 + GAP

# Stage 8: Key Findings
b9 = Box(LEFT, y, BOX_W, BOX_H)

# Adjust SVG height
H = y + BOX_H + 60

# ── Build SVG ──
parts = []
parts.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">')
parts.append('<defs>')
parts.append('<marker id="arrowhead" markerWidth="10" markerHeight="7" refX="10" refY="3.5" orient="auto">')
parts.append('<polygon points="0 0, 10 3.5, 0 7" fill="#616161"/>')
parts.append('</marker>')
parts.append('</defs>')

# Background
parts.append(f'<rect width="{W}" height="{H}" fill="white"/>')

# Title
parts.append(svg_text(W/2, 40, 'Study Design and Analytical Pipeline', size=18, weight='bold'))
parts.append(svg_subtitle(W/2, 65, 'TRIPOD-compliant | Pathway-level molecular profiling of XELOX chemoresistance in colorectal cancer', size=11, color=TXT_LIGHT))

# ═══ ROW 1: GEO Data Retrieval ═══
parts.append(svg_box(b1, COLORS['stage1']['fill'], COLORS['stage1']['stroke']))
parts.append(svg_text(b1.cx, b1.cy - 8, 'GEO Data Retrieval', size=14))
parts.append(svg_subtitle(b1.cx, b1.cy + 14, '5 Discovery Cohorts: GSE28702, GSE69657, GSE39582, GSE72970, GSE83129'))
parts.append(svg_label(b1.right + 15, b1.cy, 'n = 1,015', color=COLORS['stage1']['stroke'], size=12))

# ═══ ROW 2: Quality Control ═══
parts.append(svg_arrow(b1.cx, b1.bot, b2.top_cx[0], b2.top))
parts.append(svg_box(b2, COLORS['stage2']['fill'], COLORS['stage2']['stroke']))
parts.append(svg_text(b2.cx, b2.cy - 8, 'Quality Control & Preprocessing', size=14))
parts.append(svg_subtitle(b2.cx, b2.cy + 14, 'Platform normalization | IQR-based probe filtering | ComBat batch correction'))
parts.append(svg_label(b2.right + 15, b2.cy - 6, 'Excluded: n = 539', color=TXT_LIGHT, size=11))
parts.append(svg_label(b2.right + 15, b2.cy + 10, '(non-CRC, no chemo info)', color='#bbbbbb', size=9))

# ═══ ROW 3: XELOX mCRC ═══
parts.append(svg_arrow(b2.cx, b2.bot, b3.top_cx[0], b3.top))
parts.append(svg_box(b3, COLORS['stage3']['fill'], COLORS['stage3']['stroke']))
parts.append(svg_text(b3.cx, b3.cy - 8, 'XELOX / FOLFOX-treated mCRC', size=14))
parts.append(svg_subtitle(b3.cx, b3.cy + 14, '5 Cohorts Combined (post-QC)'))
parts.append(svg_label(b3.right + 15, b3.cy, 'n = 476', color=COLORS['stage3']['stroke'], size=12))

# ═══ ROW 4: GSE39582 ═══
parts.append(svg_arrow(b3.cx, b3.bot, b4.top_cx[0], b4.top))
parts.append(svg_box(b4, COLORS['stage4']['fill'], COLORS['stage4']['stroke'], sw=2.5))
parts.append(svg_text(b4.cx, b4.cy - 8, 'GSE39582 Discovery Cohort', size=14))
parts.append(svg_subtitle(b4.cx, b4.cy + 14, 'Primary analysis | Affymetrix HG-U133 Plus 2.0'))
parts.append(svg_label(b4.right + 15, b4.cy, 'n = 239 (79 events)', color=COLORS['stage4']['stroke'], size=12))

# ═══ ROW 5: Branch split ═══
# Arrow from GSE39582 splits to two branches
split_y = b4.bot
split_mid_x = (bl.cx + br.cx) / 2
parts.append(svg_arrow(b4.cx, b4.bot, split_mid_x, b4.bot + GAP/2))
parts.append(svg_arrow(split_mid_x, b4.bot + GAP/2, bl.cx, bl.top))
parts.append(svg_arrow(split_mid_x, b4.bot + GAP/2, br.cx, br.top))

# Left: Gene-level
parts.append(svg_box(bl, COLORS['gene']['fill'], COLORS['gene']['stroke']))
parts.append(svg_text(bl.cx, bl.cy - 8, 'Gene-Level Analysis', size=13))
parts.append(svg_subtitle(bl.cx, bl.cy + 14, 'Resistant (n=45) vs Sensitive (n=119)', size=11))
parts.append(svg_label(bl.right + 8, bl.cy + 4, 'Excluded: 75', color=TXT_LIGHT, size=10, anchor='start'))
parts.append(svg_label(bl.right + 8, bl.cy + 20, '(intermediate RFS)', color='#cccccc', size=9, anchor='start'))

# Right: Pathway-level
parts.append(svg_box(br, COLORS['pathway']['fill'], COLORS['pathway']['stroke']))
parts.append(svg_text(br.cx, br.cy - 8, 'Pathway-Level Analysis', size=13))
parts.append(svg_subtitle(br.cx, br.cy + 14, 'Relapse-Free Survival, 79 events', size=11))

# ═══ ROW 6: Methods ═══
parts.append(svg_arrow(bl.cx, bl.bot, b6l.top_cx[0], b6l.top))
parts.append(svg_arrow(br.cx, br.bot, b6r.top_cx[0], b6r.top))

parts.append(svg_box(b6l, COLORS['gene']['fill'], COLORS['gene']['stroke']))
parts.append(svg_text(b6l.cx, b6l.cy - 8, 'LASSO-CV + XGBoost + SHAP', size=13))
parts.append(svg_subtitle(b6l.cx, b6l.cy + 14, 'Ensemble feature selection & model interpretation', size=11))

parts.append(svg_box(b6r, COLORS['pathway']['fill'], COLORS['pathway']['stroke']))
parts.append(svg_text(b6r.cx, b6r.cy - 8, 'ssGSEA + LASSO-Cox + AIC', size=13))
parts.append(svg_subtitle(b6r.cx, b6r.cy + 14, '44-pathway scoring and backward selection', size=11))

# ═══ ROW 7: Outputs ═══
parts.append(svg_arrow(b6l.cx, b6l.bot, b7l.top_cx[0], b7l.top))
parts.append(svg_arrow(b6r.cx, b6r.bot, b7r.top_cx[0], b7r.top))

parts.append(svg_box(b7l, COLORS['gene']['fill'], COLORS['gene']['stroke'], sw=1.0))
parts.append(svg_text(b7l.cx, b7l.cy - 8, '10-Gene Fingerprint', size=13))
parts.append(svg_subtitle(b7l.cx, b7l.cy + 14, 'Discovery AUC = 0.899 | External AUC = 0.54\u20130.64', size=10))

parts.append(svg_box(b7r, COLORS['pathway']['fill'], COLORS['pathway']['stroke'], sw=1.0))
parts.append(svg_text(b7r.cx, b7r.cy - 8, '7-Pathway Nomogram', size=13, weight='bold'))
parts.append(svg_subtitle(b7r.cx, b7r.cy + 14, 'C-index = 0.676 | Bootstrap-corrected = 0.659', size=10, color=TXT_DARK))

# ═══ ROW 8: Validation ═══
parts.append(svg_arrow(b7l.cx, b7l.bot, b8.cx, b8.top))
parts.append(svg_arrow(b7r.cx, b7r.bot, b8.cx, b8.top))
parts.append(svg_box(b8, COLORS['output']['fill'], COLORS['output']['stroke'], sw=2.0))
parts.append(svg_text(b8.cx, b8.cy - 8, 'External Validation Cohorts', size=14))
parts.append(svg_subtitle(b8.cx, b8.cy + 14, 'TCGA-COAD/READ (n=140)  |  GSE104645 (n=113)  |  GSE72970 (n=32)  |  GSE83129 (n=26)', size=11))
parts.append(svg_label(LEFT, b8.bot + 18, 'Total: n = 311', color=COLORS['output']['stroke'], size=11, anchor='start'))

# ═══ ROW 9: Results ═══
parts.append(svg_arrow(b8.cx, b8.bot + 22, b9.cx, b9.top))
parts.append(svg_box(b9, COLORS['result']['fill'], COLORS['result']['stroke'], sw=2.0))
parts.append(svg_text(b9.cx, b9.cy - 8, 'Key Findings', size=14, weight='bold'))
parts.append(svg_subtitle(b9.cx, b9.cy + 14, 'Gene-level: overfitting (AUC drop 0.90 to 0.54)  |  Pathway-level: generalizable (C-index 0.659)  |  GDSC: 5-FU & capecitabine P<0.0001', size=11, color=TXT_DARK))

# ═══ Sidebar: Sample Flow Summary ═══
sidebar_x = b1.right + 80
sidebar_w = 170
sidebar_y = b1.y
sidebar_h = b4.bot - b1.y

parts.append(f'<rect x="{sidebar_x}" y="{sidebar_y}" width="{sidebar_w}" height="{sidebar_h}" rx="8" fill="{COLORS["sidebar"]["fill"]}" stroke="{COLORS["sidebar"]["stroke"]}" stroke-width="1.0"/>')
parts.append(svg_text(sidebar_x + sidebar_w/2, sidebar_y + 22, 'Sample Flow', size=12, weight='bold'))

# Dashed connector from main pipeline to sidebar
parts.append(f'<line x1="{b1.right}" y1="{b1.cy}" x2="{sidebar_x}" y2="{sidebar_y + sidebar_h/2}" stroke="#BDBDBD" stroke-width="1" stroke-dasharray="5,5"/>')

entries = [
    'Initial pool          n = 1,015',
    'After QC              n = 476 (47%)',
    'GSE39582              n = 239',
    '  Classification      n = 164',
    '  Survival            n = 227',
    'External validation   n = 311',
]
for i, entry in enumerate(entries):
    parts.append(svg_label(sidebar_x + 10, sidebar_y + 50 + i * 24, entry, color=TXT_MUTED, size=10, anchor='start'))

# ── Phase indicators on the left ──
phase_data = [
    (b1.y + BOX_H/2, 'Data\nAcquisition'),
    (b2.y + BOX_H/2, 'Prepro-\ncessing'),
    (b3.y + BOX_H/2, 'Cohort\nAssembly'),
    (b4.y + BOX_H/2, 'Discovery\nCohort'),
    (bl.y + BOX_H/2, 'Feature\nEngineering'),
    (b6l.y + BOX_H/2, 'Model\nTraining'),
    (b7l.y + BOX_H/2, 'Model\nOutput'),
    (b8.y + BOX_H/2, 'Validation'),
    (b9.y + BOX_H/2, 'Conclusion'),
]
parts.append('<g opacity="0.15">')
for yy, label in phase_data:
    parts.append(svg_text(60, yy, label, size=11, weight='bold', color=TXT_MUTED))
parts.append('</g>')

parts.append('</svg>')

# ── Write SVG ──
os.makedirs(OUT_DIR, exist_ok=True)
svg_path = os.path.join(OUT_DIR, 'Fig1_pipeline.svg')
with open(svg_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(parts))
print(f'SVG written: {svg_path} ({os.path.getsize(svg_path)} bytes)')

# ── Also write to HD and unified ──
for d in ['hd', 'unified']:
    dp = os.path.join(r'/path/to/xelox_project\results\figures_png', d, 'Fig1_pipeline.svg')
    os.makedirs(os.path.dirname(dp), exist_ok=True)
    with open(dp, 'w', encoding='utf-8') as f:
        f.write('\n'.join(parts))
    print(f'SVG copied to: {dp}')
