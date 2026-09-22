#!/usr/bin/env python3
"""
Figure 1 v6: Large-canvas, publication-quality study design pipeline.
1600x2200 canvas, minimum 14pt text, thick arrows, clean academic style.
"""
import os, xml.sax.saxutils as saxutils

W, H = 1800, 2200
OUT = r"/path/to/xelox_project\results\figures_png\v5.11"

def esc(t): return saxutils.escape(str(t))

def svg_box(x, y, w, h, fill, stroke, sw=2, rx=8):
    return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>'

def svg_text(x, y, text, size=20, weight='bold', color='#212121', anchor='middle'):
    return f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-family="Arial,sans-serif" font-size="{size}" font-weight="{weight}" fill="{color}">{esc(text)}</text>'

def svg_sub(x, y, text, size=16, color='#555555', anchor='middle'):
    return f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-family="Arial,sans-serif" font-size="{size}" fill="{color}">{esc(text)}</text>'

def svg_arrow(x1, y1, x2, y2, color='#444', sw=2.5):
    return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{sw}" marker-end="url(#arrowhead)"/>'

def svg_side(x, y, text, color='#555', size=16):
    return f'<text x="{x}" y="{y}" text-anchor="start" font-family="Arial,sans-serif" font-size="{size}" font-weight="bold" fill="{color}">{esc(text)}</text>'

# ── Layout constants ──
BOX_W = 800
BOX_H = 110
CX = (W - BOX_W) // 2   # 400
GAP = 80
BRANCH_W = 380
BRANCH_GAP = 40

parts = []
parts.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">')
parts.append('<defs>')
parts.append('<marker id="arrowhead" markerWidth="14" markerHeight="10" refX="14" refY="5" orient="auto">')
parts.append('<polygon points="0 0, 14 5, 0 10" fill="#444444"/>')
parts.append('</marker>')
parts.append('</defs>')
parts.append(f'<rect width="{W}" height="{H}" fill="white"/>')

# ── Title ──
parts.append(svg_text(W/2, 60, 'Study Design and Analytical Pipeline', size=30, weight='bold'))
parts.append(svg_sub(W/2, 100, 'TRIPOD-compliant | Multi-cohort pathway-level molecular profiling of XELOX chemoresistance', size=18, color='#777'))

# ── Phase sidebar (prominent labels with wide color bar) ──
phases = [
    (180, 340, '#E53935', 'Phase 1: Data Acquisition'),
    (430, 550, '#FB8C00', 'Phase 2: Quality Control'),
    (700, 910, '#1E88E5', 'Phase 3: Discovery & Analysis'),
    (980, 1190, '#43A047', 'Phase 4: Model Development'),
    (1260, 1480, '#8E24AA', 'Phase 5: Validation & Conclusion'),
]
for py1, py2, clr, label in phases:
    # Wide color bar background
    parts.append(f'<rect x="60" y="{py1}" width="20" height="{py2-py1}" rx="4" fill="{clr}" opacity="0.8"/>')
    # Bold horizontal label right of bar
    parts.append(f'<text x="95" y="{py1+18}" text-anchor="start" font-family="Arial,sans-serif" font-size="15" font-weight="bold" fill="{clr}">{esc(label)}</text>')

# ═══ ROW 1: GEO Data Retrieval ═══
y = 180
parts.append(svg_box(CX, y, BOX_W, BOX_H, '#FFEBEE', '#C62828'))
parts.append(svg_text(CX+BOX_W/2, y+38, 'GEO Data Retrieval', size=22))
parts.append(svg_sub(CX+BOX_W/2, y+70, '5 Discovery Cohorts: GSE28702, GSE69657, GSE39582, GSE72970, GSE83129', size=14, color='#666'))
parts.append(svg_side(CX+BOX_W+30, y+40, 'n = 1,015', color='#C62828', size=18))

# Arrow
arrow_y1 = y + BOX_H
arrow_y2 = arrow_y1 + GAP
parts.append(svg_arrow(CX+BOX_W/2, arrow_y1, CX+BOX_W/2, arrow_y2))

# ═══ ROW 2: Quality Control ═══
y = arrow_y2
parts.append(svg_box(CX, y, BOX_W, BOX_H, '#FFF8E1', '#F57F17'))
parts.append(svg_text(CX+BOX_W/2, y+38, 'Quality Control & Preprocessing', size=22))
parts.append(svg_sub(CX+BOX_W/2, y+70, 'Platform normalization | IQR-based probe filtering | ComBat batch correction', size=14, color='#666'))
parts.append(svg_side(CX+BOX_W+30, y+30, 'Excluded: n = 539', color='#999', size=15))
parts.append(svg_side(CX+BOX_W+30, y+55, '(non-CRC, no chemo info)', color='#bbb', size=12))

arrow_y1 = y + BOX_H
arrow_y2 = arrow_y1 + GAP
parts.append(svg_arrow(CX+BOX_W/2, arrow_y1, CX+BOX_W/2, arrow_y2))

# ═══ ROW 3: XELOX mCRC ═══
y = arrow_y2
parts.append(svg_box(CX, y, BOX_W, BOX_H, '#FFF3E0', '#E65100'))
parts.append(svg_text(CX+BOX_W/2, y+38, 'XELOX / FOLFOX-treated mCRC', size=22))
parts.append(svg_sub(CX+BOX_W/2, y+70, '5 Cohorts Combined (post-QC)', size=14, color='#666'))
parts.append(svg_side(CX+BOX_W+30, y+40, 'n = 476', color='#E65100', size=18))

arrow_y1 = y + BOX_H
arrow_y2 = arrow_y1 + GAP
parts.append(svg_arrow(CX+BOX_W/2, arrow_y1, CX+BOX_W/2, arrow_y2))

# ═══ ROW 4: GSE39582 Discovery ═══
y = arrow_y2
parts.append(svg_box(CX, y, BOX_W, BOX_H+20, '#E3F2FD', '#1565C0', sw=3))
parts.append(svg_text(CX+BOX_W/2, y+42, 'GSE39582 Discovery Cohort', size=24))
parts.append(svg_sub(CX+BOX_W/2, y+78, 'Primary analysis | Affymetrix HG-U133 Plus 2.0', size=14, color='#666'))
parts.append(svg_side(CX+BOX_W+30, y+45, 'n = 239 (79 events)', color='#1565C0', size=18))

# Arrow down then split
arrow_y1 = y + BOX_H + 20
split_y = arrow_y1 + GAP//2
parts.append(svg_arrow(CX+BOX_W/2, arrow_y1, CX+BOX_W/2, split_y))
# Horizontal bar
left_cx = CX + BRANCH_W//2
right_cx = CX + BOX_W - BRANCH_W//2
parts.append(svg_arrow(CX+BOX_W/2, split_y, left_cx, split_y + GAP//2))
parts.append(svg_arrow(CX+BOX_W/2, split_y, right_cx, split_y + GAP//2))

# ═══ ROW 5: Two Branches ═══
y = split_y + GAP//2

# Left: Gene-Level
parts.append(svg_box(CX, y, BRANCH_W, BOX_H, '#E0F7FA', '#00838F'))
parts.append(svg_text(CX+BRANCH_W//2, y+38, 'Gene-Level Analysis', size=20))
parts.append(svg_sub(CX+BRANCH_W//2, y+70, 'Resistant (n=45) vs Sensitive (n=119)', size=14, color='#666'))

# Right: Pathway-Level
parts.append(svg_box(CX+BOX_W-BRANCH_W, y, BRANCH_W, BOX_H, '#E8F5E9', '#2E7D32'))
parts.append(svg_text(CX+BOX_W-BRANCH_W//2, y+38, 'Pathway-Level Analysis', size=20))
parts.append(svg_sub(CX+BOX_W-BRANCH_W//2, y+70, 'Relapse-Free Survival, 79 events', size=14, color='#666'))

# Side notes (excluded info is in Sidebar already)

# Arrows down
arrow_y1 = y + BOX_H
arrow_y2 = arrow_y1 + GAP
parts.append(svg_arrow(CX+BRANCH_W//2, arrow_y1, CX+BRANCH_W//2, arrow_y2))
parts.append(svg_arrow(CX+BOX_W-BRANCH_W//2, arrow_y1, CX+BOX_W-BRANCH_W//2, arrow_y2))

# ═══ ROW 6: Methods ═══
y = arrow_y2
parts.append(svg_box(CX, y, BRANCH_W, BOX_H, '#E0F7FA', '#00838F'))
parts.append(svg_text(CX+BRANCH_W//2, y+38, 'LASSO-CV + XGBoost + SHAP', size=18))
parts.append(svg_sub(CX+BRANCH_W//2, y+70, 'Ensemble feature selection and model interpretation', size=13, color='#666'))

parts.append(svg_box(CX+BOX_W-BRANCH_W, y, BRANCH_W, BOX_H, '#E8F5E9', '#2E7D32'))
parts.append(svg_text(CX+BOX_W-BRANCH_W//2, y+38, 'ssGSEA + LASSO-Cox + AIC', size=18))
parts.append(svg_sub(CX+BOX_W-BRANCH_W//2, y+70, '44-pathway scoring and backward selection', size=13, color='#666'))

# Arrows down
arrow_y1 = y + BOX_H
arrow_y2 = arrow_y1 + GAP
parts.append(svg_arrow(CX+BRANCH_W//2, arrow_y1, CX+BRANCH_W//2, arrow_y2))
parts.append(svg_arrow(CX+BOX_W-BRANCH_W//2, arrow_y1, CX+BOX_W-BRANCH_W//2, arrow_y2))

# ═══ ROW 7: Model Outputs ═══
y = arrow_y2
parts.append(svg_box(CX, y, BRANCH_W, BOX_H, '#E0F7FA', '#00838F', sw=1.5))
parts.append(svg_text(CX+BRANCH_W//2, y+35, '10-Gene Fingerprint', size=20))
parts.append(svg_sub(CX+BRANCH_W//2, y+65, 'Discovery AUC = 0.899', size=14, color='#666'))
parts.append(svg_sub(CX+BRANCH_W//2, y+85, 'External AUC = 0.54-0.64', size=13, color='#999'))

parts.append(svg_box(CX+BOX_W-BRANCH_W, y, BRANCH_W, BOX_H, '#E8F5E9', '#2E7D32', sw=1.5))
parts.append(svg_text(CX+BOX_W-BRANCH_W//2, y+35, '7-Pathway Nomogram', size=20, weight='bold'))
parts.append(svg_sub(CX+BOX_W-BRANCH_W//2, y+65, 'C-index = 0.676', size=14, color='#333'))
parts.append(svg_sub(CX+BOX_W-BRANCH_W//2, y+85, 'Bootstrap-corrected = 0.659', size=13, color='#666'))

# Arrows to validation box (converge)
arrow_y1 = y + BOX_H
validate_y = arrow_y1 + GAP
parts.append(svg_arrow(CX+BRANCH_W//2, arrow_y1, CX+BOX_W//2, validate_y))
parts.append(svg_arrow(CX+BOX_W-BRANCH_W//2, arrow_y1, CX+BOX_W//2, validate_y))

# ═══ ROW 8: External Validation ═══
y = validate_y
parts.append(svg_box(CX, y, BOX_W, BOX_H+20, '#F3E5F5', '#6A1B9A', sw=3))
parts.append(svg_text(CX+BOX_W/2, y+42, 'External Validation Cohorts', size=24))
parts.append(svg_sub(CX+BOX_W/2, y+78, 'TCGA-COAD/READ (n=140)  |  GSE104645 (n=113)  |  GSE72970 (n=32)  |  GSE83129 (n=26)', size=14, color='#666'))
parts.append(svg_side(CX, y+BOX_H+20+18, 'Total: n = 311', color='#6A1B9A', size=16))

# Arrow to findings
arrow_y1 = y + BOX_H + 40
findings_y = arrow_y1 + GAP//2
parts.append(svg_arrow(CX+BOX_W//2, arrow_y1, CX+BOX_W//2, findings_y))

# ═══ ROW 9: Key Findings ═══
y = findings_y
parts.append(svg_box(CX, y, BOX_W, BOX_H, '#F1F8E9', '#558B2F', sw=3))
parts.append(svg_text(CX+BOX_W/2, y+38, 'Key Findings', size=24, weight='bold'))
parts.append(svg_sub(CX+BOX_W/2, y+65, 'Gene-level: overfitting (AUC 0.90 to 0.54)', size=14, color='#555'))
parts.append(svg_sub(CX+BOX_W/2, y+85, 'Pathway-level: generalizable (C-index 0.659) | 5-FU/Capecitabine P<0.0001', size=14, color='#333'))

# ═══ Sidebar: Sample Flow Summary (compact, further right) ═══
sx = CX + BOX_W + 120
sw_box = 300
sh_box = BOX_H * 2 + GAP * 1 + 80
sidebar_top = 250
parts.append(svg_box(sx, sidebar_top, sw_box, sh_box, '#FAFAFA', '#BDBDBD', sw=1.5, rx=12))
parts.append(svg_text(sx + sw_box//2, sidebar_top + 30, 'Sample Flow Summary', size=16, weight='bold'))
# Separator
parts.append(f'<line x1="{sx+15}" y1="{sidebar_top+42}" x2="{sx+sw_box-15}" y2="{sidebar_top+42}" stroke="#ddd" stroke-width="1"/>')

entries = [
    ('Initial pool', 'n = 1,015'),
    ('After QC', 'n = 476 (47%)'),
    ('GSE39582', 'n = 239'),
    ('  Classification', 'n = 164'),
    ('  Survival', 'n = 227'),
    ('External validation', 'n = 311'),
]
for i, (label, value) in enumerate(entries):
    ey = sidebar_top + 65 + i * 35
    parts.append(svg_text(sx + 15, ey, label, size=13, weight='normal', color='#666', anchor='start'))
    parts.append(svg_text(sx + sw_box - 15, ey, value, size=13, weight='bold', color='#333', anchor='end'))

# Dashed connector from main content to sidebar
parts.append(f'<line x1="{CX+BOX_W+5}" y1="{(180+550)/2}" x2="{sx}" y2="{sidebar_top+sh_box/2}" stroke="#BDBDBD" stroke-width="1.5" stroke-dasharray="8,6"/>')

# ── Close SVG ──
parts.append('</svg>')

os.makedirs(OUT, exist_ok=True)
svg_path = os.path.join(OUT, 'Fig1_pipeline.svg')
with open(svg_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(parts))
print(f'SVG: {svg_path} ({os.path.getsize(svg_path)} bytes)')

# Copy to hd/unified
for d in ['hd', 'unified']:
    dp = os.path.join(r'/path/to/xelox_project\results\figures_png', d, 'Fig1_pipeline.svg')
    os.makedirs(os.path.dirname(dp), exist_ok=True)
    with open(dp, 'w', encoding='utf-8') as f:
        f.write('\n'.join(parts))
    print(f'  -> {d}')
print('Done')
