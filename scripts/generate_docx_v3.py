"""Generate DOCX from manuscript_draft_v2.md with embedded figures."""
import os, re
from docx import Document
from docx.shared import Inches, Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH

PROJECT = r"/path/to/xelox_project"
MD_PATH = os.path.join(PROJECT, "reports", "manuscript_draft_v4_nature.md")
PNG_DIR = os.path.join(PROJECT, "results", "figures_png")
OUT_PATH = os.path.join(PROJECT, "reports", "manuscript_draft_v4_nature.docx")

# Read markdown
with open(MD_PATH, 'r', encoding='utf-8') as f:
    md_text = f.read()

# Create DOCX
doc = Document()

# Page setup - A4
for section in doc.sections:
    section.page_width = Cm(21.0)
    section.page_height = Cm(29.7)
    section.top_margin = Cm(2.54)
    section.bottom_margin = Cm(2.54)
    section.left_margin = Cm(2.54)
    section.right_margin = Cm(2.54)

# Style setup
style = doc.styles['Normal']
font = style.font
font.name = 'Times New Roman'
font.size = Pt(11)
style.paragraph_format.line_spacing = 1.5
style.paragraph_format.space_after = Pt(6)

def add_heading_text(text, level):
    """Add heading with correct formatting."""
    if level == 1:
        h = doc.add_heading(text, level=1)
    elif level == 2:
        h = doc.add_heading(text, level=2)
    elif level == 3:
        h = doc.add_heading(text, level=3)
    elif level == 4:
        h = doc.add_heading(text, level=4)
    else:
        p = doc.add_paragraph()
        run = p.add_run(text)
        run.bold = True
        return p
    
    # Set font
    for run in h.runs:
        run.font.name = 'Times New Roman'
        run.font.color.rgb = RGBColor(0, 0, 0)
    return h

def add_normal_para(text):
    """Add a normal paragraph."""
    p = doc.add_paragraph()
    run = p.add_run(text)
    run.font.name = 'Times New Roman'
    run.font.size = Pt(11)
    p.paragraph_format.line_spacing = 1.5
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    return p

def add_italic_para(text):
    """Add italic paragraph (e.g. abstract labels)."""
    p = doc.add_paragraph()
    run = p.add_run(text)
    run.font.name = 'Times New Roman'
    run.font.size = Pt(11)
    run.bold = True
    run.italic = True
    return p

def add_mixed_para(parts):
    """Add paragraph with mixed bold/italic/normal parts.
    parts: list of (text, bold, italic) tuples"""
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 1.5
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    for text, bold, italic in parts:
        run = p.add_run(text)
        run.font.name = 'Times New Roman'
        run.font.size = Pt(11)
        run.bold = bold
        run.italic = italic
    return p

def try_insert_figure(fig_name, caption, width_inches=5.5):
    """Insert figure if PNG exists."""
    png_path = os.path.join(PNG_DIR, fig_name)
    if os.path.exists(png_path):
        doc.add_picture(png_path, width=Inches(width_inches))
        # Caption
        cap = doc.add_paragraph()
        run = cap.add_run(caption)
        run.font.name = 'Times New Roman'
        run.font.size = Pt(9)
        run.italic = True
        cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
        return True
    return False

# ============================================================
# Parse and build DOCX
# ============================================================

# Track: figure captions to insert
lines = md_text.split('\n')
i = 0
title_done = False
in_table = False
table_rows = []
ref_mode = False

while i < len(lines):
    line = lines[i].strip()
    
    # Title
    if line.startswith('# ') and not title_done:
        title_done = True
        title_text = line[2:].strip()
        h = doc.add_heading(title_text, level=1)
        for run in h.runs:
            run.font.name = 'Times New Roman'
            run.font.color.rgb = RGBColor(0, 0, 0)
        i += 1
        continue
    
    # Running title
    if line.startswith('**Running title'):
        p = doc.add_paragraph()
        run = p.add_run(line.replace('**', ''))
        run.font.name = 'Times New Roman'
        run.font.size = Pt(10)
        run.italic = True
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        i += 1
        continue
    
    # Authors/Affiliation
    if line.startswith('**Authors**') or line.startswith('**Affiliation**') or line.startswith('**Target Journal**'):
        p = doc.add_paragraph()
        run = p.add_run(line.replace('**', ''))
        run.font.name = 'Times New Roman'
        run.font.size = Pt(11)
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        i += 1
        continue
    
    # Separator
    if line == '---':
        i += 1
        continue
    
    # Headings
    if line.startswith('#### '):
        add_heading_text(line[5:], 4)
        i += 1
        continue
    if line.startswith('### '):
        add_heading_text(line[4:], 3)
        i += 1
        continue
    if line.startswith('## '):
        add_heading_text(line[3:], 2)
        i += 1
        continue
    
    # Bold italic labels (Abstract Background/Methods/Results/Conclusions)
    bold_italic_labels = ['**Background**', '**Methods**', '**Results**', '**Conclusions**']
    matched_label = False
    for lbl in bold_italic_labels:
        if line.startswith(lbl):
            rest = line[len(lbl):].strip()
            add_mixed_para([
                (lbl.replace('**',''), True, True),
                (' ' + rest, False, False)
            ])
            matched_label = True
            break
    if matched_label:
        i += 1
        continue
    
    # Tables
    if line.startswith('|') and line.endswith('|') and '|' in line[1:-1]:
        # Check if this is a table (has pipe separators)
        in_table = True
        table_rows.append(line)
        i += 1
        # Check if next line is a separator
        continue
    
    if in_table:
        # Check if table ended
        if not line or not line.startswith('|'):
            # Process the collected table
            if len(table_rows) > 1:
                # Parse table
                headers = [h.strip() for h in table_rows[0].split('|')[1:-1]]
                # Skip separator row (|---|)
                data_rows = []
                for row in table_rows[2:]:
                    cells = [c.strip() for c in row.split('|')[1:-1]]
                    if cells:
                        data_rows.append(cells)
                
                if data_rows:
                    ncols = len(headers)
                    table = doc.add_table(rows=len(data_rows)+1, cols=ncols)
                    table.style = 'Table Grid'
                    
                    # Header
                    for j, h in enumerate(headers):
                        cell = table.rows[0].cells[j]
                        cell.text = h
                        for p in cell.paragraphs:
                            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                            for run in p.runs:
                                run.bold = True
                                run.font.name = 'Times New Roman'
                                run.font.size = Pt(9)
                    
                    # Data
                    for r_idx, row_data in enumerate(data_rows):
                        for c_idx, cell_text in enumerate(row_data):
                            if c_idx < ncols:
                                cell = table.rows[r_idx+1].cells[c_idx]
                                cell.text = cell_text
                                for p in cell.paragraphs:
                                    for run in p.runs:
                                        run.font.name = 'Times New Roman'
                                        run.font.size = Pt(9)
                
                doc.add_paragraph()  # spacing after table
            
            in_table = False
            table_rows = []
            # Don't increment i - process this line normally
            continue
        else:
            table_rows.append(line)
            i += 1
            continue
    
    # Bold parts within text (for emphasis)
    if '**' in line:
        parts = re.split(r'(\*\*[^*]+\*\*)', line)
        p = doc.add_paragraph()
        p.paragraph_format.line_spacing = 1.5
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        for part in parts:
            if part.startswith('**') and part.endswith('**'):
                run = p.add_run(part[2:-2])
                run.bold = True
            else:
                run = p.add_run(part)
            run.font.name = 'Times New Roman'
            run.font.size = Pt(11)
        i += 1
        continue
    
    # Empty lines - skip but add spacing for figure/table separators
    if not line:
        i += 1
        continue
    
    # Regular text
    if line.startswith('*(e.g.') or line.startswith('(*'):
        add_italic_para(line)
    else:
        add_normal_para(line)
    
    i += 1

# ============================================================
# Add figures at the end of the document
# ============================================================
# Add figures
figure_map = {
    'Fig1_pipeline.png': 'Figure 1. Study flowchart detailing the multi-phase analytical pipeline.',
    'Fig2_SHAP_beeswarm.png': 'Figure 2. SHAP summary plot showing the 10-gene fingerprint.',
    'Fig3A_pathway_volcano.png': 'Figure 3A. Volcano plot of pathway-level differential activity.',
    'Fig3B_pathway_heatmap.png': 'Figure 3B. Consensus pathway heatmap.',
    'Fig3C_pathway_boxplot.png': 'Figure 3C. Boxplots of top consensus pathways.',
    'Fig4A_PRS_ROC.png': 'Figure 4A. Bootstrap selection frequency of 3-pathway LASSO-PRS.',
    'Fig4B_prs_km_curve.png': 'Figure 4B. Kaplan–Meier curves stratified by median risk score.',
    'Fig4C_cox_forest.png': 'Figure 4C. Cox regression forest plot.',
    'Fig4D_nomogram.png': 'Figure 4D. Prognostic nomogram for 12-/36-/60-month RFS.',
    'Fig4E_calibration_12mo.png': 'Figure 4E. 12-month calibration plot.',
    'Fig4E_calibration_36mo.png': 'Figure 4F. 36-month calibration plot.',
    'Fig4F_calibration_60mo.png': 'Figure 4G. 60-month calibration plot.',
    'Fig5_BP_dotplot.png': 'Figure 5A. GO BP ORA dotplot.',
    'Fig5_KEGG_barplot.png': 'Figure 5B. KEGG enrichment barplot.',
    'Fig6A_ablation_forest.png': 'Figure 6A. Gene-level ablation forest plot.',
    'Fig6B_method_comparison.png': 'Figure 6B. Feature type × modeling method comparison.',
}

for fname, caption in figure_map.items():
    fpath = os.path.join(PNG_DIR, fname)
    if os.path.exists(fpath):
        try:
            doc.add_picture(fpath, width=Inches(5.5))
            cap = doc.add_paragraph()
            run = cap.add_run(caption)
            run.font.name = 'Times New Roman'
            run.font.size = Pt(9)
            run.italic = True
            cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
            doc.add_paragraph()
        except:
            doc.add_paragraph(f"[Figure: {caption}]")
    else:
        print(f"  WARNING: {fname} not found")

# Save
doc.save(OUT_PATH)
print(f"\nDOCX saved: {OUT_PATH}")
print(f"Size: {os.path.getsize(OUT_PATH)/1024:.0f} KB")
