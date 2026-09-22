"""
Full structural audit of v5.9 supplementary materials.
"""
from docx import Document

SRC = "/path/to/xelox_project/reports/manuscript_draft_v5.9_submission.docx"
doc = Document(SRC)

print(f"Total paragraphs: {len(doc.paragraphs)}")

# Scan ALL paragraphs for structural markers
supp_markers = []
table_s_markers = []
fig_s_markers = []
table_markers = []
fig_markers = []
heading_markers = []

for i, p in enumerate(doc.paragraphs):
    t = p.text.strip()
    if not t:
        continue
    
    # Supplementary Materials header
    if t in ["Supplementary Materials", "**Supplementary Materials**"]:
        supp_markers.append((i, t))
    
    # Supplementary Tables (S1-S35)
    if t.startswith("Supplementary Table") and ("S" in t or "s" in t):
        table_s_markers.append((i, t[:80]))
    
    # Supplementary Figures
    if t.startswith("Supplementary Figure") or t.startswith("Figure S"):
        fig_s_markers.append((i, t[:80]))
    
    # Main tables (Table 1, Table 2)
    if t.startswith("Table ") and "Supplementary" not in t:
        table_markers.append((i, t[:80]))
    
    # Main figures
    if t.startswith("Figure ") and "Supplementary" not in t:
        fig_markers.append((i, t[:80]))
    
    # Headings
    if p.style and p.style.name.startswith("Heading"):
        heading_markers.append((i, p.style.name, t[:80]))

print(f"\n=== Supplementary Materials headers: {len(supp_markers)} ===")
for idx, t in supp_markers:
    print(f"  Para {idx}: {t}")

print(f"\n=== Supplementary Tables: {len(table_s_markers)} ===")
seen_tables = {}
for idx, t in table_s_markers:
    # Extract table number
    key = t
    if key in seen_tables:
        print(f"  *** DUPLICATE *** Para {idx}: {t}")
    else:
        seen_tables[key] = idx
        print(f"  Para {idx}: {t}")

print(f"\n=== Supplementary Figures: {len(fig_s_markers)} ===")
seen_figs = {}
for idx, t in fig_s_markers:
    key = t
    if key in seen_figs:
        print(f"  *** DUPLICATE *** Para {idx}: {t}")
    else:
        seen_figs[key] = idx
        print(f"  Para {idx}: {t}")

print(f"\n=== Main Tables: {len(table_markers)} ===")
for idx, t in table_markers:
    print(f"  Para {idx}: {t}")

print(f"\n=== Main Figures: {len(fig_markers)} ===")
for idx, t in fig_markers:
    print(f"  Para {idx}: {t}")

print(f"\n=== Headings: {len(heading_markers)} ===")
for idx, style, t in heading_markers:
    print(f"  Para {idx} [{style}]: {t}")

# Check for tables in document tables (not just paragraphs)
print(f"\n=== Document Tables: {len(doc.tables)} ===")
for i, table in enumerate(doc.tables):
    # Get first row content as identifier
    if len(table.rows) > 0:
        first_row = ' | '.join(cell.text.strip()[:30] for cell in table.rows[0].cells)
        n_rows = len(table.rows)
        n_cols = len(table.columns)
        print(f"  Table {i}: {n_rows}x{n_cols} — {first_row[:100]}")
