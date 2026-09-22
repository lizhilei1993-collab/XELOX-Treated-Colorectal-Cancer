"""
Deep audit of v5.9 supplementary materials - every paragraph, every table.
"""
from docx import Document
from collections import Counter

SRC = "/path/to/xelox_project/reports/manuscript_draft_v5.9_submission.docx"
doc = Document(SRC)

# 1. Full paragraph inventory of supmat blocks
print("=" * 70)
print("FULL SUPPLEMENTARY MATERIALS INVENTORY")
print("=" * 70)

in_supmat = False
in_block1 = False
in_block2 = False
block1_content = []
block2_content = []

for i, p in enumerate(doc.paragraphs):
    t = p.text.strip()
    
    # Detect block 1 start (first Supplementary Materials after References)
    if t == "Supplementary Materials" and not in_supmat:
        in_supmat = True
        in_block1 = True
        continue
    
    if t == "Supplementary Materials" and in_block1:
        in_block1 = False
        in_block2 = True
        continue
    
    if in_block1 and t:
        block1_content.append((i, t[:100], 'B1'))
    
    if in_block2 and t:
        block2_content.append((i, t[:100], 'B2'))

print(f"\n--- Block 1 (Overview + Inline Tables): {len(block1_content)} items ---")
for idx, text, label in block1_content:
    print(f"  [{label}] Para {idx}: {text}")

print(f"\n--- Block 2 (Full SupMat): {len(block2_content)} items ---")

# For block 2, check for structural issues
seen_headers = Counter()
for idx, text, label in block2_content:
    # Track section headers
    if text.startswith("Supplementary") and ("Table" in text or "Figure" in text):
        seen_headers[text[:30]] += 1

# Find duplicate headers in block 2
print("\nDuplicate headers in Block 2:")
for header, count in seen_headers.items():
    if count > 1:
        print(f"  *** {count}x: {header}")

# Check block 2 structure - print first 20 and last 20 items
print("\n--- Block 2: First 30 items ---")
for idx, text, label in block2_content[:30]:
    print(f"  [{label}] Para {idx}: {text}")

print(f"\n--- Block 2: Last 30 items ---")
for idx, text, label in block2_content[-30:]:
    print(f"  [{label}] Para {idx}: {text}")

# Summary of what Block 2 contains
print("\n--- Block 2 Summary ---")
# Count unique table numbers mentioned
import re
table_nums = []
fig_nums = []
for _, text, _ in block2_content:
    m = re.search(r'Table S(\d+)', text)
    if m:
        table_nums.append(int(m.group(1)))
    m = re.search(r'Figure S(\d+)', text)
    if m:
        fig_nums.append(int(m.group(1)))

if table_nums:
    print(f"  Tables referenced: S{min(table_nums)}-S{max(table_nums)} ({len(set(table_nums))} unique)")
if fig_nums:
    print(f"  Figures referenced: S{min(fig_nums)}-S{max(fig_nums)} ({len(set(fig_nums))} unique)")

# Check table content for duplicates
print("\n--- All 7 Document Tables ---")
table_ids = []
for i, table in enumerate(doc.tables):
    if len(table.rows) > 0:
        header = ' | '.join(cell.text.strip()[:25] for cell in table.rows[0].cells)
        n_rows = len(table.rows)
        
        # Determine which block this table belongs to
        # Find the paragraph right before this table
        print(f"  Table {i} ({n_rows} rows): {header}")
        # Print first data row
        if len(table.rows) > 1:
            data_row = ' | '.join(c.text.strip()[:20] for c in table.rows[1].cells)
            print(f"    Row 1: {data_row}")
