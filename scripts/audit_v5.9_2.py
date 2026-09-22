"""
Detailed view of the duplicate region at end of block 2.
"""
from docx import Document

SRC = "/path/to/xelox_project/reports/manuscript_draft_v5.9_submission.docx"
doc = Document(SRC)

print("=== Paras 420-467 (end of block 2) ===")
for i in range(420, min(468, len(doc.paragraphs))):
    t = doc.paragraphs[i].text.strip()
    if t:
        # Show tables
        print(f"  Para {i}: {t[:120]}")
    else:
        # Check if it's a table cell
        pass

# Check what's in tables 7-11 (the duplicate ones)
print("\n=== Tables 7-11 (potential duplicates) ===")
for i in range(7, min(12, len(doc.tables))):
    table = doc.tables[i]
    first_row = ' | '.join(cell.text.strip()[:25] for cell in table.rows[0].cells)
    n_rows = len(table.rows)
    print(f"  Table {i}: {n_rows} rows — {first_row}")

# Compare Tables 2-6 (first set) with Tables 7-11 (second set)
print("\n=== Comparing first set (T2-T6) with second set (T7-T11) ===")
for t_idx in range(5):
    t1 = doc.tables[2 + t_idx]
    t2 = doc.tables[7 + t_idx]
    first1 = ' | '.join(cell.text.strip()[:20] for cell in t1.rows[0].cells)
    first2 = ' | '.join(cell.text.strip()[:20] for cell in t2.rows[0].cells)
    match = "IDENTICAL" if first1 == first2 else "DIFFERENT"
    print(f"  T{2+t_idx} vs T{7+t_idx}: {match}")
    print(f"    T{2+t_idx}: {first1}")
    print(f"    T{7+t_idx}: {first2}")
