"""
Fix v5.8 → v5.9: Apply verified corrections to manuscript.

Verified issues to fix:
1. DELETE Para 58: Duplicate Discussion paragraph (== Para 56)
2. DELETE Para 75: Conflicting Table 2 note v1 (keep Para 76+77)
3. DELETE Paras 470-792: Duplicate Supplementary Materials block 3 (V5.1, outdated)
4. KEEP everything else (Table S33 data flagged for manual review)
"""
from docx import Document
from docx.shared import Pt
import shutil
import os

SRC = "/path/to/xelox_project/reports/manuscript_draft_v5.8_submission.docx"
DST = "/path/to/xelox_project/reports/manuscript_draft_v5.9_submission.docx"

shutil.copy2(SRC, DST)
doc = Document(DST)

# Collect paragraphs to delete (in reverse order to avoid index shifting)
to_delete = []

# Issue 1: Duplicate Discussion paragraph (Para 58)
to_delete.append(58)

# Issue 2: Conflicting Table 2 note v1 (Para 75)
to_delete.append(75)

# Issue 3: Duplicate Supplementary Materials block 3 (Paras 470-792)
for i in range(470, len(doc.paragraphs)):
    to_delete.append(i)

# Sort in descending order and delete
to_delete = sorted(set(to_delete), reverse=True)

print(f"Deleting {len(to_delete)} paragraphs...")
deleted_count = 0
for idx in to_delete:
    if idx < len(doc.paragraphs):
        p = doc.paragraphs[idx]
        p._element.getparent().remove(p._element)
        deleted_count += 1
    else:
        print(f"  WARNING: Para {idx} out of range (len={len(doc.paragraphs)})")

print(f"Deleted: {deleted_count} paragraphs")

# Verify the fixes
doc2 = Document(DST)

# Check for Discussion duplicates
discussion_count = 0
for p in doc2.paragraphs:
    if "A comprehensive overfitting diagnostics summary" in p.text:
        discussion_count += 1
print(f"Discussion 'overfitting diagnostics' occurrences: {discussion_count} (expected: 1)")

# Check for Table 2 notes
table2_count = 0
for p in doc2.paragraphs:
    if p.text.strip().startswith("Note: HR = hazard ratio; CI = confidence interval."):
        table2_count += 1
        print(f"  Found: {p.text.strip()[:120]}")
print(f"Table 2 Note occurrences: {table2_count} (expected: 1)")

# Check for Supplementary Materials blocks
supp_count = 0
for p in doc2.paragraphs:
    if p.text.strip() == "Supplementary Materials":
        supp_count += 1
print(f"Supplementary Materials blocks: {supp_count} (expected: 2 - overview + full)")

# Check total paragraphs
print(f"Total paragraphs after fix: {len(doc2.paragraphs)} (was 793)")

doc.save(DST)
print(f"\nSaved to: {DST}")
