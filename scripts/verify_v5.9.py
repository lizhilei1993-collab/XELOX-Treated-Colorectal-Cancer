"""
Verify v5.9 fixes.
"""
from docx import Document

DST = "/path/to/xelox_project/reports/manuscript_draft_v5.9_submission.docx"
doc = Document(DST)

issues_fixed = {
    "discussion_dup": 0,
    "table2_note": 0,
    "supp_materials": 0,
}

for i, p in enumerate(doc.paragraphs):
    t = p.text.strip()
    
    if "comprehensive overfitting diagnostics summary" in t:
        issues_fixed["discussion_dup"] += 1
    
    if t.startswith("Note: HR = hazard ratio; CI = confidence interval."):
        issues_fixed["table2_note"] += 1
        print(f"  Table 2 note: {t[:100]}...")
    
    if t == "Supplementary Materials":
        issues_fixed["supp_materials"] += 1

print(f"Verification results:")
print(f"  'overfitting diagnostics' occurrences: {issues_fixed['discussion_dup']} (expected: 1)")
print(f"  'Table 2 Note' occurrences: {issues_fixed['table2_note']} (expected: 1)")
print(f"  'Supplementary Materials' blocks: {issues_fixed['supp_materials']} (expected: 2)")

# Quick check on the remaining Supplementary Materials content
for i, p in enumerate(doc.paragraphs):
    if t == "Supplementary Materials":
        # Look at next paragraph to see which block it is
        for j in range(1, 5):
            if i+j < len(doc.paragraphs):
                next_t = doc.paragraphs[i+j].text.strip()
                if next_t:
                    print(f"  Block at para {i}: next content = '{next_t[:80]}'")
                    break

print(f"\nTotal paragraphs in v5.9: {len(doc.paragraphs)}")
