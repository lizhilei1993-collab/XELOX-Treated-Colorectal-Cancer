"""
Diagnostic 2: Check specific paragraphs around Table 2 area.
"""
from docx import Document

SRC = "/path/to/xelox_project/reports/manuscript_draft_v5.8_submission.docx"
doc = Document(SRC)

# Check Table 2 area (paras 73-80)
print("=== Paras 73-80 (Table 2 area) ===")
for i in range(73, min(81, len(doc.paragraphs))):
    t = doc.paragraphs[i].text.strip()
    print(f"  Para {i}: {t[:150]}")

# Check what's at the end of block 1 and start of block 3
print("\n=== Paras 126-133 (block 1->2 transition) ===")
for i in range(126, min(134, len(doc.paragraphs))):
    t = doc.paragraphs[i].text.strip()
    print(f"  Para {i}: {t[:150]}")

# Check block 2->3 transition
print("\n=== Paras 465-475 (block 2->3 transition) ===")
for i in range(465, min(476, len(doc.paragraphs))):
    t = doc.paragraphs[i].text.strip()
    print(f"  Para {i}: {t[:150]}")

# Check last few paragraphs
print("\n=== Last 5 paragraphs ===")
for i in range(max(0, len(doc.paragraphs)-5), len(doc.paragraphs)):
    t = doc.paragraphs[i].text.strip()
    print(f"  Para {i}: {t[:150]}")
