"""
Fix v5.8 → v5.9 via lxml on unpacked document.xml.
Finds paragraphs by text content markers and removes them.
Then repacks to .docx.
"""
from lxml import etree
import shutil
import os
import subprocess
import sys

UNPACK_DIR = "/path/to/xelox_project/reports/unpacked_v5.8"
XML_PATH = f"{UNPACK_DIR}/word/document.xml"
DST_DOCX = "/path/to/xelox_project/reports/manuscript_draft_v5.9_submission.docx"
SRC_DOCX = "/path/to/xelox_project/reports/manuscript_draft_v5.8_submission.docx"

# Parse XML
nsmap = {
    'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
    'w14': 'http://schemas.microsoft.com/office/word/2010/wordml',
}

tree = etree.parse(XML_PATH)
root = tree.getroot()

# Find body
body = root.find('.//w:body', nsmap)
if body is None:
    print("ERROR: Could not find w:body")
    sys.exit(1)

# Get all paragraphs
all_paras = body.findall('.//w:p', nsmap)
print(f"Total paragraphs: {len(all_paras)}")

# Track which paragraphs to remove
paras_to_remove = []
marker_texts = {}

# Scan all paragraphs and collect text
for i, p in enumerate(all_paras):
    # Get all text runs in this paragraph
    texts = []
    for t in p.iter('{http://schemas.openxmlformats.org/wordprocessingml/2006/main}t'):
        if t.text:
            texts.append(t.text)
    full_text = ''.join(texts).strip()
    
    if "comprehensive overfitting diagnostics summary" in full_text:
        if "overfitting" not in marker_texts:
            marker_texts["overfitting"] = [i]
        else:
            marker_texts["overfitting"].append(i)
    
    if full_text.startswith("Note: HR = hazard ratio; CI = confidence interval."):
        if "table2_note" not in marker_texts:
            marker_texts["table2_note"] = [i]
        else:
            marker_texts["table2_note"].append(i)
    
    if "Supplementary Tables: S1" in full_text or full_text == "Supplementary Tables: S1–S35 (35 tables)":
        if "supp_tables" not in marker_texts:
            marker_texts["supp_tables"] = [i]
        else:
            marker_texts["supp_tables"].append(i)
    
    if full_text == "Supplementary Materials":
        if "supp_header" not in marker_texts:
            marker_texts["supp_header"] = [i]
        else:
            marker_texts["supp_header"].append(i)

print("\nFound markers:")
for key, indices in marker_texts.items():
    print(f"  {key}: {indices}")

# --- Determine which paragraphs to remove ---

# Fix 1: Remove 2nd overfitting paragraph (keep 1st)
if len(marker_texts.get("overfitting", [])) >= 2:
    remove_idx = marker_texts["overfitting"][1]
    paras_to_remove.append(remove_idx)
    print(f"Fix 1: Remove duplicate Discussion paragraph #{remove_idx}")

# Fix 2: Remove 1st Table 2 note (keep 2nd)
if len(marker_texts.get("table2_note", [])) >= 2:
    remove_idx = marker_texts["table2_note"][0]
    paras_to_remove.append(remove_idx)
    print(f"Fix 2: Remove conflicting Table 2 note #{remove_idx}")

# Fix 3: Remove 3rd Supplementary Materials block
# Find the 3rd "Supplementary Materials" header and remove everything from there
if len(marker_texts.get("supp_header", [])) >= 3:
    start_idx = marker_texts["supp_header"][2]
    print(f"Fix 3: Remove 3rd Supplementary Materials block starting at para #{start_idx}")
    for i in range(start_idx, len(all_paras)):
        paras_to_remove.append(i)
elif len(marker_texts.get("supp_tables", [])) >= 3:
    # The "Supplementary Tables: S1–S35" header appears in each block
    start_idx = marker_texts["supp_tables"][2]
    # Need to go back to find the paragraph with "Supplementary Materials" heading
    actual_start = start_idx - 2  # usual structure: header + title + content
    print(f"Fix 3: Remove 3rd Supplementary Materials block starting near para #{actual_start}")
    for i in range(actual_start, len(all_paras)):
        paras_to_remove.append(i)
else:
    # Fallback: use the last supplementary tables marker
    print(f"WARNING: Could not identify 3rd Supplementary Materials block")
    print(f"  supp_header at: {marker_texts.get('supp_header', [])}")
    print(f"  supp_tables at: {marker_texts.get('supp_tables', [])}")

# Deduplicate and sort in reverse
paras_to_remove = sorted(set(paras_to_remove), reverse=True)
print(f"\nRemoving {len(paras_to_remove)} paragraphs: {paras_to_remove[:10]}{'...' if len(paras_to_remove) > 10 else ''}")

# Remove paragraphs from body
for idx in paras_to_remove:
    if idx < len(all_paras):
        p = all_paras[idx]
        parent = p.getparent()
        if parent is not None:
            # Remove preceding whitespace text nodes
            prev = p.getprevious()
            if prev is not None and prev.tail and prev.tail.strip() == '':
                prev.tail = ''
            parent.remove(p)

# Save modified XML
tree.write(XML_PATH, xml_declaration=True, encoding='UTF-8', standalone=True)
print(f"\nXML saved. Remaining paragraphs: {len(body.findall('.//w:p', nsmap))}")

# Repack
print("\nRepacking to docx...")
sys.path.insert(0, r"/home/user\.workbuddy\plugins\marketplaces\codebuddy-plugins-official\plugins\docx\scripts\office")

# Run pack.py
cmd = [
    "/path/to/cache/chromium-env/6nqqlh/.workbuddy/binaries/python/envs/default/Scripts/python.exe",
    r"/home/user\.workbuddy\plugins\marketplaces\codebuddy-plugins-official\plugins\docx\scripts\office\pack.py",
    UNPACK_DIR,
    DST_DOCX,
    "--original", SRC_DOCX,
]
result = subprocess.run(cmd, capture_output=True, text=True)
if result.returncode == 0:
    print(f"Successfully created: {DST_DOCX}")
else:
    print(f"Pack error: {result.stderr}")
