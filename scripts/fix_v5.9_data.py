"""
Fix v5.9 remaining issues based on source data verification.

Fixes:
1. B=25 → B=200/1000 (actual code used B=200 gene, B=1000 pathway)
2. optimism 0.095 → 0.0945 (actual computed value)
3. Table S33 Bootstrap-corrected C-index 0.692 → 0.659
   (0.692 is bootstrap mean; 0.659 is optimism-corrected)
"""
from lxml import etree
import subprocess
import os
import sys

UNPACK_DIR = "/path/to/xelox_project/reports/unpacked_v5.9"
XML_PATH = f"{UNPACK_DIR}/word/document.xml"
DST_DOCX = "/path/to/xelox_project/reports/manuscript_draft_v5.9_submission.docx"
SRC_DOCX = "/path/to/xelox_project/reports/manuscript_draft_v5.9_submission.docx"

# Parse XML (preserve text nodes)
parser = etree.XMLParser(strip_cdata=False)
tree = etree.parse(XML_PATH, parser)
root = tree.getroot()
ns = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'

changes = 0

# Find all <w:t> elements
for t_elem in root.iter(f'{{{ns}}}t'):
    if t_elem.text is None:
        continue
    
    original = t_elem.text
    
    # Fix 2: optimism precision
    # In Abstract and Discussion, change 0.095 → 0.0945
    if "optimism: 0.011 vs. 0.095" in original:
        t_elem.text = original.replace("0.095", "0.0945")
        print(f"Fix 2a (Abstract): 0.095 → 0.0945")
        changes += 1
    
    # "optimism = 0.095" in fair comparison section
    elif "optimism = 0.095" in original and "apparent" not in original.lower():
        t_elem.text = original.replace("optimism = 0.095", "optimism = 0.0945")
        print(f"Fix 2b (Results): optimism = 0.095 → 0.0945")
        changes += 1

    # Fix 1: B=25 → actual values
    if "Bootstrap optimism correction (B = 25 resamples each)" in original:
        t_elem.text = original.replace(
            "Bootstrap optimism correction (B = 25 resamples each)",
            "Bootstrap optimism correction (B = 200 for gene level; B = 1,000 for pathway level)"
        )
        print(f"Fix 1: B=25 → B=200/1000")
        changes += 1
    
    # Fix 3: Table S33 Bootstrap-corrected C-index 0.692 → 0.659
    if t_elem.text.strip() == "0.692":
        # Check context: is this inside a table cell?
        parent_p = t_elem.getparent()
        if parent_p is not None:
            # Check if any ancestor is a table cell
            has_table = False
            ancestor = parent_p.getparent()
            for _ in range(5):
                if ancestor is None:
                    break
                if ancestor.tag == f'{{{ns}}}tc':
                    has_table = True
                    break
                ancestor = ancestor.getparent()
            
            # Also check if this is the S33 table specifically
            # by looking for nearby text "Bootstrap-corrected C-index"
            if has_table:
                # Check nearby rows for the label
                grandparent = parent_p.getparent()
                if grandparent is not None:
                    # Search backwards for "Bootstrap-corrected C-index"
                    prev_texts = []
                    for prev in grandparent.iter(f'{{{ns}}}t'):
                        if prev.text:
                            prev_texts.append(prev.text)
                    full_context = ' '.join(prev_texts)
                    if "Bootstrap-corrected C-index" in full_context:
                        t_elem.text = "0.659"
                        print(f"Fix 3: Table S33 Bootstrap-corrected C-index 0.692 → 0.659")
                        changes += 1

# Fix 4: Also fix the second copy of Table S33 (in supplementary materials block 2)
for t_elem in root.iter(f'{{{ns}}}t'):
    if t_elem.text is None:
        continue
    if "Bootstrap-corrected C-index" in t_elem.text:
        # Find the corresponding value cell in same row
        # The value 0.692 should be in the next cell
        pass

# Actually, let me do a simpler approach for Table S33:
# Find ALL "0.692" in table cells that have "Bootstrap-corrected" nearby
for t_elem in root.iter(f'{{{ns}}}t'):
    if t_elem.text and t_elem.text.strip() == "0.692":
        # Walk up to find table context
        cur = t_elem.getparent()
        table_found = False
        for _ in range(6):
            if cur is None:
                break
            if cur.tag == f'{{{ns}}}tbl':
                # Found a table. Search for "Bootstrap-corrected C-index" in this table
                table_text = ''
                for tt in cur.iter(f'{{{ns}}}t'):
                    if tt.text:
                        table_text += tt.text + ' '
                if "Bootstrap-corrected C-index" in table_text:
                    if "0.692" not in table_text.replace("0.692", "CHECKED", 1):
                        # This is the value cell
                        t_elem.text = "0.659"
                        print(f"Fix 3b: Table S33 value 0.692 → 0.659 (in table containing 'Bootstrap-corrected C-index')")
                        changes += 1
                table_found = True
                break
            cur = cur.getparent()

# Save modified XML
tree.write(XML_PATH, xml_declaration=True, encoding='UTF-8', standalone=True)
print(f"\nTotal changes: {changes}")

# Repack
sys.path.insert(0, r"/home/user\.workbuddy\plugins\marketplaces\codebuddy-plugins-official\plugins\docx\scripts\office")
cmd = [
    "C:/ProgramData/WorkBuddy/chromium-env/6nqqlh/.workbuddy/binaries/python/envs/default/Scripts/python.exe",
    r"/home/user\.workbuddy\plugins\marketplaces\codebuddy-plugins-official\plugins\docx\scripts\office\pack.py",
    UNPACK_DIR,
    DST_DOCX,
    "--original", SRC_DOCX,
]
result = subprocess.run(cmd, capture_output=True, text=True)
if result.returncode == 0:
    print(f"Successfully repacked: {DST_DOCX}")
else:
    print(f"Pack error: {result.stderr}")
