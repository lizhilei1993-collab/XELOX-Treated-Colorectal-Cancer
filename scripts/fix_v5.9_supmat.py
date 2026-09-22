"""
Clean up v5.9 supplementary materials:
1. Delete duplicate tables 7-11 (identical to tables 2-6)
2. Delete file organization / document summary / redundant GEO refs
3. Keep block 1 (inline sup tables) + block 2 (full SupMat)
"""
from lxml import etree
import subprocess
import sys

UNPACK_DIR = "/path/to/xelox_project/reports/unpacked_v5.9"
XML_PATH = f"{UNPACK_DIR}/word/document.xml"
DST_DOCX = "/path/to/xelox_project/reports/manuscript_draft_v5.9_submission.docx"
SRC_DOCX = "/path/to/xelox_project/reports/manuscript_draft_v5.9_submission.docx"

ns = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'

tree = etree.parse(XML_PATH)
root = tree.getroot()
body = root.find(f'.//{{{ns}}}body')

# Get all paragraphs and tables
all_elements = list(body)
print(f"Total body elements: {len(all_elements)}")

# Strategy: remove paragraphs and tables from block 2 that contain
# the duplicate/redundant content (File Organization, Document Summary,
# redundant tables, redundant GEO refs)
# Identified from docx audit: Para 423-467 in python-docx indexing

# But XML element ordering may differ from python-docx paragraph index.
# So I'll find elements by content.

# Markers for content to REMOVE (in block 2, after the full SupMat content)
REMOVE_MARKERS = [
    "File Organization",
    "Resource Summary",
    "Directory Structure",
    "PNG figures for manuscript",
    "PDF figures (high-resolution)",
    "CSV tables:",
    "GEO cohorts analyzed:",
    "Total cohort N:",
    "results/",
    "figures/",
    "figures_png/",
    "figures_prs/",
    "tables/",
    "scripts/",
    "reports/",
    "References (GEO Datasets)",
    "Document Summary",
]

# For tables, find tables with content "Supplementary Table S27" etc.
# in the SECOND occurrence (duplicate)
# We already know Tables 2-6 are the first set, Tables 7-11 are duplicates
# Need to find and remove all elements belonging to duplicate tables

# Find the position of "File Organization" paragraph - this marks start of junk
elements_to_remove = []
found_file_org = False
tables_removed = 0
paras_removed = 0

for elem in all_elements:
    # Check if it's a paragraph
    tag = etree.QName(elem).localname
    
    if tag == 'p':
        # Get text content
        texts = []
        for t in elem.iter(f'{{{ns}}}t'):
            if t.text:
                texts.append(t.text)
        full_text = ''.join(texts).strip()
        
        # Check for content to remove
        if full_text == "File Organization":
            found_file_org = True
        
        if found_file_org:
            # Remove everything from File Organization onwards in block 2
            # But only if we're in the second Supplementary Materials block
            # (after all the S1-S35 tables and S1-S31 figures)
            elements_to_remove.append(elem)
            paras_removed += 1
    
    elif tag == 'tbl':
        if found_file_org:
            elements_to_remove.append(elem)
            tables_removed += 1

print(f"Found elements to remove: {paras_removed} paragraphs + {tables_removed} tables")

# Remove marked elements
for elem in elements_to_remove:
    parent = elem.getparent()
    if parent is not None:
        parent.remove(elem)

print(f"Removed {len(elements_to_remove)} elements total")

# Now check for duplicate "Supplementary Tables S1-S16" text at the end
# and remove it
for elem in body:
    tag = etree.QName(elem).localname
    if tag == 'p':
        texts = []
        for t in elem.iter(f'{{{ns}}}t'):
            if t.text:
                texts.append(t.text)
        full_text = ''.join(texts).strip()
        if "Supplementary Tables: 35 (S1" in full_text:
            parent = elem.getparent()
            if parent is not None:
                parent.remove(elem)
                print(f"Removed stray 'Supplementary Tables: 35' paragraph")

# Also check for the "Supplementary Tables S1-S16 and Figures S1-S12 are provided 
# in a separate file" paragraph - keep this as it's the overview note
# Actually, let me also check for "Supplementary Tables" headers that appear
# in wrong places

# Save and repack
tree.write(XML_PATH, xml_declaration=True, encoding='UTF-8', standalone=True)

# Count remaining elements
tree2 = etree.parse(XML_PATH)
root2 = tree2.getroot()
body2 = root2.find(f'.//{{{ns}}}body')
remaining = len(list(body2))
print(f"Remaining body elements: {remaining}")

sys.path.insert(0, r"/home/user\.workbuddy\plugins\marketplaces\codebuddy-plugins-official\plugins\docx\scripts\office")
cmd = [
    "C:/ProgramData/WorkBuddy/chromium-env/6nqqlh/.workbuddy/binaries/python/envs/default/Scripts/python.exe",
    r"/home/user\.workbuddy\plugins\marketplaces\codebuddy-plugins-official\plugins\docx\scripts\office\pack.py",
    UNPACK_DIR,
    DST_DOCX,
    "--original", SRC_DOCX,
]
result = subprocess.run(cmd, capture_output=True, text=True)
print(f"Repack: {'OK' if result.returncode == 0 else 'ERROR: ' + result.stderr}")
