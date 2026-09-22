"""
Remove 10 dead references from Reference list in v5.9.
Refs to remove: [14, 22, 24, 25, 26, 29, 30, 31, 32, 33]
These are in the reference list but never cited in the manuscript text.
"""
from lxml import etree
import subprocess
import sys
import re

UNPACK_DIR = "/path/to/xelox_project/reports/unpacked_v5.9"
XML_PATH = f"{UNPACK_DIR}/word/document.xml"
DST_DOCX = "/path/to/xelox_project/reports/manuscript_draft_v5.9_submission.docx"
SRC_DOCX = "/path/to/xelox_project/reports/manuscript_draft_v5.9_submission.docx"

ns = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'

tree = etree.parse(XML_PATH)
root = tree.getroot()
body = root.find(f'.//{{{ns}}}body')

# Dead reference numbers to remove
DEAD_REFS = {14, 22, 24, 25, 26, 29, 30, 31, 32, 33}

# Find all paragraphs that look like reference entries
# References are blockquote paragraphs starting with "1. " etc.
ref_paras = []  # (ref_number, element)

for elem in body:
    tag = etree.QName(elem).localname
    if tag != 'p':
        continue
    
    # Get first text run
    first_text = ''
    for t in elem.iter(f'{{{ns}}}t'):
        if t.text:
            first_text = t.text.strip()
            break
    
    # Check if it starts with a number followed by ". " (e.g., "1. Sung H...")
    m = re.match(r'^(\d+)\.\s', first_text)
    if m:
        ref_num = int(m.group(1))
        ref_paras.append((ref_num, elem))

print(f"Found {len(ref_paras)} reference paragraphs")
print(f"Ref numbers: {[r[0] for r in ref_paras]}")

# Remove dead references
removed = 0
for ref_num, elem in ref_paras:
    if ref_num in DEAD_REFS:
        parent = elem.getparent()
        if parent is not None:
            parent.remove(elem)
            print(f"  Removed ref [{ref_num}]")
            removed += 1

print(f"Removed {removed} dead references")

# Save
tree.write(XML_PATH, xml_declaration=True, encoding='UTF-8', standalone=True)

# Repack
sys.path.insert(0, r"/home/user\.workbuddy\plugins\marketplaces\codebuddy-plugins-official\plugins\docx\scripts\office")
cmd = [
    "C:/ProgramData/WorkBuddy/chromium-env/6nqqlh/.workbuddy/binaries/python/envs/default/Scripts/python.exe",
    r"/home/user\.workbuddy\plugins\marketplaces\codebuddy-plugins-official\plugins\docx\scripts\office\pack.py",
    UNPACK_DIR,
    DST_DOCX,
    "--original", SRC_DOCX,
    "--validate", "false",
]
result = subprocess.run(cmd, capture_output=True, text=True)
print(f"Repack: {'OK' if result.returncode == 0 else 'ERROR: ' + result.stderr}")
