"""
Fix 1: Round excessive decimals in Table S30 (S27 in block1).
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
changes = 0

for t_elem in root.iter(f'{{{ns}}}t'):
    if t_elem.text is None:
        continue
    val_str = t_elem.text.strip()
    
    # Match numbers with excessive decimal places (10+ digits after decimal, or >10 chars with decimal)
    if re.match(r'^-?\d+\.\d{8,}$', val_str):
        try:
            val = float(val_str)
            # Don't round p-values (already handled with e- notation)
            rounded = f"{val:.4f}"
            if rounded != val_str:
                t_elem.text = rounded
                print(f"  {val_str} → {rounded}")
                changes += 1
        except ValueError:
            pass

print(f"Decimal rounding changes: {changes}")

tree.write(XML_PATH, xml_declaration=True, encoding='UTF-8', standalone=True)

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
print(f"Repack: {'OK' if result.returncode == 0 else 'ERROR'}")
