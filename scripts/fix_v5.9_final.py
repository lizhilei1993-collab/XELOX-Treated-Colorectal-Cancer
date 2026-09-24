"""
Final fix: remaining B=25 and optimism=0.095 issues.
"""
from lxml import etree
import subprocess
import sys

UNPACK_DIR = "/path/to/xelox_project/reports/unpacked_v5.9"
XML_PATH = f"{UNPACK_DIR}/word/document.xml"
DST_DOCX = "/path/to/xelox_project/reports/manuscript_draft_v5.9_submission.docx"
SRC_DOCX = "/path/to/xelox_project/reports/manuscript_draft_v5.9_submission.docx"

tree = etree.parse(XML_PATH)
ns = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'
changes = 0

for t_elem in tree.iter(f'{{{ns}}}t'):
    if t_elem.text is None:
        continue
    
    # Fix: "optimism = 0.095" → "0.0945" in all contexts (Results section)
    if "optimism = 0.095" in t_elem.text:
        t_elem.text = t_elem.text.replace("optimism = 0.095", "optimism = 0.0945")
        print(f"Fixed: optimism = 0.095 → 0.0945 in: ...{t_elem.text[-80:]}")
        changes += 1
    
    # Fix: Methods B=25
    if "bootstrap optimism correction (B = 25)" in t_elem.text:
        t_elem.text = t_elem.text.replace(
            "bootstrap optimism correction (B = 25)",
            "bootstrap optimism correction (B = 200 for gene level; B = 1,000 for pathway level)"
        )
        print(f"Fixed: B=25 → B=200/1000 in Methods")
        changes += 1
    
    # Fix: Supplementary materials still has "0.095"
    if "Optimism: 0.095" in t_elem.text:
        t_elem.text = t_elem.text.replace("Optimism: 0.095", "Optimism: 0.0945")
        print(f"Fixed: Optimism: 0.095 → 0.0945 in supplementary")
        changes += 1

print(f"\nTotal changes: {changes}")

tree.write(XML_PATH, xml_declaration=True, encoding='UTF-8', standalone=True)

sys.path.insert(0, r"/home/user\.workbuddy\plugins\marketplaces\codebuddy-plugins-official\plugins\docx\scripts\office")
cmd = [
    "/path/to/cache/chromium-env/6nqqlh/.workbuddy/binaries/python/envs/default/Scripts/python.exe",
    r"/home/user\.workbuddy\plugins\marketplaces\codebuddy-plugins-official\plugins\docx\scripts\office\pack.py",
    UNPACK_DIR,
    DST_DOCX,
    "--original", SRC_DOCX,
]
result = subprocess.run(cmd, capture_output=True, text=True)
print(f"Repack: {'OK' if result.returncode == 0 else 'ERROR: ' + result.stderr}")
