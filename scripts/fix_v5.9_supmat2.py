"""
Fix 4 remaining supplementary material issues in v5.9.
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

# Find all text elements
all_t_elems = list(root.iter(f'{{{ns}}}t'))

for i, t_elem in enumerate(all_t_elems):
    if t_elem.text is None:
        continue
    orig = t_elem.text

    # Fix 1: Round excessive decimal precision in Table S30
    # Oxaliplatin row values
    if len(orig) > 25 and re.match(r'^-?\d+\.\d{10,}', orig.strip()):
        try:
            val = float(orig.strip().split()[0] if ' ' in orig.strip() else orig.strip())
            if abs(val) < 0.01:
                rounded = f"{val:.4f}"
            elif abs(val) < 1:
                rounded = f"{val:.4f}"
            elif abs(val) < 100:
                rounded = f"{val:.4f}"
            else:
                rounded = f"{val:.4f}"
            
            if orig.strip() == str(val) if '.' in orig else False:
                pass  # already reasonable
            elif len(orig.strip()) > 10 and '.' in orig.strip():
                t_elem.text = rounded
                print(f"Fix 1: {orig.strip()[:20]} → {rounded}")
                changes += 1
        except ValueError:
            pass

    # Fix 1b: P-values in Table S30
    if 'e-' in orig and len(orig) < 30 and re.match(r'^[\d.]+e-\d+', orig.strip()):
        try:
            val = float(orig.strip())
            t_elem.text = f"{val:.2e}"
            print(f"Fix 1b: {orig.strip()} → {t_elem.text}")
            changes += 1
        except ValueError:
            pass

# Second pass for text replacements (need context)
for t_elem in all_t_elems:
    if t_elem.text is None:
        continue
    orig = t_elem.text

    # Fix 2: Update outdated "separate file" note in Para 112
    if "Supplementary Tables S1-S16 and Figures S1-S12 are provided in a separate file" in orig:
        new_text = orig.replace(
            "Supplementary Tables S1-S16 and Figures S1-S12 are provided in a separate file (supplementary_materials.docx). Key additions:",
            "A comprehensive catalog of Supplementary Tables S1\u2013S35 and Figures S1\u2013S31 is provided below. Key tables referenced in the main text:"
        )
        t_elem.text = new_text
        print("Fix 2: Updated 'separate file' note in Para 112")
        changes += 1

    # Fix 3: n=536 clarification in Figure S30/S31 descriptions
    if "n = 536 complete cases after filtering time=0" in orig:
        t_elem.text = orig.replace(
            "n = 536 complete cases after filtering time=0",
            "n = 536 (all complete cases after excluding zero-time records; main Cox nomogram uses n = 227 XELOX subset)"
        )
        print("Fix 3: Clarified n=536 in Figure S30")
        changes += 1

    if "Data: GSE39582 XELOX (n = 536). Elastic Net" in orig:
        t_elem.text = orig.replace(
            "Data: GSE39582 XELOX (n = 536). Elastic Net",
            "Data: GSE39582 XELOX (n = 536 complete cases; main Cox nomogram uses n = 227 XELOX subset). Elastic Net"
        )
        print("Fix 3b: Clarified n=536 in Figure S31")
        changes += 1

    # Fix 4: Table S33 annotation precision
    if "Comparison between gene-level (500-gene LASSO-Cox) and pathway-level (44-pathway AIC-Cox nomogram) models" in orig:
        t_elem.text = orig.replace(
            "Comparison between gene-level (500-gene LASSO-Cox) and pathway-level (44-pathway AIC-Cox nomogram) models",
            "Diagnostics for the 8-variable pathway-based Cox nomogram (7 pathway scores + tumor location) in GSE39582 XELOX subset (n = 227, 79 events)"
        )
        print("Fix 4: Corrected Table S33 annotation")
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
    "--validate", "false",
]
result = subprocess.run(cmd, capture_output=True, text=True)
print(f"Repack: {'OK' if result.returncode == 0 else 'ERROR: ' + result.stderr}")
