"""
Fix v5.8 → v5.9 via XML manipulation.
Directly edit the unpacked document.xml and repack.
"""
import shutil
import os
import re
import sys

# Add the docx scripts directory to path
sys.path.insert(0, r"/home/user\.workbuddy\plugins\marketplaces\codebuddy-plugins-official\plugins\docx\scripts\office")

UNPACK_DIR = "/path/to/xelox_project/reports/unpacked_v5.8"
XML_PATH = f"{UNPACK_DIR}/word/document.xml"
DST = "/path/to/xelox_project/reports/manuscript_draft_v5.9_submission.docx"

# Read the XML
with open(XML_PATH, 'r', encoding='utf-8') as f:
    xml = f.read()

# Backup original
with open(XML_PATH + '.backup', 'w', encoding='utf-8') as f:
    f.write(xml)

# Count <w:p> elements (paragraphs) to verify
p_count_before = xml.count('<w:p ') + xml.count('<w:p>')
print(f"Paragraphs before: {p_count_before}")

# --- Fix 1: Remove 2nd occurrence of Discussion duplicate ---
# The two paragraphs are identical. Find both and remove the 2nd one.
disc_text = "comprehensive overfitting diagnostics summary is provided in Supplementary Table S33"
parts = xml.split(disc_text)
print(f"Occurrences of 'overfitting diagnostics': {len(parts) - 1}")

if len(parts) >= 2:
    # We need to find the second occurrence and remove its enclosing <w:p>...</w:p>
    # Find position of 2nd occurrence
    first_pos = xml.find(disc_text)
    second_pos = xml.find(disc_text, first_pos + 1)
    
    # Find the enclosing <w:p> for the 2nd occurrence
    # Search backwards for <w:p and forward for </w:p>
    p_start = xml.rfind('<w:p', 0, second_pos)
    p_end = xml.find('</w:p>', second_pos) + len('</w:p>')
    
    if p_start != -1 and p_end != -1:
        # Check if this is a paragraph that starts with <w:p (not <w:pPr)
        # Extract the paragraph
        dup_para = xml[p_start:p_end]
        # Remove it (plus preceding whitespace/newline)
        # Find the start of preceding whitespace
        trim_start = p_start
        while trim_start > 0 and xml[trim_start - 1] in ' \t\n\r':
            trim_start -= 1
        xml = xml[:trim_start] + xml[p_end:]
        print(f"  Removed duplicate Discussion paragraph ({p_end - p_start} chars)")
    else:
        print("  ERROR: Could not find enclosing <w:p> tags for duplicate")
else:
    print("  ERROR: Discussion duplicate not found")
