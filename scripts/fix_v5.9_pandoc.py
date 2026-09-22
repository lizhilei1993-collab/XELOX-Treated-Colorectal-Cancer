"""
Fix v5.8 → v5.9 via pandoc roundtrip + targeted fixes.

Approach:
1. Convert docx → markdown (already done as manuscript_draft_v5.8_submission.md)
2. Apply fixes in markdown  
3. Use pandoc reference-docx to preserve styles
4. Convert back to docx

This preserves more formatting than raw pandoc.
"""
import subprocess
import os

WORK = "/path/to/xelox_project"
SRC_MD = f"{WORK}/reports/manuscript_draft_v5.8_submission.md"
DST_DOCX = f"{WORK}/reports/manuscript_draft_v5.9_submission.docx"
SRC_DOCX = f"{WORK}/reports/manuscript_draft_v5.8_submission.docx"

# Read the markdown
with open(SRC_MD, 'r', encoding='utf-8') as f:
    content = f.read()

lines = content.split('\n')
print(f"Total lines: {len(lines)}")

# --- Go through the markdown and fix issues ---

# Fix 1: Remove the duplicated Discussion paragraph
# The duplicate starts at line ~646 and goes to ~668
# Actually lines 646-668 are the same as 602-624
# Let me identify the exact lines

# First, find all occurrences of the marker
marker = "A comprehensive overfitting diagnostics summary is provided in Supplementary Table S33."
occurrences = []
for i, line in enumerate(lines):
    if marker in line:
        occurrences.append(i)
print(f"'{marker[:60]}...' found at lines: {occurrences}")

if len(occurrences) >= 2:
    # The first occurrence is the original, second is the duplicate
    # Find paragraph boundaries around the second occurrence
    first = occurrences[0]
    second = occurrences[1]
    
    # Search backwards from the second occurrence for paragraph break
    para_start = second
    for j in range(second - 1, max(0, second - 30), -1):
        if lines[j].strip() == '':
            para_start = j + 1
            break
    
    # Search forward for end of this paragraph
    para_end = second
    for j in range(second + 1, min(len(lines), second + 50)):
        if lines[j].strip() == '':
            para_end = j
            break
    
    print(f"  Duplicate paragraph: lines {para_start}-{para_end}")
    print(f"  First 80 chars: {lines[para_start][:80]}")
    
    # Remove these lines (and preceding empty line)
    del_start = max(0, para_start - 1)  # include empty line before
    del lines[del_start:para_end]
    print(f"  Removed lines {del_start}-{para_end-1}")
else:
    print("  ERROR: Could not find duplicate")

# Fix 2: Remove the conflicting Table 2 note
# Find "Note: HR = hazard ratio; CI = confidence interval." lines
table2_marker = "Note: HR = hazard ratio; CI = confidence interval."
table2_occurrences = []
for i, line in enumerate(lines):
    if line.startswith(table2_marker):
        table2_occurrences.append(i)
print(f"\nTable 2 'Note:' found at lines: {table2_occurrences}")

# The first occurrence has VIF info and conflicts with second (leave-one-variable-out)
# We keep the second (more accurate) and remove the first
if len(table2_occurrences) >= 2:
    first_note = table2_occurrences[0]
    second_note = table2_occurrences[1]
    
    # Find paragraph boundary for first note
    note_start = first_note
    for j in range(first_note - 1, max(0, first_note - 10), -1):
        if lines[j].strip() == '':
            note_start = j + 1
            break
    
    note_end = first_note
    for j in range(first_note + 1, min(len(lines), first_note + 20)):
        if lines[j].strip() == '':
            note_end = j
            break
    
    print(f"  Removing conflicting note: lines {note_start}-{note_end}")
    del lines[note_start:note_end]
    print(f"  Removed {note_end - note_start} lines")
else:
    print(f"  Found {len(table2_occurrences)} occurrences")

# Fix 3: Remove the duplicate Supplementary Materials (3rd copy)
# Find all "Supplementary Materials" headers (bold, level-1)
supp_markers = []
for i, line in enumerate(lines):
    if line.strip() == '**Supplementary Materials**':
        supp_markers.append(i)
print(f"\n'**Supplementary Materials**' found at lines: {supp_markers}")

# The original markdown had 3 occurrences but pandoc may render them differently
# Let's also check for non-bold versions
for i, line in enumerate(lines):
    if line.strip() == 'Supplementary Materials' and i not in supp_markers:
        supp_markers.append(i)
print(f"All 'Supplementary Materials' lines: {supp_markers}")

# Find the last Supplementary Materials block and remove everything after it
if len(supp_markers) >= 3:
    last_supp = supp_markers[-1]
    print(f"  Removing everything from line {last_supp} onwards (3rd Supplementary Materials block)")
    lines = lines[:last_supp]
    print(f"  Truncated to {len(lines)} lines")
elif len(supp_markers) == 2:
    # Check if the 2nd block ends and a 3rd starts elsewhere
    print(f"  Found 2 '**Supplementary Materials**' (expected 3). Checking for plain version...")
    # The 3rd could be a non-markedup version
    for i in range(max(supp_markers) + 1, len(lines)):
        if lines[i].strip() == '**Supplementary Tables:**' or 'Supplementary Tables: S1' in lines[i]:
            print(f"  Found supplementary tables header at line {i}, checking if duplicate")
    print(f"  Keeping both blocks for now")
else:
    print(f"  Found {len(supp_markers)} blocks, nothing to remove")

# Write the fixed markdown
fixed_md = f"{WORK}/reports/manuscript_draft_v5.9_fixed.md"
with open(fixed_md, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))
print(f"\nFixed markdown saved to: {fixed_md}")

# Convert to docx using pandoc with reference-docx
print("\nConverting to docx...")
cmd = [
    'pandoc',
    fixed_md,
    '-o', DST_DOCX,
    '--from', 'markdown',
    '--reference-doc', SRC_DOCX,
]
result = subprocess.run(cmd, capture_output=True, text=True)
if result.returncode == 0:
    print(f"Successfully created: {DST_DOCX}")
else:
    print(f"Error: {result.stderr}")
    # Try without reference-doc
    cmd2 = ['pandoc', fixed_md, '-o', DST_DOCX, '--from', 'markdown']
    result2 = subprocess.run(cmd2, capture_output=True, text=True)
    if result2.returncode == 0:
        print(f"Successfully created (without reference): {DST_DOCX}")
    else:
        print(f"Error: {result2.stderr}")
