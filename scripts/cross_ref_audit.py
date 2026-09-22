"""
Cross-reference audit: Extract all key numbers and verify consistency
across Abstract, Results, Discussion, Tables.
"""
import re

with open("/path/to/xelox_project/reports/manuscript_draft_v5.9.md", 'r', encoding='utf-8') as f:
    content = f.read()

# Split into sections
lines = content.split('\n')

# Find section boundaries
sections = {}
current_section = "Header"
for i, line in enumerate(lines):
    if line.startswith('# **Abstract**'):
        current_section = "Abstract"
        sections[current_section] = i
    elif line.startswith('# **Introduction**'):
        current_section = "Introduction"
        sections[current_section] = i
    elif line.startswith('# **Results**'):
        current_section = "Results"
        sections[current_section] = i
    elif line.startswith('# **Discussion**'):
        current_section = "Discussion"
        sections[current_section] = i
    elif line.startswith('# **Methods**'):
        current_section = "Methods"
        sections[current_section] = i
    elif line.startswith('# **Tables**'):
        current_section = "Tables"
        sections[current_section] = i
    elif line.startswith('# **References**'):
        current_section = "References"
        sections[current_section] = i
    elif line.startswith('# **Supplementary'):
        current_section = "Supplementary"
        sections[current_section] = i

print("=== Section boundaries ===")
for sec, line_num in sections.items():
    print(f"  {sec}: line {line_num}")

# Extract key numbers with context
key_patterns = [
    "n = 164", "n = 227", "n = 476", "n = 1,015",
    "EPV = 3.4", "EPV = 9.9", "EPV = 1.1", "EPV = 11.3",
    "C-index = 0.676", "C-index = 0.834", "C-index = 0.489", "C-index = 0.406",
    "corrected C-index = 0.659", "corrected C-index = 0.489", "corrected C-index = 0.406",
    "optimism = 0.017", "optimism = 0.011", "optimism = 0.0945",
    "AUC = 0.899", "AUC = 0.54", "AUC = 0.49", "AUC = 0.60",
    "B = 200", "B = 1,000", "B = 25",
    "66.5%", "59.0%", "49.0%",
    "79 events", "79 events",
    "23 genes", "23-gene",
    "8-variable", "8 variable",
    "10-gene",
    "253 genes",
    "5-fluorouracil: p", "capecitabine: p",
    "p = 0.600", "p < 0.0001",
]

print("\n=== Key number occurrences ===")
for pattern in key_patterns:
    matches = []
    for i, line in enumerate(lines):
        if pattern in line:
            # Determine which section
            sec = "HDR"
            for s_name, s_line in sorted(sections.items(), key=lambda x: x[1]):
                if i >= s_line:
                    sec = s_name
            matches.append((i, sec, line.strip()[:100]))
    
    if len(matches) > 1 or len(matches) == 0:
        if len(matches) > 1:
            print(f"  '{pattern}': {len(matches)}x")
        elif "25" not in pattern and "XELOX-like" not in pattern:
            pass  # skip some non-critical patterns

# Check for reference number consistency
print("\n=== Reference number gaps ===")
refs_used = set()
for m in re.finditer(r'\\\[(\d+(?:,\d+)*)\\\]', content):
    for num in m.group(1).split(','):
        refs_used.add(int(num.strip()))
all_refs = sorted(refs_used)
print(f"  References used: {all_refs}")
print(f"  Range: [{min(all_refs)}-{max(all_refs)}]")
missing = [i for i in range(min(all_refs), max(all_refs)+1) if i not in all_refs]
if missing:
    print(f"  GAPS: {missing}")

# Check for consistent drug name usage
print("\n=== Drug name consistency ===")
for drug in ["XELOX", "FOLFOX", "5-fluorouracil", "5-FU", "capecitabine", "oxaliplatin"]:
    count = content.count(drug)
    if count > 0:
        print(f"  {drug}: {count}x")
