"""
v6.0: Update manuscript v5.9 → v6.0 with streamlined supplementary numbering.
Map: 35 tables + 31 figures → 13 tables + 16 figures
"""
import re

# ============================================================
# OLD → NEW mapping for Supplementary Tables
# ============================================================
TABLE_MAP = {
    # Keep same: S1 (meta-analysis)
    # S2-S6 merged into new S2 (5-cohort DEG summary)
    'S3': 'S2', 'S4': 'S2', 'S5': 'S2', 'S6': 'S2',
    'S7': 'S3',   # WGCNA hub genes
    'S8': 'S4',   # Univariate Cox + BH
    'S9': 'S5',   # Pathway consistency
    'S10': None,  # DELETED: per-cohort pathway diff (redundant with S5)
    'S11': 'S4',  # Merged into S4 (bootstrap freq)
    'S12': None,  # DELETED: merged PRS not final model
    'S13': None,  # DELETED: 21K gene univariate (noise)
    'S14': 'S6',  # 23-gene AIC Cox
    'S15': 'S7',  # Ablation comparison
    'S16': 'S8',  # Merged into S8 (PH)
    'S17': 'S8',  # Merged into S8 (dfbeta)
    'S18': 'S8',  # Merged into S8 (LOVO)
    'S19': 'S8',  # Merged into S8 (stratified)
    'S20': 'S8',  # Merged into S8 (cutoff) -- kept in S8
    'S21': 'S8',  # Merged into S8 (bootstrap CI)
    'S22': None,  # DELETED: GSE83129 n=29
    'S23': 'S9',  # TCGA validation
    'S24': None,  # DELETED: 2.5K gene votes
    'S25': None,  # DELETED: GDSC pre-filter
    'S26': 'S1',  # Merged into S1 (direction matrix)
    'S27': 'S10', # Fair comparison
    'S28': 'S11', # AIC bootstrap freq
    'S29': 'S11', # Merged into S11 (coef distribution)
    'S30': 'S12', # Drug sensitivity
    'S31': 'S8',  # Merged into S8 (calibration)
    'S32': 'S13', # DCA + TDROC
    'S33': None,  # DELETED: duplicate of S10
    'S34': None,  # DELETED: Ridge
    'S35': None,  # DELETED: Elastic Net
}

# ============================================================
# OLD → NEW mapping for Supplementary Figures
# ============================================================
FIGURE_MAP = {
    'S1': 'S1',   # WGCNA
    'S2': 'S2',   # ComBat PCA
    'S3': None,   # DELETED
    'S4': None,   # DELETED
    'S5': None,   # DELETED
    'S6': None,   # DELETED
    'S7': None,   # DELETED
    'S8': 'S3',   # Pathway analysis
    'S9': 'S3',   # Merged into S3
    'S10': 'S4',  # Bootstrap stability
    'S11': None,  # DELETED
    'S12': None,  # DELETED
    'S13': None,  # DELETED
    'S14': None,  # DELETED
    'S15': 'S5',  # Sensitivity diagnostics
    'S16': None,  # DELETED
    'S17': None,  # DELETED
    'S18': None,  # DELETED
    'S19': 'S6',  # Gene-level Cox forest
    'S20': None,  # DELETED
    'S21': 'S7',  # TCGA validation
    'S22': 'S8',  # Meta-analysis
    'S23': 'S9',  # TGF-beta paradox
    'S24': 'S10', # Fair comparison
    'S25': 'S11', # Bootstrap coef forest
    'S26': 'S12', # Drug sensitivity
    'S27': 'S13', # Calibration
    'S28': 'S14', # DCA
    'S29': 'S7',  # Merged into S7 (TDROC)
    'S30': None,  # DELETED
    'S31': None,  # DELETED
}

# ============================================================
# Patterns for replacement in manuscript body
# ============================================================
# We use regex to find and replace references like:
# "Supplementary Table S27" → "Supplementary Table S10"
# "Table S27" → "Table S10"
# "Supplementary Figure S10" → "Supplementary Figure S4"
# etc.

def build_replacement_rules():
    """Build ordered replacement rules (longer patterns first to avoid partial matches)."""
    rules = []
    
    # Table replacements
    for old_num, new in sorted(TABLE_MAP.items(), key=lambda x: -int(x[0][1:])):
        old_num_stripped = old_num  # e.g., "S27"
        
        # Pattern: "Supplementary Table SXX"
        if new is not None:
            rules.append((
                re.compile(r'Supplementary Table ' + old_num_stripped + r'\b'),
                'Supplementary Table ' + new
            ))
            rules.append((
                re.compile(r'Table ' + old_num_stripped + r'\b'),
                'Table ' + new
            ))
            # For S27-S29 style ranges
            # Handle separately below
    
    # Figure replacements
    for old_num, new in sorted(FIGURE_MAP.items(), key=lambda x: -int(x[0][1:])):
        old_num_stripped = old_num
        if new is not None:
            rules.append((
                re.compile(r'Supplementary Figure ' + old_num_stripped + r'\b'),
                'Supplementary Figure ' + new
            ))
            rules.append((
                re.compile(r'Figure ' + old_num_stripped + r'\b'),
                'Figure ' + new
            ))
    
    return rules

# Handle ranges specially:
# "Supplementary Table S27-S29" → "Supplementary Tables S10-S11"
SPECIAL_RULES = [
    # Table ranges
    (re.compile(r'Supplementary Table S27[\–\-]S29\b'), 'Supplementary Tables S10 and S11'),
    (re.compile(r'Table S27[\–\-]S29\b'), 'Tables S10 and S11'),
    (re.compile(r'Table S1[\–\-]S6\b'), 'Tables S1 and S2'),
    
    # Removed entries - add note
    (re.compile(r'Supplementary Table S33\b'), 'Supplementary Table S10'),
    (re.compile(r'Table S33\b'), 'Table S10'),
    (re.compile(r'Supplementary Table S34\b'), '(see Table S10 for overfitting diagnostics)'),
    (re.compile(r'Table S34\b'), '(see Table S10)'),
    (re.compile(r'Supplementary Table S35\b'), '(see Table S10 for overfitting diagnostics)'),
    (re.compile(r'Table S35\b'), '(see Table S10)'),
    (re.compile(r'Table S11\b(?![\–\-])'), 'Table S3'),  # Open Targets → WGCNA context; handle carefully
]

def replace_in_manuscript(text):
    """Apply all replacement rules to manuscript text."""
    rules = build_replacement_rules()
    
    # Apply special range rules first
    for pattern, replacement in SPECIAL_RULES:
        text = pattern.sub(replacement, text)
    
    # Apply individual table/figure rules
    for pattern, replacement in rules:
        text = pattern.sub(replacement, text)
    
    return text

# ============================================================
# Test
# ============================================================
if __name__ == '__main__':
    PROJECT = r"/path/to/xelox_project"
    md_path = f"{PROJECT}/reports/manuscript_draft_v5.9_submission.md"
    
    with open(md_path, 'r', encoding='utf-8') as f:
        text = f.read()
    
    # Apply replacements
    updated = replace_in_manuscript(text)
    
    # Check changes
    old_refs = set()
    new_refs = set()
    for m in re.finditer(r'(?:Supplementary )?(?:Table|Figure) S\d+', text):
        old_refs.add(m.group())
    for m in re.finditer(r'(?:Supplementary )?(?:Table|Figure) S\d+', updated):
        new_refs.add(m.group())
    
    print("=== OLD references in body ===")
    for r in sorted(old_refs):
        print(f"  {r}")
    print()
    print("=== NEW references in body ===")
    for r in sorted(new_refs):
        print(f"  {r}")
    
    # Save preview
    out_path = f"{PROJECT}/reports/manuscript_draft_v6.0_submission.md"
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write(updated)
    print(f"\nSaved: {out_path}")
    print(f"Size: {len(updated)} chars")
