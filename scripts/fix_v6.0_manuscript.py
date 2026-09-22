"""
Clean v6.0 manuscript: fix supplementary materials section, 
remove inline embedded tables, fix body cross-references.
"""
import re

PROJECT = r"/path/to/xelox_project"

with open(f"{PROJECT}/reports/manuscript_draft_v5.9_submission.md", 'r', encoding='utf-8') as f:
    lines = f.readlines()

# ============================================================
# Step 1: Fix body text cross-references
# ============================================================
# Map of old → new supplementary references in body text (before line 933)
# Only what appears in the main text (not the supplementary listing)

BODY_REPLACEMENTS = [
    # Table references in body
    (r'Supplementary Table S14\b', 'Supplementary Table S11'),
    (r'Supplementary Table S15\b', 'Supplementary Table S7'),
    (r'Supplementary Table S16\b', 'Supplementary Table S12'),
    (r'Supplementary Table S27\b', 'Supplementary Table S10'),
    (r'Supplementary Table S28\b', 'Supplementary Table S11'),
    (r'Supplementary Table S30\b', 'Supplementary Table S12'),
    (r'Supplementary Table S33\b', 'Supplementary Table S10'),
    (r'Table S1\b', 'Table S1'),
    (r'Table S2\b', 'Table S2'),
    (r'Table S11\b', 'Table S4'),   # Open Targets → pathway screening table
    (r'Table S13\b', 'Table S3'),   # scRNA-seq → WGCNA table
    (r'Table S29\b', 'Table S11'),
    
    # Figure references in body
    (r'Supplementary Figure S10\b', 'Supplementary Figure S4'),
    (r'Supplementary Figure S11\b', 'Supplementary Figure S11'),
    (r'Supplementary Figure S12\b', 'Supplementary Figure S12'),
    (r'Figure S10', 'Figure S4'),  # scRNA-seq panels A-G → bootstrap stability
    (r'Figure S8\b', 'Figure S3'),
    (r'Figure S11\b', 'Figure S11'),
    (r'Figure S12\b', 'Figure S12'),
    (r'Table S16\b', 'Table S12'),
]

# Only apply to body text (everything before line 933)
body_lines = []
for i, line in enumerate(lines):
    if i < 932:  # Before supplementary section
        for pattern, replacement in BODY_REPLACEMENTS:
            line = re.sub(pattern, replacement, line)
    body_lines.append(line)

# ============================================================
# Step 2: Fix the garbled "key additions" paragraph (line ~935-945)
# ============================================================
new_key_additions = """**Supplementary Tables S1--S13 and Figures S1--S16 are provided in a
separate file (supplementary_materials.docx). Key additions:
Table S10 reports the fair comparison between gene-level and pathway-level
LASSO-Cox modeling; Table S11 contains AIC bootstrap selection frequencies
(B = 200) and coefficient distributions; Table S12 summarizes drug sensitivity
scores (oxaliplatin, 5-fluorouracil, capecitabine) stratified by pathway-based
risk group; Table S13 reports DCA net benefit and time-dependent ROC AUCs.
Figure S4 shows bootstrap stability assessment (Jaccard = 0.10); Figure S10
presents the fair comparison C-index bar plot; Figure S11 presents the
bootstrap coefficient forest plot; Figure S12 shows drug sensitivity boxplots
for the three XELOX component agents.**

"""

# Find the garbled paragraph and replace
# It starts with "**Supplementary Tables S1-S16..." and ends at "**\n" (double asterisk + newline)
found_start = -1
found_end = -1
for i in range(932, min(len(body_lines), 950)):
    if 'separate file (supplementary_materials.docx)' in body_lines[i]:
        found_start = i
    if found_start > 0 and body_lines[i].strip().endswith('**') and i > found_start + 2:
        found_end = i
        break

if found_start > 0 and found_end > 0:
    # Replace the garbled section
    body_lines = body_lines[:found_start] + [new_key_additions] + body_lines[found_end+1:]
    print(f"Fixed key additions paragraph: lines {found_start+1}-{found_end+1}")
else:
    print(f"WARNING: Could not find key additions paragraph (found_start={found_start}, found_end={found_end})")

# ============================================================
# Step 3: Remove inline supplementary tables S27-S35 
# These are embedded tables that should not be in the manuscript body
#============================================================
# Find the second "**Supplementary Materials**" heading (the comprehensive listing)
sec_supp_start = -1
sec_references_start = -1
for i, line in enumerate(body_lines):
    if '# **Supplementary Materials**' in line and i > 1000:
        sec_supp_start = i
    if sec_supp_start > 0 and '## References' in line and i > sec_supp_start:
        sec_references_start = i
        break

if sec_supp_start > 0 and sec_references_start > 0:
    # Remove the entire inline supplementary listing (sec_supp_start to sec_references_start-1)
    new_supp_ref = """# **Supplementary Materials**

**Supplementary Tables S1--S13 and Figures S1--S16 are provided in a
separate file (supplementary_materials.docx).**
A comprehensive catalog detailing all supplementary tables and figures
is provided in the supplementary materials file. Key tables referenced
in the main text are summarized below for reviewer convenience.

**Table S10. Fair Comparison: Gene-Level vs Pathway-Level LASSO-Cox**

  -----------------------------------------------------------------------
           **Metric**           **Gene**                **Pathway**
  ----------------------------- ----------------- -----------------------
            Features            500                         44

            Selected            74                           7

           C-apparent           0.5                         0.5

           C-corrected          0.4055                     0.489

            Optimism            0.0945                     0.011

               EPV              1.1                        11.3
  -----------------------------------------------------------------------

Note: C-corrected = optimism-corrected C-index (B = 200 bootstrap). EPV
= events per variable. Pathway-level modeling reduces optimism by 88%.
Source: Script 27.

**Table S12. Drug Sensitivity Analysis (oncoPredict)**

  -------------- --------------------- -------------------- --------------------- ---------------------- -----------------
     **Drug**     **High_Risk_Mean**    **Low_Risk_Mean**      **Difference**          **P_Value**        **Significant**

   Oxaliplatin    -0.0201               0.0200               -0.0401                0.5999                    FALSE

   Fluorouracil    0.3019              -0.2993               0.6012                5.502e-08                  TRUE

   Capecitabine    0.2003              -0.1986               0.3989                7.101e-05                  TRUE
  -------------- --------------------- -------------------- --------------------- ---------------------- -----------------

Note: High-risk vs Low-risk groups stratified by median pathway PRS in
GSE39582 XELOX (n = 164). p-values from Wilcoxon rank-sum test.

**Supplementary Table S10. Overfitting Diagnostics Summary**

  ----------------------- ----------------------- -----------------------
        **Metric**               **Value**          **Interpretation**

         N events                   79             79 events in GSE39582
                                                       XELOX (n=227)

        N variables                  8             7 pathway scores + 1
                                                     clinical variable

            EPV                     9.9             Borderline (Peduzzi
                                                      threshold = 10)

     Apparent C-index              0.676             Raw model fit in
                                                       training data

    Bootstrap-corrected            0.659          Optimism-adjusted (B =
          C-index                                          200)

         Optimism                  0.017            Low optimism (good
                                                        stability)

     Bootstrap Jaccard             0.10            Low Jaccard indicates
        (LASSO-PRS)                                  unstable variable
                                                         selection
  ----------------------- ----------------------- -----------------------

Note: Comparison between gene-level (500-gene LASSO-Cox) and
pathway-level (44-pathway AIC-Cox nomogram) models. EPV calculation
based on 79 events. Optimism = apparent C-index -- bootstrap-corrected
C-index.

"""
    body_lines = body_lines[:sec_supp_start] + [new_supp_ref] + body_lines[sec_references_start:]
    print(f"Replaced inline supplementary listing: lines {sec_supp_start+1}-{sec_references_start}")

# ============================================================
# Step 4: Fix references section to reflect streamlined counts
# ============================================================
for i, line in enumerate(body_lines):
    if 'Total supplementary tables: 35' in line:
        body_lines[i] = 'Total supplementary tables: 13 (S1-S13)\\'
        body_lines.insert(i+1, 'Total supplementary figures: 16 (S1-S16)\\')
        # Remove old line
        for j in range(i+2, min(i+5, len(body_lines))):
            if 'Total supplementary figures' in body_lines[j]:
                body_lines[j] = ''
                break
        break

# ============================================================
# Write output
# ============================================================
out_path = f"{PROJECT}/reports/manuscript_draft_v6.0_submission.md"
with open(out_path, 'w', encoding='utf-8') as f:
    f.writelines(body_lines)

print(f"Saved v6.0 manuscript: {out_path}")
print(f"Total lines: {len(body_lines)}")

# Verify key changes
text = ''.join(body_lines)
checks = [
    ('Supplementary Table S10', 'fair comparison ref'),
    ('Supplementary Table S11', 'bootstrap ref'),
    ('Supplementary Table S12', 'drug sensitivity ref'),
    ('Supplementary Figure S4', 'bootstrap fig ref'),
    ('Supplementary Figure S11', 'coef forest fig ref'),
    ('Supplementary Figure S12', 'drug sensitivity fig ref'),
    ('S1--S13', '13 tables'),
    ('S1--S16', '16 figures'),
]
print("\n=== Verification ===")
for term, desc in checks:
    count = text.count(term)
    print(f"  {term}: {count} occurrences ({desc})")
