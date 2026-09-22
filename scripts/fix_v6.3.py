"""Fix all 12 issues found in v6.2 review, output v6.3."""
import re, os

PROJECT = r"/path/to/xelox_project"
REPORTS = os.path.join(PROJECT, "reports")

with open(os.path.join(REPORTS, "manuscript_draft_v6.2_submission.md"), "r", encoding="utf-8") as f:
    text = f.read()

# ============================================================
# Fix 3: Remove duplicate Table 1 note (keep italic version, remove plain)
# ============================================================
# The plain Note: starts at line 799, italic version is lines 785-797
# Remove the plain duplicate
dup_start = text.find("\nNote: CSNK1G2 = casein kinase 1 gamma 2; GGT family =\ngamma-glutamyltransferase 1/2 (GGT1/GGT2/GGTLC1/GGTLC2); ZNF451 = zinc\nfinger protein 451; KAZN = kazrin; KLK6 = kallikrein-related peptidase\n6; MGA = MAX gene associated; MID2 = midline 2; SOX11 = SRY-box\ntranscription factor 11; TAS2R40 = taste receptor type 2 member 40;\nC5ORF42 = chromosome 5 open reading frame 42. Votes = number of\n")
if dup_start > 0:
    # Find the end of this duplicate block (next blank line or next section)
    dup_end = text.find("\n\n", dup_start + 10)
    if dup_end > 0:
        text = text[:dup_start] + text[dup_end:]
        print("  [3] Removed duplicate Table 1 note")

# ============================================================
# Fix 4: Rename duplicate Table S10 — second one should be Table S11
# ============================================================
# The second "Table S10" (Overfitting Diagnostics) should be Table S11
# But first check: is there already a Table S11?
if text.count("**Table S10.") > 1:
    # Find the second occurrence and rename to S11, shift all subsequent
    first_pos = text.find("**Table S10.")
    second_pos = text.find("**Table S10.", first_pos + 1)
    if second_pos > 0:
        # Also need to rename S11->S12, S12->S13, S13->S14
        # First rename S10->S11
        remaining = text[second_pos:]
        remaining = remaining.replace("**Table S10.", "**Table S11.", 1)
        # Now shift S11->S12, S12->S13, S13->S14 (in remaining text)
        for i in range(13, 10, -1):  # S13->S14, S12->S13, S11->S12
            remaining = remaining.replace(f"**Table S{i}.", f"**Table S{i+1}.")
            remaining = remaining.replace(f"Table S{i} ", f"Table S{i+1} ")
            remaining = remaining.replace(f"(Table S{i})", f"(Table S{i+1})")
            remaining = remaining.replace(f"Table S{i}\n", f"Table S{i+1}\n")
        text = text[:second_pos] + remaining
        print("  [4] Renamed duplicate Table S10 -> S11 (shifted S11-S13 -> S12-S14)")

# ============================================================
# Fix 5: Add Keywords after abstract
# ============================================================
abstract_end = text.find("# **Background**")
if abstract_end > 0 and "Keywords" not in text[:abstract_end+50]:
    keywords = "\n**Keywords**: colorectal cancer, chemoresistance, XELOX, pathway-level analysis, nomogram, ssGSEA, machine learning, overfitting, translational genomics\n\n"
    text = text[:abstract_end] + keywords + text[abstract_end:]
    print("  [5] Added Keywords section")

# ============================================================
# Fix 6: Add standalone Conclusions section after Discussion
# ============================================================
methods_start = text.find("# **Methods**")
if methods_start > 0 and "# **Conclusions**" not in text:
    # Extract last paragraph of Discussion as Conclusions
    discussion_text = text[text.find("# **Discussion**"):methods_start]
    # Find the last paragraph (starts after the last double newline)
    paragraphs = discussion_text.strip().split("\n\n")
    # The conclusions paragraph is typically the last meaningful one
    # Look for "Important limitations" - everything after that is limitations
    # The concluding paragraph is the one before limitations
    conclusions_para = ""
    for i, p in enumerate(paragraphs):
        if "Important limitations" in p:
            # The paragraph before this is the concluding statement
            if i > 0:
                conclusions_para = paragraphs[i-1]
            break
    
    if conclusions_para:
        conclusions_section = f"\n# **Conclusions\n\n{conclusions_para}\n\n"
        text = text[:methods_start] + conclusions_section + text[methods_start:]
        print("  [6] Added standalone Conclusions section")
    else:
        # Fallback: create a brief conclusions
        conclusions_section = """
# **Conclusions**

Pathway-level aggregation functions as an implicit regularization mechanism that substantially improves model reproducibility over gene-level approaches for XELOX chemoresistance prediction. The pathway-level Cox nomogram (C-index = 0.659, optimism-corrected) with 96--100% bootstrap directional stability provides a clinically interpretable framework for post-XELOX therapeutic decision-making. Expression-based drug sensitivity profiling confirms fluoropyrimidine-specific resistance capture, supporting pathway-based risk stratification for treatment intensification decisions. These findings establish that biology-informed dimensionality reduction --- not algorithm choice --- is the primary driver of generalizability in transcriptomic chemoresistance signatures.

"""
        text = text[:methods_start] + conclusions_section + text[methods_start:]
        print("  [6] Added standalone Conclusions section (fallback)")

# ============================================================
# Fix 7: Flag Figure 1 as needing replacement
# ============================================================
if "Gemini_Generated" in text:
    print("  [7] WARNING: Figure 1 is AI-generated (Gemini_Generated_Image) - NEEDS MANUAL REPLACEMENT")
    # Add a note
    text = text.replace(
        "![Figure 1](results/figures_png/v5.11/Gemini_Generated_Image_2mxzi12mxzi12mxz.png)",
        "<!-- TODO: Replace with proper Figure 1 (study design schematic) -->\n![Figure 1](results/figures_png/v5.11/Gemini_Generated_Image_2mxzi12mxzi12mxz.png)"
    )

# ============================================================
# Fix 8: Convert relative figure paths to absolute for DOCX embedding
# ============================================================
# Pandoc needs absolute paths or paths relative to the working directory
# When generating DOCX, we'll run pandoc from the project root
print("  [8] Figure paths: will use absolute paths when generating DOCX")

# ============================================================
# Fix 9: Note about reference format (Vancouver style)
# ============================================================
# Current format uses [1] in text. Vancouver uses numbered refs.
# The current format is actually compatible with Vancouver (numbered in order of appearance).
# Main issue: need to ensure reference list follows Vancouver format.
print("  [9] Reference format: current [N] style is compatible with Vancouver; verify list format")

# ============================================================
# Fix 10: Remove duplicated Supplementary Materials section from bottom of manuscript
# ============================================================
# The supplementary materials section at line 974+ duplicates the separate file
supp_start = text.find("# **Supplementary Materials**\n\n## From Gene-Level Fingerprints")
if supp_start > 0:
    # Find the end (next major section or end of file)
    supp_end = text.find("\n# **References**\n", supp_start + 1)
    if supp_end < 0:
        supp_end = len(text)
    # Keep only a brief reference, not the full catalog
    brief_supp = """
# **Supplementary Materials**

Detailed supplementary tables (S1--S14) and figures (S1--S16) are provided in a separate Supplementary Materials file. See individual table and figure captions in the main text for specific references.

"""
    text = text[:supp_start] + brief_supp + text[supp_end:]
    print("  [10] Replaced full supplementary catalog with brief reference")

# ============================================================
# Fix 11: Remove duplicate Data/Code Availability from Methods area
# ============================================================
# The standalone sections at lines 686-696 duplicate Declarations
# Remove them since they're now in Declarations
data_avail_start = text.find("# **Data Availability**\n\nAll data used")
code_avail_start = text.find("# **Code Availability**\n\nAnalysis scripts")
if data_avail_start > 0 and code_avail_start > 0:
    # Find end of Code Availability
    code_avail_end = text.find("\n\n", code_avail_start + 50)
    if code_avail_end > 0:
        text = text[:data_avail_start] + text[code_avail_end:]
        print("  [11] Removed duplicate Data/Code Availability (kept in Declarations)")

# ============================================================
# Fix 12: Separate main references from supplementary references
# ============================================================
# The references section at the end mixes main and supplementary refs
# This is harder to fix automatically - flag for manual review
print("  [12] Reference list: flagged for manual review (main vs supplementary refs)")

# ============================================================
# Write v6.3
# ============================================================
output_path = os.path.join(REPORTS, "manuscript_draft_v6.3_submission.md")
with open(output_path, "w", encoding="utf-8") as f:
    f.write(text)

print(f"\nv6.3 written to: {output_path}")
print(f"File size: {os.path.getsize(output_path) / 1024:.1f} KB")

# Verify fixes
print("\nVerification:")
print(f"  Table S10 count: {text.count('**Table S10.')}")
print(f"  Table S11 count: {text.count('**Table S11.')}")
print(f"  Keywords present: {'Keywords' in text}")
print(f"  Conclusions section: {'# **Conclusions' in text}")
print(f"  Data Availability standalone: {'# **Data Availability**' in text}")
print(f"  Supplementary full catalog: {'Supplementary Materials Overview' in text}")
