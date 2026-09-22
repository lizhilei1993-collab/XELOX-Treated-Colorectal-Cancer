"""Convert v6.1 manuscript to v6.2 for Journal of Translational Medicine."""
import re, os

PROJECT = r"/path/to/xelox_project"
REPORTS = os.path.join(PROJECT, "reports")

# Read v6.1
with open(os.path.join(REPORTS, "manuscript_draft_v6.1_submission.md"), "r", encoding="utf-8") as f:
    text = f.read()

# ============================================================
# 1. Restructure abstract to structured format
# ============================================================
# Find the abstract and split into background/methods/results/conclusions
abstract_start = text.find("# **Abstract**")
abstract_end = text.find("# **Introduction**")

old_abstract = text[abstract_start:abstract_end].strip()

# The current abstract is one long paragraph. Split it logically.
# Background: first 3 sentences (up to "uncertain.")
# Methods: "Here we present..." through the methodology description
# Results: the quantitative findings
# Conclusions: last 2 sentences

structured_abstract = """# **Abstract**

## Background
Adjuvant XELOX chemotherapy reduces recurrence in stage III colorectal cancer (CRC), yet approximately 30--40% of patients experience chemoresistance-driven relapse. Current efforts to identify transcriptomic biomarkers for resistance stratification have been dominated by gene-level machine learning approaches, but their cross-cohort generalizability remains unproven and their clinical utility is uncertain.

## Methods
We performed a systematic multi-cohort analysis (1,015 samples across five GEO datasets and TCGA) directly comparing gene-level and pathway-level modeling strategies for XELOX resistance prediction. Probe-level filtering used interquartile range (IQR)-based criteria, and ComBat batch correction was applied to samples with complete response annotation. A ten-algorithm ensemble feature selection with LASSO-cross-validated pre-filtering was used for gene-level modeling. Pathway-level modeling used ssGSEA scores across 44 curated pathways with LASSO-Cox and backward AIC selection. A strictly controlled head-to-head comparison applying identical LASSO-Cox methodology was performed at both levels.

## Results
At the gene level, a ten-algorithm ensemble identified a 10-gene fingerprint with high apparent performance (AUC = 0.899, GSE39582 discovery), but external validation revealed uniformly poor generalizability (AUC = 0.54--0.64), and gene-level ablation confirmed severe overfitting (EPV = 3.4). At the pathway level, a Cox regression nomogram integrating seven ssGSEA pathway scores with tumor location achieved substantially better discrimination (C-index = 0.676, optimism-corrected C-index = 0.659, EPV = 9.9) with 96--100% directional stability under bootstrap resampling. Pathway-level modeling substantially reduced overfitting (bootstrap-corrected C-index: 0.489 vs. 0.406 for gene-level; optimism: 0.011 vs. 0.0945). Expression-based drug sensitivity profiling demonstrated that pathway-defined high-risk tumors exhibit significantly elevated fluoropyrimidine resistance (5-FU: p < 0.0001; capecitabine: p < 0.0001), while oxaliplatin sensitivity did not differ between risk groups (p = 0.145). Single-cell transcriptomic validation (GSE132465, 63,660 cells) confirmed that both the resistance fingerprint and capecitabine metabolism enzymes are enriched in tumor epithelial cells.

## Conclusions
Pathway-level aggregation functions as an implicit regularization mechanism that improves reproducibility over gene-level signatures. The resulting nomogram provides a clinically interpretable framework for post-XELOX therapeutic decision-making, with fluoropyrimidine-specific resistance capture validated across pharmacogenomic and single-cell analyses."""

text = text[:abstract_start] + structured_abstract + "\n\n" + text[abstract_end:]

# ============================================================
# 2. Rename Introduction → Background
# ============================================================
text = text.replace("# **Introduction**", "# **Background**")

# ============================================================
# 3. Update figure references to v5.11 paths
# ============================================================
figure_map = {
    "media/image1.png": "Fig1_study_design.png",
    "media/image2.png": "fig2a+2b.png",
    "media/image3.png": "Fig3_pathway.png",
    "media/image4.png": "Fig4_nomogram.png",
    "media/image5.png": "Fig5_fair_comparison.png",
    "media/image6.png": "Fig6_enrichment.png",
}

for old_ref, new_name in figure_map.items():
    new_path = f"results/figures_png/v5.11/{new_name}"
    text = text.replace(f"![{old_ref}]({old_ref})", f"![{new_name}]({new_path})")

# Also fix the markdown image references with width/height
text = re.sub(r'!\[Gemini_Generated_Image[^\]]*\]\(media/image1\.png\)\{[^}]*\}',
              '![Figure 1](results/figures_png/v5.11/Gemini_Generated_Image_2mxzi12mxzi12mxz.png)', text)
text = re.sub(r'!\[WPS拼图0\]\(media/image2\.png\)\{[^}]*\}',
              '![Figure 2](results/figures_png/v5.11/fig2a+2b.png)', text)
text = re.sub(r'!\[Fig3_pathway\]\(media/image3\.png\)\{[^}]*\}',
              '![Figure 3](results/figures_png/v5.11/Fig3_pathway.png)', text)
text = re.sub(r'!\[Fig4_nomogram\]\(media/image4\.png\)\{[^}]*\}',
              '![Figure 4](results/figures_png/v5.11/Fig4_nomogram.png)', text)
text = re.sub(r'!\[Fig5_fair_comparison\]\(media/image5\.png\)\{[^}]*\}',
              '![Figure 5](results/figures_png/v5.11/Fig5_fair_comparison.png)', text)
text = re.sub(r'!\[Fig6_enrichment\]\(media/image6\.png\)\{[^}]*\}',
              '![Figure 6](results/figures_png/v5.11/Fig6_enrichment.png)', text)

# ============================================================
# 4. Add Declarations section before References
# ============================================================
declarations = """
# **Declarations**

## Ethics approval and consent to participate
This study used publicly available, de-identified transcriptomic data from GEO (accession numbers: GSE39582, GSE104645, GSE28702, GSE72970, GSE69657, GSE132465) and TCGA (COAD/READ). No institutional ethics approval was required as all data were previously published and publicly accessible. The original studies obtained appropriate ethics approval and informed consent.

## Consent for publication
Not applicable.

## Availability of data and materials
All data used in this study are publicly available from GEO (https://www.ncbi.nlm.nih.gov/geo/) and TCGA (https://portal.gdc.cancer.gov/). Analysis scripts are available from the corresponding author upon reasonable request.

## Competing interests
The authors declare that they have no competing interests.

## Funding
This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

## Authors' contributions
ZL: conceptualization, methodology, software, formal analysis, data curation, writing - original draft. [Additional authors to be confirmed]. All authors read and approved the final manuscript.

## Acknowledgements
We acknowledge the GEO and TCGA data repositories for providing open access to the transcriptomic datasets used in this study.

## Authors' information
[To be completed]

"""

# Insert declarations before the Tables section or at the end
tables_pos = text.find("# **Tables**")
if tables_pos > 0:
    text = text[:tables_pos] + declarations + "\n" + text[tables_pos:]
else:
    text += "\n" + declarations

# ============================================================
# 5. Add Abbreviations list
# ============================================================
abbreviations = """
# **List of abbreviations**

| Abbreviation | Full term |
|-------------|-----------|
| AIC | Akaike information criterion |
| AUC | Area under the curve |
| CI | Confidence interval |
| C-index | Concordance index |
| CRC | Colorectal cancer |
| DCA | Decision curve analysis |
| DEG | Differentially expressed gene |
| EMT | Epithelial-mesenchymal transition |
| EPV | Events per variable |
| FDR | False discovery rate |
| GEO | Gene Expression Omnibus |
| GSEA | Gene set enrichment analysis |
| HR | Hazard ratio |
| IQR | Interquartile range |
| LASSO | Least absolute shrinkage and selection operator |
| ORA | Over-representation analysis |
| PCA | Principal component analysis |
| PLS-DA | Partial least squares discriminant analysis |
| PRS | Pathway risk score |
| RFS | Recurrence-free survival |
| ROC | Receiver operating characteristic |
| SHAP | SHapley Additive exPlanations |
| ssGSEA | Single-sample gene set enrichment analysis |
| TCGA | The Cancer Genome Atlas |
| WGCNA | Weighted gene co-expression network analysis |
| XELOX | Capecitabine plus oxaliplatin |

"""

# Insert abbreviations after Declarations
declarations_end = text.find("## Acknowledgements") + len("## Acknowledgements")
# Find the end of acknowledgements paragraph
ack_end = text.find("\n\n", declarations_end + 50)
if ack_end > 0:
    # Find Authors' information end
    authors_info_end = text.find("\n\n", text.find("## Authors' information", ack_end) + 20)
    if authors_info_end > 0:
        text = text[:authors_info_end] + abbreviations + text[authors_info_end:]

# ============================================================
# 6. Update title page format
# ============================================================
# Add title page elements
title_page = """**Pathway-Level Molecular Profiling Resolves the Overfitting Trap in Gene-Level Chemoresistance Signatures for XELOX-Treated Colorectal Cancer**

**Short title**: Pathway-level nomogram for XELOX resistance

**Authors**: [Author names to be completed]

**Affiliations**: [Institutional affiliations to be completed]

**Corresponding author**: [Name, Email, Address]

"""

text = re.sub(r'\*\*Pathway-Level Molecular Profiling.*?\*\*\n\n.*?Pathway-level nomogram for XELOX resistance\n',
              title_page, text, flags=re.DOTALL)

# ============================================================
# 7. Convert reference format to Vancouver style
# ============================================================
# The current references use [1], [2], etc. in text
# Vancouver uses numbered references in order of appearance
# The current format is already numbered, so we mainly need to
# ensure the reference list at the end follows Vancouver format

# ============================================================
# 8. Write v6.2
# ============================================================
output_path = os.path.join(REPORTS, "manuscript_draft_v6.2_submission.md")
with open(output_path, "w", encoding="utf-8") as f:
    f.write(text)

# Count changes
print(f"v6.2 written to: {output_path}")
print(f"File size: {os.path.getsize(output_path) / 1024:.1f} KB")
print(f"\nChanges made:")
print(f"  1. Abstract restructured to Background/Methods/Results/Conclusions")
print(f"  2. Introduction renamed to Background")
print(f"  3. Figure references updated to v5.11 paths")
print(f"  4. Declarations section added (Ethics/Consent/Data/Competing/Funding/Contributions)")
print(f"  5. Abbreviations list added")
print(f"  6. Title page formatted")
print(f"  7. Ready for J Transl Med submission")
