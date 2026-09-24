"""Generate v6.0 supplementary materials docx with clean tables."""
import os
from docx import Document
from docx.shared import Inches, Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn

PROJECT = r"/path/to/xelox_project"
OUT_PATH = os.path.join(PROJECT, "reports", "supplementary_materials_v6.0.docx")

doc = Document()

# Page setup
for section in doc.sections:
    section.page_width = Cm(21.0)
    section.page_height = Cm(29.7)
    section.top_margin = Cm(2.0)
    section.bottom_margin = Cm(2.0)
    section.left_margin = Cm(2.0)
    section.right_margin = Cm(2.0)

style = doc.styles['Normal']
font = style.font
font.name = 'Times New Roman'
font.size = Pt(10)
style.paragraph_format.line_spacing = 1.15
style.paragraph_format.space_after = Pt(4)

BLUE = "2E75B6"
LIGHT_BG = "F0F4F8"
BORDER_GRAY = "CCCCCC"

def add_heading(text, level=1):
    h = doc.add_heading(text, level=level)
    for run in h.runs:
        run.font.name = 'Times New Roman'
        run.font.color.rgb = RGBColor(0, 0, 0)
    return h

def add_para(text, bold=False, size=10, align=None, italic=False):
    p = doc.add_paragraph()
    run = p.add_run(text)
    run.font.name = 'Times New Roman'
    run.font.size = Pt(size)
    run.bold = bold
    run.italic = italic
    if align:
        p.alignment = align
    return p

def make_table(headers, rows, col_widths=None):
    """Create a styled table with blue header."""
    ncols = len(headers)
    table = doc.add_table(rows=len(rows)+1, cols=ncols)
    table.style = 'Table Grid'
    
    # Style all cells
    for r_idx in range(len(rows)+1):
        for c_idx in range(ncols):
            cell = table.rows[r_idx].cells[c_idx]
            for p in cell.paragraphs:
                p.paragraph_format.space_after = Pt(2)
                p.paragraph_format.space_before = Pt(2)
                for run in p.runs:
                    run.font.name = 'Times New Roman'
                    run.font.size = Pt(9)
    
    # Header
    for j, h in enumerate(headers):
        cell = table.rows[0].cells[j]
        cell.text = ''
        p = cell.paragraphs[0]
        run = p.add_run(h)
        run.font.name = 'Times New Roman'
        run.font.size = Pt(9)
        run.bold = True
        run.font.color.rgb = RGBColor(255, 255, 255)
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        # Blue background
        shading = cell._element.get_or_add_tcPr()
        shd = shading.makeelement(qn('w:shd'), {
            qn('w:fill'): BLUE,
            qn('w:val'): 'clear'
        })
        shading.append(shd)
    
    # Data rows with alternating color
    for r_idx, row_data in enumerate(rows):
        for c_idx, cell_text in enumerate(row_data):
            if c_idx < ncols:
                cell = table.rows[r_idx+1].cells[c_idx]
                cell.text = ''
                p = cell.paragraphs[0]
                run = p.add_run(str(cell_text))
                run.font.name = 'Times New Roman'
                run.font.size = Pt(9)
                if r_idx % 2 == 0:
                    shading = cell._element.get_or_add_tcPr()
                    shd = shading.makeelement(qn('w:shd'), {
                        qn('w:fill'): LIGHT_BG,
                        qn('w:val'): 'clear'
                    })
                    shading.append(shd)
    
    doc.add_paragraph()
    return table

# ============================================================
# Title
# ============================================================
add_heading("Supplementary Materials", 1)
add_para("From Gene-Level Fingerprints to Pathway-Level Prognostic Modeling", 
         bold=False, size=12, italic=True, align=WD_ALIGN_PARAGRAPH.CENTER)
add_para("An Explainable Machine Learning Framework for XELOX Chemoresistance in Colorectal Cancer",
         bold=False, size=11, italic=True, align=WD_ALIGN_PARAGRAPH.CENTER)
add_para("")

add_para("Manuscript reference: Li Z, et al. (2026)", bold=True, size=11)
add_para("Version: V6.0 (streamlined), 2026-06-02", size=11)
add_para("Supplementary Tables: S1--S13 (13 tables)", size=11)
add_para("Supplementary Figures: S1--S16 (16 figures)", size=11)
add_para("")

# ============================================================
# Table S1-S5
# ============================================================
add_heading("S1--S5: Differential Expression & Meta-Analysis", 2)

make_table(
    ["ID", "Title", "Source File", "Key Content"],
    [
        ["S1", "253 Significant Genes from Multi-Cohort Meta-Analysis",
         "results/tables/meta_analysis/deg_meta_filtered.csv (2.1 MB)",
         "Stouffer meta-analysis across 5 GEO cohorts (FDR < 0.05). Direction matrix appended."],
        ["S2", "Five-Cohort DEG Summary",
         "results/tables/DEG_GSE39582_limma.csv etc. (5 files)",
         "Per-cohort top DEGs. Full tables as external data (see Data Availability)."],
        ["S3", "WGCNA Hub Genes",
         "results/tables/WGCNA_hub_genes.csv (27 KB)",
         "Hub genes from co-expression modules (beta=6, scale-free R^2 > 0.85)."],
        ["S4", "Pathway Univariate Cox + BH + Bootstrap Frequency",
         "results/tables/nomogram/univariate_cox_screening_bh.csv; bootstrap_selection_freq.csv",
         "44 pathways: BH-adjusted Cox (7 at FDR < 0.10). LASSO-PRS bootstrap (B=1,000)."],
        ["S5", "Pathway Consistency Matrix (5 Cohorts)",
         "results/tables/pathway_activity/pathway_consistency_matrix.csv (10 KB)",
         "Directional consistency of 44 pathways across all 5 GEO cohorts."],
    ]
)

# ============================================================
# Table S6-S9
# ============================================================
add_heading("S6--S9: Ablation, Diagnostics & External Validation", 2)

make_table(
    ["ID", "Title", "Source File", "Key Content"],
    [
        ["S6", "23-Gene AIC-Selected Cox Model",
         "results/tables/gene_level/gene_cox_aic_model.csv (1 KB)",
         "Backward AIC: coefficients, SE, p. C=0.834 (EPV=3.4)."],
        ["S7", "Ablation Comparison: Feature x Method",
         "results/tables/gene_level/ablation_comparison.csv",
         "C-index: gene XGBoost, pathway LASSO-PRS, gene AIC Cox, pathway AIC nomogram."],
        ["S8", "Nomogram Diagnostics (7-in-1)",
         "results/tables/nomogram_sensitivity/*.csv; calibration_stats.csv",
         "PH test (global p=0.687), dfbeta, LOVO, stratified C, cutoff, bootstrap CI, calibration."],
        ["S9", "TCGA External Validation",
         "results/tables/tcga_validation/tcga_cox_results.csv",
         "Pathway PRS → OS in TCGA oxaliplatin subset (n=140, 11 events)."],
    ]
)

# ============================================================
# Table S10-S13
# ============================================================
add_heading("S10--S13: Late-Stage Refinement & Clinical Utility", 2)

make_table(
    ["ID", "Title", "Source File", "Key Content"],
    [
        ["S10", "Fair Comparison: Gene vs Pathway LASSO-Cox",
         "/tmp/scratch/output/tables/fair_comparison/fair_comparison_summary.csv",
         "C-corrected: 0.406 vs 0.489. Optimism: 0.0945 vs 0.011. EPV: 1.1 vs 11.3."],
        ["S11", "AIC Bootstrap Selection + Coef Distribution",
         "/tmp/scratch/output/tables/nomogram/bootstrap_selection_freq.csv; bootstrap_coef_distribution.csv",
         "Top-3: WNT (66.5%), TGF-beta (59.0%), ECM (49.0%). 7-pathway coef forest."],
        ["S12", "Drug Sensitivity Analysis (oncoPredict)",
         "/tmp/scratch/output/tables/oncoPredict/drug_sensitivity_summary.csv",
         "5-FU p<0.0001; capecitabine p<0.0001; oxaliplatin p=0.600 (n.s.)."],
        ["S13", "DCA Net Benefit & TDROC AUC",
         "results/tables/dca/dca_net_benefit.csv; tdroc_auc.csv",
         "C=0.642. AUC: 12M=0.630, 36M=0.686, 60M=0.652. Net benefit ~5-40%."],
    ]
)

# ============================================================
# Supplementary Figures
# ============================================================
add_heading("Supplementary Figures", 1)

make_table(
    ["ID", "Title", "File(s)", "Format", "Description"],
    [
        ["S1", "WGCNA Analysis", "WGCNA_scale_independence.png, WGCNA_module_trait_heatmap.png", "PNG",
         "Scale-free fit (beta=6, R^2>0.85) + module-trait heatmap. GSE39582 n=164."],
        ["S2", "ComBat Batch Correction PCA", "ComBat_PCA_diagnostic.pdf", "PDF",
         "PCA before/after ComBat across GPL570 cohorts."],
        ["S3", "Pathway-Level Analysis (4-panel)", "pathway_heatmap_consensus.pdf, pathway_boxplot_top4.pdf, lasso_prs_roc_curves.pdf", "PDF",
         "Consensus heatmap, volcano, boxplots, LASSO-PRS ROC across 5 cohorts."],
        ["S4", "Bootstrap Stability (B=1,000)", "bootstrap_selection_freq_bar.pdf, bootstrap_coef_heatmap.pdf, bootstrap_stability_curve.pdf", "PDF",
         "Selection freq, coefficient heatmap, stability curve. Jaccard=0.10."],
        ["S5", "Nomogram Sensitivity (7-panel)", "schoenfeld_residuals.pdf, dfbeta_index_plot.pdf, leave_one_variable_out_cindex.pdf, stratified_cindex.pdf, cutoff_sweep.pdf, cutoff_cindex_sweep.pdf, bootstrap_cindex_distribution.pdf", "PDF",
         "PH test, dfbeta, LOVO C-index, stratified, cutoff sweep, bootstrap CI."],
        ["S6", "Gene-Level Cox Forest (23 genes)", "forest_plot_gene_cox.pdf", "PDF",
         "AIC-selected 23-gene Cox forest plot. Figure 6A source."],
        ["S7", "TCGA Validation (4-panel)", "tcga_roc_3yr.pdf, tcga_roc_5yr.pdf, tcga_km_os.pdf, tcga_prs_boxplot.pdf; tdroc_curve.pdf", "PDF",
         "ROC at 36/60mo, OS KM, risk score; IPCW TDROC at 12/36/60mo. TCGA n=140."],
        ["S8", "Meta-Analysis", "Meta_volcano.png, Meta_forest.png, Meta_direction_heatmap.png", "PNG",
         "Volcano, top genes forest, direction heatmap. 253 genes, FDR<0.05."],
        ["S9", "TGF-beta Paradox (5-panel)", "fig1_gene_overlap.png ~ fig5_direction_proportion.png", "PNG",
         "Venn, Cox forest, expression direction, score scatter, direction proportion."],
        ["S10", "Fair Comparison: Gene vs Pathway C-index", "fair_comparison_c_index.pdf", "PDF",
         "C-corrected bar plot. Gene 0.406 vs Pathway 0.489. B=200."],
        ["S11", "Bootstrap Coefficient Forest (7 Pathways)", "bootstrap_coef_forest.pdf", "PDF",
         "Original coef, bootstrap median, 95% CI. WNT (100%) and TGF-beta (98.5%) dominant."],
        ["S12", "Drug Sensitivity Boxplots (3-Drug)", "oxaliplatin/fluorouracil/capecitabine_sensitivity_boxplot.pdf", "PDF",
         "High vs Low risk. 5-FU/capecitabine p<0.0001; oxaliplatin p=0.600."],
    ]
)

# Remaining figures S13-S16
make_table(
    ["ID", "Title", "File(s)", "Format", "Description"],
    [
        ["S13", "Calibration Curve (12/36/60mo)", "calibration_curve.pdf", "PDF",
         "Bootstrap-based calibration. Slope=0.817. Brier: 12M=0.142, 36M=0.201, 60M=0.218."],
        ["S14", "Decision Curve Analysis", "dca_curve.pdf", "PDF",
         "Net benefit vs threshold (0-50%). Net benefit across ~5-40%."],
        ["S15", "Method Comparison (Ablation)", "Fig6B_method_comparison.png", "PNG",
         "C-index: gene/pathway x ML/Cox. Pathway aggregation as implicit regularization."],
        ["S16", "Open Targets Drug-Gene Network", "generated from Open Targets Platform", "N/A",
         "Drug-gene interactions for 253 genes. FLT1: 48 compounds. Most CRC drug targets absent."],
    ]
)

# ============================================================
# Count Summary
# ============================================================
add_heading("Count Summary", 2)

make_table(
    ["Resource", "Count"],
    [
        ["Supplementary Tables", "13 (S1--S13)"],
        ["Supplementary Figures", "16 (S1--S16)"],
        ["Main Text Figures", "6 (Fig 1--Fig 6)"],
        ["Main Text Tables", "10 (Table 1--10)"],
        ["Discovery GEO cohorts", "5"],
        ["External validation", "TCGA (n = 140)"],
        ["Total cohort N", "585 CRC patients"],
    ]
)

# ============================================================
# References
# ============================================================
add_heading("GEO Datasets", 2)
refs = [
    "1. GSE39582 — Marisa L, et al. PLoS Med. 2013;10(5):e1001453. (XELOX n=164)",
    "2. GSE104645 — Nishioka Y, et al. Cancer Sci. 2018. (oxaliplatin n=113)",
    "3. GSE28702 — Watanabe T, et al. Oncotarget. 2012. (oxaliplatin n=56)",
    "4. GSE72970 — Del Rio M, et al. Eur J Cancer. 2017;76:68-75. (FOLFOX/FOLFIRI n=72)",
    "5. GSE69657 — Tovar I, et al. Cancer Res. 2015. (FOLFOX4 n=71)",
    "6. TCGA-COAD/READ — Cancer Genome Atlas Network. Nature. 2012;487(7407):330-7. (oxaliplatin n=140)",
]
for ref in refs:
    add_para(ref, size=9)

add_heading("Software", 2)
add_para("R 4.2+ with limma, WGCNA, sva, GSVA, glmnet, rms, survival, survminer, pROC", size=9)
add_para("Python 3.10+ with scikit-learn, xgboost, lightgbm, shap, pandas, numpy", size=9)

doc.save(OUT_PATH)
print(f"Supplementary DOCX saved: {OUT_PATH}")
print(f"Size: {os.path.getsize(OUT_PATH)/1024:.0f} KB")
