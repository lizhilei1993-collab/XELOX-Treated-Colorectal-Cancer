"""Convert core PDF figures to PNG for Word document embedding."""
import fitz  # PyMuPDF
import os

BASE = r"/path/to/xelox_project\results"
OUT = os.path.join(BASE, "figures_png")
os.makedirs(OUT, exist_ok=True)

DPI = 200  # High resolution for Word embedding

def convert_pdf(pdf_path, out_name, zoom=2.0):
    """Convert first page of PDF to PNG."""
    if not os.path.exists(pdf_path):
        print(f"  [NOT FOUND] {pdf_path}")
        return None
    doc = fitz.open(pdf_path)
    page = doc[0]
    mat = fitz.Matrix(zoom, zoom)
    pix = page.get_pixmap(matrix=mat)
    out_path = os.path.join(OUT, out_name)
    pix.save(out_path)
    doc.close()
    print(f"  [OK] {out_path}")
    return out_path

# === Figure 2: SHAP ===
convert_pdf(
    os.path.join(BASE, "figures", "ml_phase2", "SHAP_beeswarm.pdf"),
    "Fig2_SHAP_beeswarm.png", zoom=2.5
)
convert_pdf(
    os.path.join(BASE, "figures", "ml_phase2", "SHAP_summary_bar_top20.pdf"),
    "Fig2_SHAP_bar.png", zoom=2.5
)

# === Figure 3: Pathway Activity ===
convert_pdf(
    os.path.join(BASE, "figures", "pathway_diff", "pathway_volcano_5panel.pdf"),
    "Fig3A_pathway_volcano.png", zoom=2.0
)
convert_pdf(
    os.path.join(BASE, "figures", "pathway_diff", "pathway_heatmap_consensus.pdf"),
    "Fig3B_pathway_heatmap.png", zoom=2.0
)
convert_pdf(
    os.path.join(BASE, "figures", "pathway_diff", "pathway_boxplot_top4.pdf"),
    "Fig3C_pathway_boxplot.png", zoom=2.0
)

# === Figure 4: PRS ===
convert_pdf(
    os.path.join(BASE, "figures", "lasso_prs_roc_curves.pdf"),
    "Fig4A_PRS_ROC.png", zoom=2.0
)
convert_pdf(
    os.path.join(BASE, "figures", "sigscore_roc_curves.pdf"),
    "Fig4A_sigscore_ROC.png", zoom=2.0
)

# PRS figures from Script 16 (Figure 4B-F)
import shutil
PRS_SRC = os.path.join(BASE, "figures_prs")
os.makedirs(PRS_SRC, exist_ok=True)

prs_files = {
    "prs_km_curve.pdf": "Fig4B_prs_km_curve.png",
    "cox_prs_forest.pdf": "Fig4C_cox_forest.png",
    "nomogram_prs.pdf": "Fig4D_nomogram.png",
    "calibration_12mo_prs.pdf": "Fig4E_calibration_12mo.png",
    "calibration_36mo_prs.pdf": "Fig4E_calibration_36mo.png",
    "calibration_60mo_prs.pdf": "Fig4F_calibration_60mo.png",
}
alt_prs_dir = r"C:\xelox_work\results\figures\prs"
for src_name, dst_name in prs_files.items():
    # Try project dir first, then fallback to /path/to/xelox_project
    src_path = os.path.join(PRS_SRC, src_name)
    if not os.path.exists(src_path):
        alt_path = os.path.join(alt_prs_dir, src_name)
        if os.path.exists(alt_path):
            shutil.copy(alt_path, src_path)
            print(f"  [COPY] {alt_path} -> {PRS_SRC}")
    # Convert
    if os.path.exists(src_path):
        convert_pdf(src_path, dst_name, zoom=2.0)
    else:
        print(f"  [NOT FOUND] {src_name} (checked {PRS_SRC} and {alt_prs_dir})")

# === ML Phase 2 supplementary ===
convert_pdf(
    os.path.join(BASE, "figures", "ml_phase2", "ROC_curves_phase2.pdf"),
    "ML_ROC_phase2.png", zoom=2.0
)
convert_pdf(
    os.path.join(BASE, "figures", "ml_phase2", "Validation_multi_level_auc.pdf"),
    "ML_validation_auc.png", zoom=2.0
)
convert_pdf(
    os.path.join(BASE, "figures", "ml_phase2", "TCGA_validation_roc.pdf"),
    "ML_TCGA_validation.png", zoom=2.0
)
convert_pdf(
    os.path.join(BASE, "figures", "ml_phase2", "SHAP_dependence_plots.pdf"),
    "SHAP_dependence.png", zoom=2.0
)

# === Enrichment (Figure 5) ===
convert_pdf(
    os.path.join(BASE, "figures", "All_BP_dotplot.pdf"),
    "Fig5_BP_dotplot.png", zoom=2.0
)
convert_pdf(
    os.path.join(BASE, "figures", "All_BP_barplot.pdf"),
    "Fig5_BP_barplot.png", zoom=2.0
)
convert_pdf(
    os.path.join(BASE, "figures", "KEGG_All_barplot.pdf"),
    "Fig5_KEGG_barplot.png", zoom=2.0
)
convert_pdf(
    os.path.join(BASE, "figures", "KEGG_HighPriority_barplot.pdf"),
    "Fig5_KEGG_high.png", zoom=2.0
)

# === Meta-analysis supplementary ===
convert_pdf(
    os.path.join(BASE, "figures", "meta_analysis", "volcano_combined.pdf"),
    "Meta_volcano.png", zoom=2.0
)
convert_pdf(
    os.path.join(BASE, "figures", "meta_analysis", "top_genes_forest.pdf"),
    "Meta_forest.png", zoom=2.0
)
convert_pdf(
    os.path.join(BASE, "figures", "meta_analysis", "direction_heatmap.pdf"),
    "Meta_direction_heatmap.png", zoom=2.0
)

# === GSEA ===
convert_pdf(
    os.path.join(BASE, "figures", "gsea", "GSEA_Meta_BP_running_antigen_processing_and_presentation_of_p.pdf"),
    "GSEA_antigen_processing.png", zoom=2.0
)
convert_pdf(
    os.path.join(BASE, "figures", "gsea", "GSEA_KEGG_running_Natural_killer_cell_mediated_cytotoxicit.pdf"),
    "GSEA_NK_cytotoxicity.png", zoom=2.0
)

# === WGCNA (Supplementary) ===
convert_pdf(
    os.path.join(BASE, "figures", "WGCNA_module_trait_heatmap.png").replace("results\\figures\\", "results\\figures\\"),
    "WGCNA_module_trait.png", zoom=2.0
)
# The WGCNA is already PNG, just copy it
import shutil
src = os.path.join(BASE, "figures", "WGCNA_module_trait_heatmap.png")
if os.path.exists(src):
    shutil.copy(src, os.path.join(OUT, "WGCNA_module_trait.png"))
    print(f"  [COPY] {src} -> PNG")

# === Batch correction supp ===
convert_pdf(
    os.path.join(BASE, "figures", "batch_correction", "ComBat_PCA_diagnostic.pdf"),
    "ComBat_PCA.png", zoom=2.0
)

# === Same-platform ROC ===
convert_pdf(
    os.path.join(BASE, "figures", "same_platform_roc_curves.pdf"),
    "Same_platform_ROC.png", zoom=2.0
)

print("\n=== ALL CONVERSIONS COMPLETE ===")
print(f"Output directory: {OUT}")
for f in sorted(os.listdir(OUT)):
    sz = os.path.getsize(os.path.join(OUT, f))
    print(f"  {f}: {sz/1024:.1f} KB")
