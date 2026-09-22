#!/usr/bin/env python3
"""
Unify all figure dimensions to consistent size for HTML gallery.
Target: 1950px width for all single-column figures.
"""
import os
from PIL import Image
import numpy as np

BASE = r'/path/to/xelox_project'
SRC = os.path.join(BASE, 'results', 'figures_png', 'hd')
OUT = os.path.join(BASE, 'results', 'figures_png', 'unified')

os.makedirs(OUT, exist_ok=True)

TARGET_WIDTH = 1950
TARGET_DPI = 300

print('=' * 60)
print('Unifying Figure Dimensions')
print('=' * 60)

# Figures that should be resized
figures = [
    'Fig1_pipeline.png',
    'Fig2_SHAP_beeswarm.png',
    'Fig2_SHAP_bar.png',
    'Fig3A_pathway_volcano.png',
    'Fig3B_pathway_heatmap.png',
    'Fig3C_pathway_boxplot.png',
    'Fig4A_PRS_ROC.png',
    'Fig4B_prs_km_curve.png',
    'Fig4C_cox_forest.png',
    'Fig4D_nomogram.png',
    'Fig4E_calibration_12mo.png',
    'Fig4E_calibration_36mo.png',
    'Fig4F_calibration_60mo.png',
    'Fig5_BP_dotplot.png',
    'Fig5_KEGG_barplot.png',
    'Fig6A_ablation_forest.png',
    'Fig6B_method_comparison.png',
]

for fname in figures:
    src_path = os.path.join(SRC, fname)
    out_path = os.path.join(OUT, fname)
    
    if not os.path.exists(src_path):
        print(f'  [SKIP] {fname} not found')
        continue
    
    try:
        img = Image.open(src_path)
        w, h = img.size
        
        # Calculate new height maintaining aspect ratio
        new_w = TARGET_WIDTH
        new_h = int(h * (new_w / w))
        
        # Resize with LANCZOS
        img_resized = img.resize((new_w, new_h), Image.LANCZOS)
        
        # Save with 300 DPI
        img_resized.save(out_path, dpi=(TARGET_DPI, TARGET_DPI))
        
        kb_out = os.path.getsize(out_path) / 1024
        print(f'  {fname:<35} {w}x{h} -> {new_w}x{new_h} ({kb_out:.0f}KB)')
        
    except Exception as e:
        print(f'  [ERROR] {fname}: {e}')

print()
print('=' * 60)
print('Done! Unified figures saved to:', OUT)
print('=' * 60)
