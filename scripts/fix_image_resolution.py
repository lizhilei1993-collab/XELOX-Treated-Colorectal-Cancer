#!/usr/bin/env python3
"""Fix image DPI and resolution for V5.1 manuscript.
Upscales low-res PNGs to 300 DPI at 6.5 inch target width."""
import sys; sys.stdout.reconfigure(encoding='utf-8')
import os, struct
from PIL import Image

BASE_PNG = r'/path/to/xelox_project\results\figures_png'
OUT_DIR = os.path.join(BASE_PNG, 'hd')  # high-DPI output
TMP_DIR = r'E:/tmp/output/figures_hd'
TARGET_WIDTH_INCH = 6.5
TARGET_DPI = 300
MIN_PIXELS = int(TARGET_WIDTH_INCH * TARGET_DPI)  # 1950

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(TMP_DIR, exist_ok=True)

def fix_png_dpi(input_path, output_path):
    """Upscale to target DPI if needed, then save with correct DPI metadata."""
    with Image.open(input_path) as img:
        w, h = img.size
        dpi = img.info.get('dpi', (96, 96))[0] if 'dpi' in img.info else 96
        
        # Check if image needs upscaling
        if w < MIN_PIXELS:
            scale = MIN_PIXELS / w
            new_w = MIN_PIXELS
            new_h = int(h * scale)
            img = img.resize((new_w, new_h), Image.LANCZOS)
            w, h = new_w, new_h
            action = f'upscaled {scale:.2f}x'
        else:
            action = 'kept native'
        
        # Save with 300 DPI metadata
        img.save(output_path, dpi=(TARGET_DPI, TARGET_DPI), quality=95)
        
        print(f'  {os.path.basename(input_path):35s} {w:5d}x{h:<5d} -> {output_path} [{action}]')

# Process main figures
for f in sorted(os.listdir(BASE_PNG)):
    if not f.endswith('.png'):
        continue
    input_path = os.path.join(BASE_PNG, f)
    output_path = os.path.join(OUT_DIR, f)
    fix_png_dpi(input_path, output_path)

# Process supplementary figure PDFs -> PNG if possible
# Check for PDFs that need conversion
results_fig = r'/path/to/xelox_project\results\figures'
tmp_output_fig = r'E:/tmp/output/figures'

pdf_to_convert = [
    ('fair_comparison', 'fair_comparison_c_index.pdf'),
    ('nomogram', 'bootstrap_coef_forest.pdf'),
    ('tgfb_paradox', 'fig2_gene_cox_forest.png'),  # already PNG
]

for subdir, fname in pdf_to_convert:
    src = os.path.join(results_fig, subdir, fname)
    if os.path.exists(src) and fname.endswith('.png'):
        # Already PNG - just copy
        out = os.path.join(TMP_DIR, f'figS_{subdir}.png')
        fix_png_dpi(src, out)

# Also copy Ridge/EN figures from tmp output
for src_name in ['ridge_cox_path.pdf', 'elastic_net_stability.pdf']:
    src = os.path.join(tmp_output_fig, 'overfitting', src_name)
    if os.path.exists(src):
        out_png = os.path.join(TMP_DIR, src_name.replace('.pdf', '.png'))
        print(f'  [PDF] {src_name} -> cannot convert without pdf2image')

print(f'\nAll HD PNGs saved to: {OUT_DIR}')
print(f'  -> Total files: {len(os.listdir(OUT_DIR))}')
