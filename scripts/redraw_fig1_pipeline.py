#!/usr/bin/env python3
"""
Redraw Figure 1: Study design and analytical pipeline
Following TRIPOD/REMARK reporting guidelines with explicit sample flow.

Key requirements:
1. Explicit sample counts at each stage
2. Clear exclusion criteria and counts
3. Show分流 of intermediate samples (RFS 12-36 months)
4. Both classification and survival analysis pathways
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np

# Publication settings
DPI = 300
FONT_FAMILY = 'Arial'
plt.rcParams.update({
    'font.family': FONT_FAMILY,
    'font.size': 8,
    'figure.dpi': DPI,
    'savefig.dpi': DPI,
    'savefig.bbox': 'tight',
})

# Color scheme (colorblind-friendly)
COLORS = {
    'discovery': '#E64B35',     # Vermillion (red)
    'validation': '#4DBBD5',    # Teal
    'excluded': '#999999',      # Gray
    'analysis': '#00A087',      # Green
    'box_bg': '#F5F5F5',        # Light gray background
    'arrow': '#333333',         # Dark gray
    'text': '#000000',          # Black
    'highlight': '#FFC000',     # Gold for key numbers
}

def draw_box(ax, x, y, w, h, text, color='white', fontsize=7, bold=False, 
             border_color='#333333', border_width=1, alpha=1.0):
    """Draw a rounded rectangle with text."""
    box = FancyBboxPatch((x, y), w, h, 
                         boxstyle="round,pad=0.02",
                         facecolor=color, edgecolor=border_color,
                         linewidth=border_width, alpha=alpha,
                         transform=ax.transData)
    ax.add_patch(box)
    
    weight = 'bold' if bold else 'normal'
    ax.text(x + w/2, y + h/2, text, ha='center', va='center',
            fontsize=fontsize, fontweight=weight, color=COLORS['text'],
            wrap=True)

def draw_arrow(ax, x1, y1, x2, y2, color='#333333', width=1.5, style='->'):
    """Draw an arrow between two points."""
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=width))

def draw_arrow_label(ax, x1, y1, x2, y2, label, label_pos=0.5, color='#333333'):
    """Draw an arrow with a label."""
    draw_arrow(ax, x1, y1, x2, y2, color)
    # Calculate midpoint
    mx = x1 + (x2 - x1) * label_pos
    my = y1 + (y2 - y1) * label_pos
    ax.text(mx, my, label, ha='center', va='center', fontsize=6,
            color=color, fontstyle='italic')

def create_pipeline():
    """Create the study pipeline figure."""
    fig, ax = plt.subplots(figsize=(7.5, 9))
    ax.set_xlim(0, 7.5)
    ax.set_ylim(0, 9)
    ax.axis('off')
    
    # Title
    ax.text(3.75, 8.7, 'Figure 1. Study Design and Analytical Pipeline',
            ha='center', va='center', fontsize=11, fontweight='bold')
    ax.text(3.75, 8.5, 'Following TRIPOD/REMARK reporting guidelines',
            ha='center', va='center', fontsize=8, color='gray', fontstyle='italic')
    
    # ======== LEFT COLUMN: Discovery Cohorts ========
    col1_x = 0.3
    col2_x = 2.7
    col3_x = 5.1
    
    # Row 1: GEO Data Retrieval
    y = 8.0
    draw_box(ax, col1_x, y, 2.2, 0.5, 
             'GEO Data Retrieval\n5 Cohorts', 
             color=COLORS['discovery'], fontsize=8, bold=True)
    
    # Row 2: Cohort Details
    y = 7.0
    cohorts_text = (
        'GSE39582 (n=585)\n'
        'GSE104645 (n=193)\n'
        'GSE28702 (n=83)\n'
        'GSE72970 (n=124)\n'
        'GSE69657 (n=30)\n'
        'Total: 1,015'
    )
    draw_box(ax, col1_x, y, 2.2, 0.9, cohorts_text, fontsize=6.5)
    
    # Row 3: Quality Control
    y = 5.6
    draw_box(ax, col1_x, y, 2.2, 0.5, 
             'Quality Control\n& Preprocessing',
             color='#FFE4B5', fontsize=7, bold=True)
    
    # QC exclusion arrow
    draw_arrow_label(ax, col1_x + 2.3, y + 0.25, col1_x + 2.8, y + 0.25,
                     'Excluded:\n• Missing expression\n• Non-CRC samples\n• No chemotherapy info',
                     label_pos=0.5, color=COLORS['excluded'])
    
    # Row 4: XELOX Subgroup
    y = 4.6
    draw_box(ax, col1_x, y, 2.2, 0.5,
             'XELOX/FOLFOX-treated\nmCRC Subgroup',
             color=COLORS['discovery'], fontsize=7, bold=True)
    
    # Exclusion count
    ax.text(col1_x + 2.3, y + 0.5, 'n = 476 (47%)\nExcluded: 539',
            fontsize=6, ha='left', color=COLORS['excluded'], fontstyle='italic')
    
    # ======== CENTER COLUMN: GSE39582 Focus ========
    
    # Row 5: GSE39582 Focus
    y = 3.5
    draw_box(ax, col2_x, y, 2.2, 0.5,
             'GSE39582 XELOX\nDiscovery Cohort\n(n = 239)',
             color=COLORS['highlight'], fontsize=7, bold=True, border_width=2)
    
    # Arrow from left
    draw_arrow(ax, col1_x + 2.2, y + 0.25, col2_x, y + 0.25, COLORS['arrow'])
    ax.text(col1_x + 2.35, y + 0.35, 'Primary\nAnalysis', fontsize=6, ha='center')
    
    # ======== RIGHT COLUMN: Two Analysis Pathways ========
    
    # Path A: Binary Classification (top)
    y_class = 3.8
    draw_box(ax, col3_x, y_class, 2.0, 0.7,
             'Binary Classification\n(Resistant vs Sensitive)',
             color=COLORS['analysis'], fontsize=7, bold=True)
    
    # Arrow to classification
    draw_arrow(ax, col2_x + 2.2, y_class + 0.35, col3_x, y_class + 0.35, COLORS['arrow'])
    
    # Classification sample flow
    y_class_detail = 2.8
    class_text = (
        'GSE39582 (n=164)\n'
        '• Resistant: 45\n'
        '• Sensitive: 119\n'
        '• Intermediate: Excluded (n=75)\n\n'
        'XGBoost → SHAP → 10-gene\n'
        'fingerprint'
    )
    draw_box(ax, col3_x, y_class_detail, 2.0, 0.9, class_text, fontsize=6)
    
    # Path B: Survival Analysis (bottom)
    y_surv = 1.5
    draw_box(ax, col3_x, y_surv, 2.0, 0.7,
             'Survival Analysis\n(RFS Endpoint)',
             color=COLORS['analysis'], fontsize=7, bold=True)
    
    # Arrow to survival
    draw_arrow(ax, col2_x + 2.2, y_surv + 0.35, col3_x, y_surv + 0.35, COLORS['arrow'])
    
    # Survival sample flow
    y_surv_detail = 0.5
    surv_text = (
        'GSE39582 (n=227)\n'
        '• Events: 79 (34.8%)\n'
        '• Intermediate RFS: Included\n'
        '  (n=66, 12-36 months)\n\n'
        'ssGSEA → AIC-Cox → Nomogram'
    )
    draw_box(ax, col3_x, y_surv_detail, 2.0, 0.9, surv_text, fontsize=6)
    
    # ======== SAMPLE FLOW ANNOTATIONS ========
    
    # GSE39582 sample split annotation
    y_note = 2.4
    ax.text(col2_x + 0.1, y_note,
            'GSE39582 Sample Split:\n'
            '  Intermediate RFS (n=66):\n'
            '    Excluded from classification\n'
            '    Included in survival analysis\n'
            '  Classification: n=164 (R vs S)',
            fontsize=6, va='top', ha='left',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    # ======== EXTERNAL VALIDATION ========
    
    # Row 6: External Validation
    y = -0.8
    draw_box(ax, 1.5, y, 5.0, 0.5,
             'External Validation Cohorts',
             color=COLORS['validation'], fontsize=8, bold=True)
    
    # Validation details
    y = -1.5
    val_text = (
        'TCGA-COAD/READ (oxaliplatin-treated, n=140)  •  GSE83129 (FOLFOX, n=26)  •  '
        'GSE104645 (oxaliplatin, n=113)  •  GSE72970 (FOLFOX, n=32)'
    )
    draw_box(ax, 0.5, y, 7.0, 0.4, val_text, fontsize=6)
    
    # Arrow from validation cohorts to analysis
    draw_arrow(ax, 3.75, -0.3, 3.75, -0.8, COLORS['arrow'])
    ax.text(3.75, -0.55, 'External validation of 10-gene fingerprint\nand nomogram performance',
            fontsize=6, ha='center', va='center', fontstyle='italic')
    
    # ======== SAMPLE FLOW SUMMARY BOX ========
    
    # Summary box
    y = 0.8
    summary_text = (
        '┌─────────────────────────────────────────┐\n'
        '│  Sample Flow Summary (TRIPOD)           │\n'
        '├─────────────────────────────────────────┤\n'
        '│  Discovery pool: 1,015                  │\n'
        '│  After QC: 476 (47%)                    │\n'
        '│  GSE39582 XELOX: 239                    │\n'
        '│  ├─ Classification: 164 (45R + 119S)    │\n'
        '│  │  Excluded: 75 intermediate           │\n'
        '│  └─ Survival: 227 (79 events)           │\n'
        '│     RFS 12-36mo: 66 (included)          │\n'
        '└─────────────────────────────────────────┘'
    )
    # Use a box instead
    draw_box(ax, 0.3, 0.2, 2.3, 1.5, 
             'Sample Flow Summary\n\n'
             'Discovery pool: 1,015\n'
             'After QC: 476 (47%)\n'
             'GSE39582: 239\n'
             '  Classification: 164\n'
             '  Survival: 227 (79 events)',
             fontsize=6, color='#F0F0F0', border_color='#666666')
    
    plt.tight_layout()
    return fig


if __name__ == '__main__':
    print('Creating Figure 1: Study Pipeline...')
    
    fig = create_pipeline()
    
    # Save to both locations
    import os
    base = r'/path/to/xelox_project'
    hd_path = os.path.join(base, 'results', 'figures_png', 'hd', 'Fig1_pipeline.png')
    fig_path = os.path.join(base, 'results', 'figures_png', 'Fig1_pipeline.png')
    
    fig.savefig(hd_path, dpi=DPI, bbox_inches='tight', facecolor='white')
    fig.savefig(fig_path, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    
    print(f'  Saved HD: {hd_path}')
    print(f'  Saved: {fig_path}')
    print('Done!')
