#!/usr/bin/env python3
"""
Figure 2: SHAP Analysis from XGBoost model data.
Panel A: Beeswarm plot (per-sample SHAP values)
Panel B: Mean |SHAP| bar plot (feature importance)

Uses XGBoost's built-in predcontrib for SHAP values (no shap library needed).
Loads training data from RDS using pandas (via rpy2/pyreadr) or CSV fallback.
"""
import os, sys, warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

# Paths
BASE = r'/path/to/xelox_project'
MODEL_PATH = os.path.join(BASE, 'results/tables/ml_phase2/XGBoost_model_phase2.xgb')
OUT_DIR = os.path.join(BASE, 'results/figures_png/v5.11')
SHAP_CSV = os.path.join(BASE, 'results/tables/ml_phase2/SHAP_values_per_sample.csv')
TRAIN_RDS = os.path.join(BASE, 'results/tables/ml_phase2/XELOX_training_data_prepared.rds')

# Settings
DPI = 300
FAMILY = 'Arial'
W_INCH = 3.35  # single column = 85mm

plt.rcParams.update({
    'font.family': FAMILY, 'font.size': 9,
    'axes.labelsize': 9, 'axes.titlesize': 0,
    'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'legend.fontsize': 7, 'figure.dpi': DPI,
    'savefig.dpi': DPI, 'figure.facecolor': 'white',
    'savefig.facecolor': 'white', 'savefig.bbox': 'tight',
})

# ── Step 1: Try to compute SHAP values from model ──
def compute_shap_from_model():
    """Compute per-sample SHAP values using XGBoost predcontrib."""
    import xgboost as xgb
    
    # Load model
    model = xgb.Booster()
    model.load_model(MODEL_PATH)
    fnames = model.feature_names
    print(f'Model features: {fnames}')
    
    # Try to load training data from RDS using pyreadr
    try:
        import pyreadr
        result = pyreadr.read_r(TRAIN_RDS)
        df = list(result.values())[0]
        print(f'Training data loaded: {df.shape}')
    except Exception as e:
        print(f'Cannot load RDS: {e}')
        # Try to reconstruct training data from gene expression
        # Fall back to SHAP summary only
        return None, fnames
    
    # Ensure feature columns match
    available = [f for f in fnames if f in df.columns]
    if len(available) < len(fnames):
        print(f'Warning: only {len(available)}/{len(fnames)} features available')
        # Try to find matching columns
        matched = {}
        for feat in fnames:
            if feat in df.columns:
                matched[feat] = feat
            else:
                # Try alternative names
                for col in df.columns:
                    if feat.replace(' /// ', '') in col or col.replace(' /// ', '') == feat:
                        matched[feat] = col
                        break
        available = list(matched.values())
    
    if len(available) == 0:
        print('ERROR: No feature columns found!')
        return None, fnames
    
    # Prepare data matrix
    x = df[[matched.get(f, f) for f in fnames if matched.get(f, f) in df.columns]].values
    valid_fnames = [f for f in fnames if matched.get(f, f) in [matched.get(f2, f2) for f2 in fnames if matched.get(f2, f2) in df.columns]]
    
    # Create DMatrix
    dmat = xgb.DMatrix(x, feature_names=valid_fnames)
    
    # Compute SHAP values (predcontrib gives marginal contribution)
    shap_contrib = model.predict(dmat, predcontrib=True)  # n_samples x (n_features + 1)
    
    # Last column is bias
    bias = shap_contrib[:, -1]
    shap_vals = shap_contrib[:, :-1]  # Remove bias column
    print(f'SHAP matrix: {shap_vals.shape}')
    
    # Also get labels for coloring
    labels = None
    for col in ['group', 'status', 'label', 'target']:
        if col in df.columns:
            labels = df[col].values
            break
    
    return shap_vals, valid_fnames, labels


# ── Step 2: Check for pre-computed SHAP CSV ──
def load_existing_shap():
    """Load pre-computed SHAP values if available."""
    # Check if we have per-sample SHAP CSV
    shap_csv = os.path.join(BASE, 'results/tables/ml_phase2/SHAP_values_per_sample.csv')
    if os.path.exists(shap_csv):
        df = pd.read_csv(shap_csv)
        print(f'Loaded existing SHAP values: {df.shape}')
        return df
    
    # Check for SHAP summary
    summary_csv = os.path.join(BASE, 'results/tables/ml_phase2/SHAP_summary_phase2.csv')
    if os.path.exists(summary_csv):
        print(f'Only SHAP summary available: {summary_csv}')
        return None
    return None


# ── Step 3: Generate figures ──
def create_figure2():
    print('\n=== Creating Figure 2: SHAP Analysis ===')
    
    # Try to compute SHAP values
    result = compute_shap_from_model()
    
    if result is not None and result[0] is not None:
        shap_vals, fnames, labels = result
        
        # Panel A: Beeswarm plot
        print('Creating beeswarm plot...')
        fig, ax = plt.subplots(figsize=(W_INCH, 3.5))
        
        # Sort features by mean |SHAP|
        mean_abs = np.mean(np.abs(shap_vals), axis=0)
        sorted_idx = np.argsort(mean_abs)[::-1]
        sorted_fnames = [fnames[i] for i in sorted_idx]
        
        # Create beeswarm
        positions = []
        for feat_idx, fi in enumerate(sorted_idx):
            vals = shap_vals[:, fi]
            jitter = np.random.normal(0, 0.15, size=len(vals))
            positions.extend([feat_idx] * len(vals))
            
            # Color by value (blue=negative, red=positive)
            for v, j in zip(vals, jitter):
                color = '#C0392B' if v > 0 else '#3498DB'
                ax.plot(v, feat_idx + j, 'o', color=color, markersize=2.5, alpha=0.4, markeredgecolor='none')
        
        ax.set_yticks(range(len(sorted_fnames)))
        ax.set_yticklabels(sorted_fnames, fontsize=7)
        ax.set_xlabel('SHAP value (impact on model output)', fontsize=9)
        ax.axvline(0, color='gray', linestyle='--', linewidth=0.5)
        ax.invert_yaxis()
        
        # Legend
        ax.plot([], [], 'o', color='#C0392B', label='High risk', markersize=6)
        ax.plot([], [], 'o', color='#3498DB', label='Low risk', markersize=6)
        ax.legend(fontsize=6, loc='lower right', framealpha=0.8)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        fig.text(-0.02, 1.06, 'A', fontsize=13, fontweight='bold', transform=ax.transAxes, va='top')
        plt.tight_layout()
        
        out_a = os.path.join(OUT_DIR, 'Fig2A_SHAP_beeswarm.png')
        fig.savefig(out_a, dpi=DPI, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print(f'  Saved: {out_a} ({os.path.getsize(out_a)//1024} KB)')
        
        # Panel B: Mean |SHAP| bar plot
        print('Creating bar plot...')
        fig, ax = plt.subplots(figsize=(W_INCH, 3.5))
        
        # Compute stats
        mean_abs = np.mean(np.abs(shap_vals), axis=0)
        std_abs = np.std(np.abs(shap_vals), axis=0)
        sorted_idx = np.argsort(mean_abs)[::-1]
        
        y_pos = range(len(sorted_idx))
        vals_sorted = mean_abs[sorted_idx]
        stds_sorted = std_abs[sorted_idx]
        names_sorted = [fnames[i] for i in sorted_idx]
        
        # Color top 3 differently
        colors = ['#C0392B' if i < 3 else '#4DBBD5' for i in range(len(sorted_idx))]
        bars = ax.barh(y_pos, vals_sorted, xerr=stds_sorted, color=colors, 
                       edgecolor='white', linewidth=0.5, capsize=3,
                       error_kw={'linewidth': 0.8, 'color': '#666'})
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels(names_sorted, fontsize=7)
        ax.set_xlabel('Mean |SHAP value| ± SD', fontsize=9)
        ax.invert_yaxis()
        
        # Add value labels
        for i, v in enumerate(vals_sorted):
            ax.text(v + 0.01, i, f'{v:.3f}', va='center', fontsize=6.5, color='#333')
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        fig.text(-0.02, 1.06, 'B', fontsize=13, fontweight='bold', transform=ax.transAxes, va='top')
        plt.tight_layout()
        
        out_b = os.path.join(OUT_DIR, 'Fig2B_SHAP_bar.png')
        fig.savefig(out_b, dpi=DPI, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print(f'  Saved: {out_b} ({os.path.getsize(out_b)//1024} KB)')
        
        # Combined figure
        print('Creating combined figure...')
        fig, axes = plt.subplots(1, 2, figsize=(6.7, 3.5))
        plt.subplots_adjust(wspace=0.3)
        
        for ax_idx, (ax, vals, names) in enumerate(zip(axes,
            [sorted_idx, sorted_idx],
            [sorted_fnames, sorted_fnames]
        )):
            if ax_idx == 0:  # Beeswarm
                for feat_idx, fi in enumerate(vals):
                    sv = shap_vals[:, fi]
                    jitter = np.random.normal(0, 0.15, size=len(sv))
                    for v, j in zip(sv, jitter):
                        color = '#C0392B' if v > 0 else '#3498DB'
                        ax.plot(v, feat_idx + j, 'o', color=color, markersize=2, alpha=0.35, markeredgecolor='none')
                ax.set_yticks(range(len(names)))
                ax.set_yticklabels(names, fontsize=7)
                ax.set_xlabel('SHAP value', fontsize=9)
                ax.axvline(0, color='gray', linestyle='--', linewidth=0.5)
                ax.invert_yaxis()
                ax.plot([], [], 'o', color='#C0392B', label='Resistant', markersize=5)
                ax.plot([], [], 'o', color='#3498DB', label='Sensitive', markersize=5)
                ax.legend(fontsize=6, loc='lower right', framealpha=0.8)
            else:  # Bar
                mean_abs_v = np.mean(np.abs(shap_vals), axis=0)
                std_abs_v = np.std(np.abs(shap_vals), axis=0)
                sorted_i = np.argsort(mean_abs_v)[::-1]
                colors = ['#C0392B' if i < 3 else '#4DBBD5' for i in range(len(sorted_i))]
                ax.barh(range(len(sorted_i)), mean_abs_v[sorted_i], xerr=std_abs_v[sorted_i],
                       color=colors, edgecolor='white', linewidth=0.5, capsize=3,
                       error_kw={'linewidth': 0.8, 'color': '#666'})
                ax.set_yticks(range(len(sorted_i)))
                ax.set_yticklabels([fnames[i] for i in sorted_i], fontsize=7)
                ax.set_xlabel('Mean |SHAP| ± SD', fontsize=9)
                ax.invert_yaxis()
                for i, v in enumerate(mean_abs_v[sorted_i]):
                    ax.text(v + 0.01, i, f'{v:.3f}', va='center', fontsize=6.5, color='#333')
            
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_linewidth(0.5)
            ax.spines['bottom'].set_linewidth(0.5)
            fig.text(0.02 if ax_idx == 0 else 0.52, 1.06, chr(65 + ax_idx), fontsize=13, fontweight='bold', transform=ax.transAxes, va='top')
        
        out_combined = os.path.join(OUT_DIR, 'Fig2_SHAP.png')
        fig.savefig(out_combined, dpi=DPI, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print(f'  Saved: {out_combined} ({os.path.getsize(out_combined)//1024} KB)')
        
    else:
        print('WARNING: Cannot compute per-sample SHAP values (RDS data not accessible)')
        print('Using SHAP summary data only for bar plot')
        
        # Load summary
        summary = pd.read_csv(os.path.join(BASE, 'results/tables/ml_phase2/SHAP_summary_phase2.csv'))
        summary = summary.sort_values('MeanAbsSHAP', ascending=True)
        
        fig, ax = plt.subplots(figsize=(W_INCH, 3.5))
        colors = ['#C0392B' if i >= len(summary) - 3 else '#4DBBD5' for i in range(len(summary))]
        ax.barh(range(len(summary)), summary['MeanAbsSHAP'], color=colors, edgecolor='white')
        ax.set_yticks(range(len(summary)))
        ax.set_yticklabels(summary['Gene'].str.replace(' /// ', '/'), fontsize=7)
        ax.set_xlabel('Mean |SHAP| Value', fontsize=9)
        ax.invert_yaxis()
        for i, v in enumerate(summary['MeanAbsSHAP']):
            ax.text(v + 0.01, i, f'{v:.3f}', va='center', fontsize=6.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        fig.text(-0.02, 1.06, 'B', fontsize=13, fontweight='bold', transform=ax.transAxes, va='top')
        plt.tight_layout()
        
        out = os.path.join(OUT_DIR, 'Fig2_SHAP.png')
        fig.savefig(out, dpi=DPI, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print(f'  Saved: {out} ({os.path.getsize(out)//1024} KB)')


if __name__ == '__main__':
    create_figure2()
    print('\nDone!')
