"""
Script 21: PRS Multi-Cohort Meta-Analysis
==========================================
Integrate PRS classification AUC across validation cohorts
using random-effects meta-analysis (DerSimonian-Laird).

Outputs:
  - Forest plot (meta-analysis)
  - Funnel plot (publication bias)
  - Subgroup analysis (XELOX vs other oxaliplatin)
  - Numerical results table (pooled AUC + I²)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import os, sys, json
from math import log, exp, sqrt

# ============================================================
# 0. Setup paths
# ============================================================
PROJECT_ROOT = r"/path/to/xelox_project"
XELOX_WORK = r"C:\xelox_work"
RESULTS_DIR = os.path.join(XELOX_WORK, "results", "meta_prs")
os.makedirs(RESULTS_DIR, exist_ok=True)

# ============================================================
# 1. Input data — from PRS_validation_summary.csv
# ============================================================
cohorts = {
    "Training (merged)": {"n": 277, "auc": 0.639, "ci_l": 0.571, "ci_u": 0.706, "group": "Training"},
    "GSE39582":          {"n": 164, "auc": 0.644, "ci_l": 0.552, "ci_u": 0.735, "group": "XELOX-like"},
    "GSE104645":         {"n": 113, "auc": 0.574, "ci_l": 0.465, "ci_u": 0.682, "group": "Oxaliplatin"},
    "GSE28702":          {"n": 83,  "auc": 0.544, "ci_l": 0.418, "ci_u": 0.670, "group": "mFOLFOX6"},
    "GSE69657":          {"n": 30,  "auc": 0.624, "ci_l": 0.404, "ci_u": 0.845, "group": "FOLFOX4/XELOX"},
    "GSE72970":          {"n": 36,  "auc": 0.311, "ci_l": 0.167, "ci_u": 0.455, "group": "FOLFOX"},
}

# LOCO CV results (from Script 22 — held-out validation, training never sees this cohort)
loco_cohorts = {
    "GSE28702 (LOCO)": {"n": 83,  "auc": 0.530, "group": "mFOLFOX6"},
    "GSE69657 (LOCO)": {"n": 30,  "auc": 0.638, "group": "FOLFOX4/XELOX"},
    "GSE72970 (LOCO)": {"n": 36,  "auc": 0.311, "group": "FOLFOX"},
}

# ============================================================
# 2. Meta-analysis functions
# ============================================================

def auc_to_fisher_z(auc):
    """Convert AUC to Fisher's z."""
    return 0.5 * log((1 + auc) / (1 - auc))

def fisher_z_to_auc(z):
    """Convert Fisher's z back to AUC."""
    return (exp(2 * z) - 1) / (exp(2 * z) + 1)

def meta_analysis_random_effects(cohort_data, label="Random-effects meta-analysis"):
    """
    DerSimonian-Laird random-effects meta-analysis.
    
    Parameters
    ----------
    cohort_data : list of dict
        Each dict: {"name": str, "y": float (Fisher's z), "v": float (variance = 1/(n-3))}
    label : str
    
    Returns
    -------
    dict with pooled results
    """
    k = len(cohort_data)
    names = [c["name"] for c in cohort_data]
    y = np.array([c["y"] for c in cohort_data])
    v = np.array([c["v"] for c in cohort_data])
    w_fixed = 1.0 / v  # fixed-effect weights
    
    # Fixed-effect estimate
    y_fixed = np.sum(y * w_fixed) / np.sum(w_fixed)
    
    # Q statistic (Cochran's Q)
    Q = np.sum(w_fixed * (y - y_fixed) ** 2)
    
    # DerSimonian-Laird tau²
    denom = np.sum(w_fixed) - np.sum(w_fixed ** 2) / np.sum(w_fixed)
    tau2 = max(0, (Q - (k - 1)) / denom)
    
    # Random-effects weights
    w_random = 1.0 / (v + tau2)
    
    # Random-effects pooled estimate
    y_pooled = np.sum(y * w_random) / np.sum(w_random)
    se_pooled = sqrt(1.0 / np.sum(w_random))
    
    # 95% CI for pooled z
    z_alpha = 1.96
    z_ci_l = y_pooled - z_alpha * se_pooled
    z_ci_u = y_pooled + z_alpha * se_pooled
    
    # Convert back to AUC
    auc_pooled = fisher_z_to_auc(y_pooled)
    auc_ci_l = fisher_z_to_auc(z_ci_l)
    auc_ci_u = fisher_z_to_auc(z_ci_u)
    
    # Heterogeneity
    I2 = max(0, (Q - (k - 1)) / Q * 100) if Q > 0 else 0
    
    # Q-test p-value
    p_hetero = 1 - stats.chi2.cdf(Q, df=k - 1)
    
    # Test if pooled AUC ≠ 0.5
    z_test = y_pooled / se_pooled
    p_pooled = 2 * (1 - stats.norm.cdf(abs(z_test)))  # two-sided
    
    return {
        "k": k,
        "Q": Q,
        "tau2": tau2,
        "I2": I2,
        "p_hetero": p_hetero,
        "y_pooled": y_pooled,
        "se_pooled": se_pooled,
        "auc_pooled": auc_pooled,
        "auc_ci_l": auc_ci_l,
        "auc_ci_u": auc_ci_u,
        "z_test": z_test,
        "p_pooled": p_pooled,
        "names": names,
        "y": y,
        "v": v,
        "w_random": w_random,
    }


def ci_from_n_auc(n, auc):
    """
    Estimate 95% CI for AUC using normal approximation:
    SE(AUC) ~ sqrt(AUC * (1-AUC) / n)
    """
    se = sqrt(auc * (1 - auc) / n)
    ci_l = max(0.5, auc - 1.96 * se)
    ci_u = min(1.0, auc + 1.96 * se)
    return ci_l, ci_u


# ============================================================
# 3. Run analyses
# ============================================================

all_analyses = {}

# --------------------------------------------------
# Analysis A: All 4 validation cohorts (excl. training)
# --------------------------------------------------
validation_cohorts = ["GSE39582", "GSE104645", "GSE28702", "GSE69657"]
analysis_A_data = []
for name in validation_cohorts:
    c = cohorts[name]
    z = auc_to_fisher_z(c["auc"])
    v = 1.0 / (c["n"] - 3)
    analysis_A_data.append({"name": name, "y": z, "v": v})

result_A = meta_analysis_random_effects(analysis_A_data, "4 Validation Cohorts (REML)")
all_analyses["A_4_validation"] = result_A

print("=" * 60)
print("ANALYSIS A: 4 Validation Cohorts (excl. training)")
print("=" * 60)
for name in validation_cohorts:
    c = cohorts[name]
    print(f"  {name:15s}  N={c['n']:3d}  AUC={c['auc']:.3f}  [{c['ci_l']:.3f}-{c['ci_u']:.3f}]")
print(f"  {'─' * 50}")
print(f"  Pooled AUC = {result_A['auc_pooled']:.3f}  "
      f"[{result_A['auc_ci_l']:.3f}-{result_A['auc_ci_u']:.3f}]")
print(f"  Heterogeneity: Q={result_A['Q']:.2f}, I²={result_A['I2']:.1f}%, "
      f"p_hetero={result_A['p_hetero']:.4f}")
print(f"  Z-test (pooled AUC ≠ 0.5): z={result_A['z_test']:.2f}, "
      f"p={result_A['p_pooled']:.4f}")

# --------------------------------------------------
# Analysis B: All 6 cohorts (incl. training + GSE72970, sensitivity)
# --------------------------------------------------
analysis_B_data = []
for name in ["Training (merged)", "GSE39582", "GSE104645", "GSE28702", "GSE69657", "GSE72970"]:
    c = cohorts[name]
    z = auc_to_fisher_z(c["auc"])
    v = 1.0 / (c["n"] - 3)
    analysis_B_data.append({"name": name, "y": z, "v": v})

result_B = meta_analysis_random_effects(analysis_B_data, "All 6 Cohorts (sensitivity)")
all_analyses["B_all_6"] = result_B

print("\n")
print("=" * 60)
print("ANALYSIS B: All 6 Cohorts (incl. training + GSE72970, sensitivity)")
print("=" * 60)
for d in analysis_B_data:
    print(f"  {d['name']:20s}  z={d['y']:.4f}  v={d['v']:.6f}")
print(f"  Pooled AUC = {result_B['auc_pooled']:.3f}  "
      f"[{result_B['auc_ci_l']:.3f}-{result_B['auc_ci_u']:.3f}]")
print(f"  Heterogeneity: I²={result_B['I2']:.1f}%, p={result_B['p_hetero']:.4f}")

# --------------------------------------------------
# Analysis C: Subgroup — Pure XELOX vs Other Oxaliplatin
# --------------------------------------------------
xelox_cohorts = ["GSE39582", "GSE69657"]  # XELOX-like or FOLFOX4/XELOX
other_cohorts = ["GSE104645", "GSE28702", "GSE72970"]  # Oxaliplatin, mFOLFOX6, FOLFOX

def subgroup_analysis(name_list, group_label):
    data = []
    for name in name_list:
        c = cohorts[name]
        z = auc_to_fisher_z(c["auc"])
        v = 1.0 / (c["n"] - 3)
        data.append({"name": name, "y": z, "v": v})
    result = meta_analysis_random_effects(data, group_label)
    return result

result_C_xelox = subgroup_analysis(xelox_cohorts, "Pure XELOX / XELOX-like")
result_C_other = subgroup_analysis(other_cohorts, "Other Oxaliplatin Regimens")
all_analyses["C_xelox"] = result_C_xelox
all_analyses["C_other"] = result_C_other

print("\n")
print("=" * 60)
print("ANALYSIS C: Subgroup Analysis")
print("=" * 60)
print(f"  XELOX/XELOX-like: pooled AUC={result_C_xelox['auc_pooled']:.3f} "
      f"[{result_C_xelox['auc_ci_l']:.3f}-{result_C_xelox['auc_ci_u']:.3f}], "
      f"I²={result_C_xelox['I2']:.1f}%")
print(f"  Other Oxaliplatin: pooled AUC={result_C_other['auc_pooled']:.3f} "
      f"[{result_C_other['auc_ci_l']:.3f}-{result_C_other['auc_ci_u']:.3f}], "
      f"I²={result_C_other['I2']:.1f}%")

# Test for subgroup difference
z_diff = (result_C_xelox["y_pooled"] - result_C_other["y_pooled"]) / \
         sqrt(result_C_xelox["se_pooled"]**2 + result_C_other["se_pooled"]**2)
p_diff = 2 * (1 - stats.norm.cdf(abs(z_diff)))
print(f"  Subgroup difference: z={z_diff:.2f}, p={p_diff:.4f}")


# --------------------------------------------------
# Analysis D: LOCO-based meta-analysis (3 held-out cohorts)
# --------------------------------------------------
analysis_D_data = []
loco_auc_ci_estimates = {}
for name, lc in loco_cohorts.items():
    z = auc_to_fisher_z(lc["auc"])
    v = 1.0 / (lc["n"] - 3)
    analysis_D_data.append({"name": name, "y": z, "v": v})
    # Estimate CI for display from Fisher-z
    se_z = sqrt(v)
    ci_l_z = z - 1.96 * se_z
    ci_u_z = z + 1.96 * se_z
    loco_auc_ci_estimates[name] = {
        "auc": lc["auc"], "ci_l": fisher_z_to_auc(ci_l_z), "ci_u": fisher_z_to_auc(ci_u_z)
    }

result_D = meta_analysis_random_effects(analysis_D_data, "3 LOCO Held-Out Cohorts")
all_analyses["D_loco_3"] = result_D

print("\n")
print("=" * 60)
print("ANALYSIS D: LOCO Independent Validation (3 held-out cohorts)")
print("=" * 60)
for name, ci in loco_auc_ci_estimates.items():
    print(f"  {name:20s}  N={loco_cohorts[name]['n']:2d}  "
          f"AUC={ci['auc']:.3f} [{ci['ci_l']:.3f}-{ci['ci_u']:.3f}]")
print(f"  {'─' * 50}")
print(f"  Pooled AUC = {result_D['auc_pooled']:.3f}  "
      f"[{result_D['auc_ci_l']:.3f}-{result_D['auc_ci_u']:.3f}]")
print(f"  Heterogeneity: Q={result_D['Q']:.2f}, I²={result_D['I2']:.1f}%, "
      f"p_hetero={result_D['p_hetero']:.4f}")

# --------------------------------------------------
# Analysis E: Leave-One-Out sensitivity analysis
# --------------------------------------------------
print("\n")
print("=" * 60)
print("ANALYSIS E: Leave-One-Out Sensitivity Analysis (4-val)")
print("=" * 60)
loo_results = {}
loo_cohort_list = ["GSE39582", "GSE104645", "GSE28702", "GSE69657", "GSE72970"]
for left_out in loo_cohort_list:
    loo_data = []
    for name in loo_cohort_list:
        if name == left_out:
            continue
        c = cohorts[name]
        z = auc_to_fisher_z(c["auc"])
        v = 1.0 / (c["n"] - 3)
        loo_data.append({"name": name, "y": z, "v": v})
    res = meta_analysis_random_effects(loo_data, f"LOO (excl. {left_out})")
    loo_results[left_out] = {
        "pooled_auc": res["auc_pooled"],
        "ci_l": res["auc_ci_l"],
        "ci_u": res["auc_ci_u"],
        "I2": res["I2"],
    }
    print(f"  Excluding {left_out:15s}: pooled AUC={res['auc_pooled']:.3f} "
          f"[{res['auc_ci_l']:.3f}-{res['auc_ci_u']:.3f}], I²={res['I2']:.1f}%")
all_analyses["E_loo"] = loo_results

# ============================================================
# 4. Forest plot
# ============================================================
def plot_forest(result, title, filename, show_training=False):
    """
    Draw forest plot for meta-analysis.
    """
    k = result["k"]
    names = result["names"]
    y = result["y"]
    v = result["v"]
    w_random = result["w_random"]
    
    # Individual AUC and CI for display
    aucs = [fisher_z_to_auc(yi) for yi in y]
    auc_ci_l = [fisher_z_to_auc(yi - 1.96 * sqrt(vi)) for yi, vi in zip(y, v)]
    auc_ci_u = [fisher_z_to_auc(yi + 1.96 * sqrt(vi)) for yi, vi in zip(y, v)]
    
    # Pooled values for diamond
    pooled_auc = result["auc_pooled"]
    pooled_ci_l = result["auc_ci_l"]
    pooled_ci_u = result["auc_ci_u"]
    I2 = result["I2"]
    p_hetero = result["p_hetero"]
    p_pooled = result["p_pooled"]
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 3.5 + 0.5 * k))
    
    # y-positions
    y_pos = list(range(k, 0, -1))  # top to bottom
    pooled_y = 0
    y_offset = 1
    
    # Plot each study
    for i in range(k):
        yi = y_pos[i] * y_offset
        xi = aucs[i]
        xi_l = auc_ci_l[i]
        xi_u = auc_ci_u[i]
        
        # Point estimate
        marker_size = 8 + 20 * w_random[i] / max(w_random) if max(w_random) > 0 else 8
        ax.plot(xi, yi, 's', color='#2c3e50', markersize=np.sqrt(marker_size),
                zorder=5)
        
        # CI line
        ax.plot([xi_l, xi_u], [yi, yi], '-', color='#2c3e50', linewidth=1.5, zorder=4)
        
        # Study label
        weight_str = f"{w_random[i] / sum(w_random) * 100:.1f}%"
        label_text = f"{names[i]}  {aucs[i]:.3f} [{xi_l:.3f}, {xi_u:.3f}]  {weight_str}"
        ax.text(0.42, yi, label_text, va='center', ha='right', fontsize=9,
                fontfamily='monospace')
    
    # Pooled diamond
    pooled_y = (k + 1) * y_offset
    diamond_center = pooled_auc
    diamond_left = pooled_ci_l
    diamond_right = pooled_ci_u
    diamond_height = 0.3
    
    # Diamond shape
    ax.fill_between([diamond_left, diamond_center],
                     [pooled_y, pooled_y + diamond_height],
                     [pooled_y, pooled_y - diamond_height],
                     color='#c0392b', alpha=0.8, zorder=6)
    ax.fill_between([diamond_center, diamond_right],
                     [pooled_y + diamond_height, pooled_y],
                     [pooled_y - diamond_height, pooled_y],
                     color='#c0392b', alpha=0.8, zorder=6)
    
    # Pooled label
    pooled_text = (f"Random-effects (DL) pooled AUC  "
                   f"{pooled_auc:.3f} [{pooled_ci_l:.3f}, {pooled_ci_u:.3f}]")
    ax.text(0.42, pooled_y, pooled_text, va='center', ha='right', fontsize=10,
            fontweight='bold', fontfamily='monospace', color='#c0392b')
    
    # Reference line at AUC = 0.5
    ax.axvline(x=0.5, color='gray', linestyle='--', linewidth=0.8, alpha=0.6, zorder=2)
    ax.axvline(x=pooled_auc, color='#c0392b', linestyle=':', linewidth=0.8, alpha=0.4, zorder=2)
    
    # Heterogeneity annotation
    het_text = (f"Heterogeneity: Q={result['Q']:.2f}, I²={I2:.1f}%, "
                f"p_hetero={p_hetero:.4f}\n"
                f"Test pooled AUC ≠ 0.5: z={result['z_test']:.2f}, p={p_pooled:.4f}")
    ax.text(0.98, 0.02, het_text, transform=ax.transAxes, fontsize=8,
            ha='right', va='bottom', fontfamily='monospace', color='#555555',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#f9f9f9', alpha=0.8))
    
    # Axis labels
    ax.set_xlabel('AUC (95% CI)', fontsize=11)
    ax.set_title(title, fontsize=12, fontweight='bold', pad=10)
    
    # Y-axis: hide ticks
    ax.set_ylim(0, (k + 2) * y_offset + 0.5)
    ax.set_yticks([])
    
    # X-axis
    ax.set_xlim(0.3, 1.0)
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
    ax.grid(axis='x', alpha=0.3, linestyle=':')
    
    # Border styling
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    for spine in ['bottom', 'left']:
        ax.spines[spine].set_color('#cccccc')
    
    plt.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, filename), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  → Saved {filename}")
    
    return fig

# Generate forest plots
plot_forest(result_A, "PRS Multi-Cohort Meta-Analysis\n4 Validation Cohorts (Random Effects)",
            "prs_meta_forest_4val.png")
plot_forest(result_B, "PRS Multi-Cohort Meta-Analysis (Sensitivity)\nAll 6 Cohorts (Random Effects)",
            "prs_meta_forest_all5.png")

# LOCO forest plot (D)
plot_forest(result_D, "PRS LOCO Independent Validation\n3 Held-Out Cohorts (Random Effects)",
            "prs_meta_forest_loco.png")

# --------------------------------------------------
# Leave-One-Out sensitivity forest plot (E)
# --------------------------------------------------
def plot_loo_forest(loo_results, full_result, title, filename):
    """Draw leave-one-out sensitivity forest plot."""
    k_loo = len(loo_results)
    fig, ax = plt.subplots(1, 1, figsize=(10, 3.0 + 0.5 * (k_loo + 1)))
    
    y_offset = 1
    y_pos = list(range(k_loo, 0, -1))
    
    # LOO results
    for i, (left_out, res) in enumerate(loo_results.items()):
        yi = y_pos[i] * y_offset
        ax.plot(res["pooled_auc"], yi, 'o', color='#2980b9', markersize=8, zorder=5)
        ax.plot([res["ci_l"], res["ci_u"]], [yi, yi], '-', color='#2980b9', linewidth=1.5, zorder=4)
        label_text = f"Excluding {left_out:15s}  {res['pooled_auc']:.3f} [{res['ci_l']:.3f}, {res['ci_u']:.3f}]  I²={res['I2']:.1f}%"
        ax.text(0.42, yi, label_text, va='center', ha='right', fontsize=9, fontfamily='monospace')
    
    # Full analysis reference
    full_y = (k_loo + 1) * y_offset
    ax.plot(full_result["auc_pooled"], full_y, 'D', color='#c0392b', markersize=9, zorder=6)
    ax.plot([full_result["auc_ci_l"], full_result["auc_ci_u"]], [full_y, full_y], '-', color='#c0392b', linewidth=2, zorder=5)
    label_full = f"All 5 cohorts  {full_result['auc_pooled']:.3f} [{full_result['auc_ci_l']:.3f}, {full_result['auc_ci_u']:.3f}]  I²={full_result['I2']:.1f}%"
    ax.text(0.42, full_y, label_full, va='center', ha='right', fontsize=10, fontweight='bold', fontfamily='monospace', color='#c0392b')
    
    # Reference line at AUC = 0.5
    ax.axvline(x=0.5, color='gray', linestyle='--', linewidth=0.8, alpha=0.6, zorder=2)
    
    ax.set_xlabel('Pooled AUC (95% CI)', fontsize=11)
    ax.set_title(title, fontsize=12, fontweight='bold', pad=10)
    ax.set_ylim(0, (k_loo + 2) * y_offset + 0.5)
    ax.set_yticks([])
    ax.set_xlim(0.3, 1.0)
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
    ax.grid(axis='x', alpha=0.3, linestyle=':')
    
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    for spine in ['bottom', 'left']:
        ax.spines[spine].set_color('#cccccc')
    
    plt.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, filename), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  → Saved {filename}")

plot_loo_forest(loo_results, result_A,
                "PRS Leave-One-Out Sensitivity Analysis\n(5 Validation Cohorts)",
                "prs_meta_forest_loo.png")


# ============================================================
# 5. Funnel plot (publication bias)
# ============================================================
def plot_funnel(result, title, filename):
    """Draw funnel plot with pseudo-95% CI lines."""
    y = result["y"]
    v = result["v"]
    y_pooled = result["y_pooled"]
    se_pooled = result["se_pooled"]
    k = result["k"]
    
    se = np.sqrt(v)
    
    fig, ax = plt.subplots(1, 1, figsize=(7, 6))
    
    # Pseudo 95% CI funnel
    se_grid = np.linspace(0, max(se) * 1.15, 100)
    upper = y_pooled + 1.96 * se_grid
    lower = y_pooled - 1.96 * se_grid
    
    ax.fill_betweenx(se_grid, lower, upper, alpha=0.1, color='#3498db', 
                     label='Pseudo 95% CI')
    ax.plot(upper, se_grid, '--', color='#3498db', linewidth=0.8, alpha=0.5)
    ax.plot(lower, se_grid, '--', color='#3498db', linewidth=0.8, alpha=0.5)
    
    # Pooled estimate line
    ax.axvline(x=y_pooled, color='#c0392b', linestyle='-', linewidth=1, alpha=0.7)
    
    # Individual studies
    colors = plt.cm.Set2(np.linspace(0, 1, k))
    for i in range(k):
        ax.scatter(y[i], se[i], s=50, color=colors[i], edgecolors='#333333',
                   linewidth=0.5, zorder=5, label=result["names"][i])
    
    # Invert y-axis (SE increases downward)
    ax.invert_yaxis()
    ax.set_ylabel('Standard Error (SE)', fontsize=11)
    ax.set_xlabel("Fisher's z", fontsize=11)
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.legend(fontsize=8, loc='lower right')
    
    # Border styling
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    for spine in ['bottom', 'left']:
        ax.spines[spine].set_color('#cccccc')
    
    plt.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, filename), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  → Saved {filename}")


plot_funnel(result_A, "Funnel Plot — 4 Validation Cohorts", "prs_meta_funnel_4val.png")
plot_funnel(result_D, "Funnel Plot — 3 LOCO Held-Out Cohorts", "prs_meta_funnel_loco.png")


# ============================================================
# 6. Summary table
# ============================================================
print("\n")
print("=" * 70)
print("SUMMARY TABLE: PRS Multi-Cohort Meta-Analysis")
print("=" * 70)
print(f"{'Analysis':30s} {'k':>3s} {'Pooled AUC':>10s} {'95% CI':>16s} {'I²':>6s} {'p_hetero':>9s} {'p_auc=0.5':>10s}")
print(f"{'─' * 84}")
for key, res in all_analyses.items():
    if key == "E_loo":
        continue  # handled separately below
    if key == "B_all_6":
        label = "All 6 Cohorts (sensitivity)"
    elif key == "D_loco_3":
        label = "D: LOCO (3 held-out)"
    else:
        label = {
            "A_4_validation": "4 Validation Cohorts",
            "C_xelox": "Subgroup: XELOX/XELOX-like",
            "C_other": "Subgroup: Other Oxaliplatin",
        }.get(key, key)
    ci_str = f"[{res['auc_ci_l']:.3f}, {res['auc_ci_u']:.3f}]"
    print(f"{label:30s} {res['k']:3d} {res['auc_pooled']:8.3f}  {ci_str:14s} "
          f"{res['I2']:5.1f}% {res['p_hetero']:8.4f} {res['p_pooled']:9.4f}")

# LOO entries (custom formatting, D handled in main loop)
for left_out, loo_res in loo_results.items():
    print(f"{'E: LOO (excl. ' + left_out + ')':30s} {'4':>3s} {loo_res['pooled_auc']:8.3f}  "
          f"[{loo_res['ci_l']:.3f}, {loo_res['ci_u']:.3f}]  "
          f"{loo_res['I2']:5.1f}%")

# ============================================================
# 7. Save results as JSON
# ============================================================
json_results = {}
for key, res in all_analyses.items():
    if key == "E_loo":
        json_results[key] = {}
        for left_out, loo_res in res.items():
            json_results[key][left_out] = {
                "pooled_auc": round(loo_res["pooled_auc"], 4),
                "ci_lower": round(loo_res["ci_l"], 4),
                "ci_upper": round(loo_res["ci_u"], 4),
                "I2_pct": round(loo_res["I2"], 2),
            }
    else:
        json_results[key] = {
            "k": res["k"],
            "Q": round(res["Q"], 4),
            "tau2": round(res["tau2"], 6),
            "I2_pct": round(res["I2"], 2),
            "p_heterogeneity": round(res["p_hetero"], 6),
            "pooled_auc": round(res["auc_pooled"], 4),
            "pooled_ci_lower": round(res["auc_ci_l"], 4),
            "pooled_ci_upper": round(res["auc_ci_u"], 4),
            "z_test": round(res["z_test"], 4),
            "p_pooled": round(res["p_pooled"], 6),
        }

with open(os.path.join(RESULTS_DIR, "prs_meta_results.json"), "w") as f:
    json.dump(json_results, f, indent=2)
print(f"\n  → Saved prs_meta_results.json")

print("\n✓ P0-2 Multi-cohort PRS meta-analysis complete!")
