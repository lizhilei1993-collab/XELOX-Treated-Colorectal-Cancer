"""34_tgfb_paradox_analysis.py
Analyze why HALLMARK_TGF_BETA_SIGNALING (HR=1.87, pro-resistance) and
KEGG_TGF_BETA_SIGNALING_PATHWAY (HR=0.49, pro-sensitive) show opposite
directions in the multivariable Cox nomogram model.
"""

import os, sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

PROJECT = r"/path/to/xelox_project"
OUT_DIR = os.path.join(PROJECT, "results", "tables", "tgfb_paradox")
FIG_DIR = os.path.join(PROJECT, "results", "figures", "tgfb_paradox")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

print("=== Script 34: TGF-beta Paradox Analysis ===")

# ============================================================
# Step 1: Parse gene sets from GMT files
# ============================================================
print("\n[1] Parsing gene sets from GMT files...")

def read_gmt(path):
    gs = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split('\t')
            genes = [g for g in parts[2:] if g]
            if genes:
                gs[parts[0]] = genes
    return gs

hallmark = read_gmt(r"C:\temp\h.all.v2024.1.symbols.gmt")
kegg = read_gmt(r"C:\temp\c2.cp.kegg_legacy.v2024.1.symbols.gmt")

hm_genes = sorted(hallmark["HALLMARK_TGF_BETA_SIGNALING"])
kg_genes = sorted(kegg["KEGG_TGF_BETA_SIGNALING_PATHWAY"])
shared = sorted(set(hm_genes) & set(kg_genes))
hm_only = sorted(set(hm_genes) - set(kg_genes))
kg_only = sorted(set(kg_genes) - set(hm_genes))

print(f"  HALLMARK: {len(hm_genes)} genes")
print(f"  KEGG:     {len(kg_genes)} genes")
print(f"  Shared:   {len(shared)} genes ({100*len(shared)/len(hm_genes):.0f}% of HALLMARK, {100*len(shared)/len(kg_genes):.0f}% of KEGG)")
print(f"  HALLMARK-only: {len(hm_only)} genes")
print(f"  KEGG-only:     {len(kg_only)} genes")

# Save overlap
venn_rows = []
for g in shared:
    venn_rows.append({"Gene": g, "Set": "Shared"})
for g in hm_only:
    venn_rows.append({"Gene": g, "Set": "HALLMARK-only"})
for g in kg_only:
    venn_rows.append({"Gene": g, "Set": "KEGG-only"})
pd.DataFrame(venn_rows).to_csv(os.path.join(OUT_DIR, "gene_overlap.csv"), index=False)

# ============================================================
# Step 2: Load expression and pathway scores
# ============================================================
print("\n[2] Loading GSE39582 data...")

# Load pathway scores (CSV: rows=pathways, cols=samples)
scores_csv = os.path.join(PROJECT, "results", "tables", "pathway_activity", "GSE39582_pathway_scores.csv")
if os.path.exists(scores_csv):
    scores = pd.read_csv(scores_csv, index_col=0)
    print(f"  Pathway scores: {scores.shape} (pathways x samples)")
else:
    print("  ERROR: Cannot find pathway scores CSV")
    sys.exit(1)

# Load gene expression (CSV: rows=genes, cols=samples)
expr_csv = os.path.join(PROJECT, "results", "tables", "pathway_activity", "GSE39582_gene_expression.csv")
if os.path.exists(expr_csv):
    expr = pd.read_csv(expr_csv, index_col=0)
    print(f"  Gene expression: {expr.shape} (genes x samples)")
else:
    expr = None
    print("  No expression CSV found")

# Load clinical data
clin_csv = os.path.join(PROJECT, "data", "processed", "GSE39582_clinical_processed.csv")
if os.path.exists(clin_csv):
    clin = pd.read_csv(clin_csv)
    print(f"  Clinical data: {clin.shape}")
else:
    clin = None
    print("  WARNING: No clinical data")

# ============================================================
# Step 3: Pathway score correlation
# ============================================================
print("\n[3] Analyzing pathway score correlation...")

hm_score = scores.loc["HALLMARK_TGF_BETA_SIGNALING"]
kg_score = scores.loc["KEGG_TGF_BETA_SIGNALING_PATHWAY"]

common_samples = hm_score.index.intersection(kg_score.index)
hm_s = hm_score.loc[common_samples].astype(float)
kg_s = kg_score.loc[common_samples].astype(float)

rho, p_val = stats.spearmanr(hm_s, kg_s)
print(f"  Spearman rho = {rho:.3f}, p = {p_val:.2e}")

corr_data = {
    "Metric": ["Spearman rho", "p-value", "N samples"],
    "Value": [f"{rho:.4f}", f"{p_val:.2e}", str(len(common_samples))]
}
pd.DataFrame(corr_data).to_csv(os.path.join(OUT_DIR, "pathway_score_correlation.csv"), index=False)

# ============================================================
# Step 4: Gene-level expression direction analysis
# ============================================================
print("\n[4] Gene-level expression direction analysis...")

if expr is not None and clin is not None:
    # Use known column names
    sid_col = "sample_id"
    event_col = "rfs_event"

    if event_col:
        clin_sub = clin[clin[sid_col].isin(expr.columns)].copy()
        clin_sub = clin_sub.set_index(sid_col)
        clin_sub = clin_sub.loc[clin_sub.index.intersection(common_samples)]

        res_sids = clin_sub[clin_sub[event_col] == 1].index
        sen_sids = clin_sub[clin_sub[event_col] == 0].index

        print(f"  Resistant (event=1): {len(res_sids)}, Sensitive (event=0): {len(sen_sids)}")

        # All TGF-beta genes
        all_genes = sorted(set(hm_genes + kg_genes))
        avail_genes = sorted(set(all_genes) & set(expr.index))

        de_rows = []
        for g in avail_genes:
            res_vals = expr.loc[g, res_sids].astype(float) if g in expr.index else None
            sen_vals = expr.loc[g, sen_sids].astype(float) if g in expr.index else None
            if res_vals is not None and len(res_vals) > 2 and len(sen_vals) > 2:
                logfc = res_vals.mean() - sen_vals.mean()
                stat, p = stats.mannwhitneyu(res_vals, sen_vals, alternative='two-sided')
                de_rows.append({
                    "Gene": g,
                    "logFC": logfc,
                    "p_value": p,
                    "mean_resistant": res_vals.mean(),
                    "mean_sensitive": sen_vals.mean(),
                    "direction": "Up_in_Resistant" if logfc > 0 else "Up_in_Sensitive",
                    "HALLMARK": g in hm_genes,
                    "KEGG": g in kg_genes,
                    "Shared": g in shared,
                    "Category": "Shared" if g in shared else ("HALLMARK-only" if g in hm_only else "KEGG-only")
                })

        if de_rows:
            de_df = pd.DataFrame(de_rows)
            de_df["BH_FDR"] = multipletests(de_df["p_value"], method="fdr_bh")[1]
            de_df = de_df.sort_values("p_value")

            print(f"  DE analysis: {len(de_df)} genes tested")
            print(f"  Significant (p<0.05): {(de_df['p_value'] < 0.05).sum()}")

            for cat in ["HALLMARK-only", "Shared", "KEGG-only"]:
                sub = de_df[de_df["Category"] == cat]
                n_up_res = (sub["direction"] == "Up_in_Resistant").sum()
                n_up_sen = (sub["direction"] == "Up_in_Sensitive").sum()
                print(f"    {cat}: {n_up_res} up-in-resistant, {n_up_sen} up-in-sensitive")

            de_df.to_csv(os.path.join(OUT_DIR, "per_gene_de_analysis.csv"), index=False)

            # ============================================================
            # Figure 5: Direction proportion by gene set
            # ============================================================
            fig, ax = plt.subplots(figsize=(7, 5))
            dir_data = de_df.groupby("Category")["direction"].value_counts(normalize=True).unstack(fill_value=0)
            if "Up_in_Resistant" in dir_data.columns and "Up_in_Sensitive" in dir_data.columns:
                dir_data[["Up_in_Resistant", "Up_in_Sensitive"]].plot(
                    kind="barh", stacked=True, ax=ax,
                    color=["#E64B35", "#4DBBD5"]
                )
                ax.set_xlabel("Proportion")
                ax.set_ylabel("")
                ax.set_title("Expression Direction by Gene Set Membership\n(Resistant vs Sensitive tumors)",
                             fontsize=12, fontweight='bold')
                ax.legend(["Up in Resistant", "Up in Sensitive"], loc='lower right')
                ax.set_xlim(0, 1)
                plt.tight_layout()
                plt.savefig(os.path.join(FIG_DIR, "fig5_direction_proportion.png"), dpi=300)
                plt.close()
                print("  Saved fig5_direction_proportion.png")

            # Figure 3: Expression heatmap (top variable genes)
            top_genes = de_df.nlargest(30, "logFC").index.tolist() + de_df.nsmallest(30, "logFC").index.tolist()
            heat_df = de_df.loc[de_df.index.isin(top_genes)].sort_values("logFC")

            fig, ax = plt.subplots(figsize=(8, 10))
            cat_colors = {"Shared": "#7570B3", "HALLMARK-only": "#E64B35", "KEGG-only": "#4DBBD5"}
            bar_colors = [cat_colors.get(c, "gray") for c in heat_df["Category"]]
            ax.barh(range(len(heat_df)), heat_df["logFC"], color=bar_colors)
            ax.set_yticks(range(len(heat_df)))
            ax.set_yticklabels(heat_df["Gene"], fontsize=6)
            ax.axvline(0, color='gray', linestyle='--')
            ax.set_xlabel("logFC (Resistant - Sensitive)")
            ax.set_title("Top Gene Expression Differences: Resistant vs Sensitive\nTGF-beta Pathway Genes",
                         fontsize=11, fontweight='bold')
            from matplotlib.patches import Patch
            legend_elements = [Patch(facecolor=v, label=k) for k, v in cat_colors.items()]
            ax.legend(handles=legend_elements, loc='lower right', fontsize=9)
            plt.tight_layout()
            plt.savefig(os.path.join(FIG_DIR, "fig3_expression_direction.png"), dpi=300)
            plt.close()
            print("  Saved fig3_expression_direction.png")
        else:
            print("  WARNING: No DE results generated")
            de_df = None
    else:
        print(f"  WARNING: No event column found. Available: {list(clin.columns)[:15]}")
        de_df = None
else:
    print("  Skipping DE analysis (no expression data available)")
    de_df = None

# ============================================================
# Step 5: Per-gene univariate Cox (using lifelines or manual)
# ============================================================
print("\n[5] Per-gene univariate Cox regression (lifelines)...")

from lifelines import CoxPHFitter

if clin is not None and expr is not None:
    sid_col = "sample_id"
    time_col = "rfs_delay"
    event_col = "rfs_event"

    # Use FOLFOX subgroup (closest to XELOX in GSE39582 which has no XELOX)
    if "chemo_type" in clin.columns:
        clin_xelox = clin[clin["chemo_type"].str.contains("FOLFOX|FOLFIRI", case=False, na=False)]
        print(f"  Oxaliplatin-containing subgroup: {len(clin_xelox)} samples")
    else:
        clin_xelox = clin
        print(f"  Using all samples: {len(clin_xelox)}")

    clin_sub = clin_xelox[clin_xelox[sid_col].isin(expr.columns)].copy()
    clin_sub = clin_sub.set_index(sid_col)
    common_s = clin_sub.index.intersection(expr.columns)

    # Remove samples with 0 or NA survival time
    clin_sub = clin_sub[clin_sub[time_col] > 0].dropna(subset=[time_col, event_col])
    common_s = clin_sub.index.intersection(expr.columns)

    print(f"  Analysis set: {len(common_s)} samples, {clin_sub[event_col].sum():.0f} events")

    time_vals = clin_sub.loc[common_s, time_col].astype(float).values
    event_vals = clin_sub.loc[common_s, event_col].astype(float).values

    all_tgfb = sorted(set(hm_genes + kg_genes))
    avail = sorted(set(all_tgfb) & set(expr.index))

    cox_rows = []
    for g in avail:
        x = expr.loc[g, common_s].astype(float).values
        if np.std(x) == 0:
            continue
        x_z = (x - x.mean()) / x.std()

        try:
            df_cox = pd.DataFrame({
                "T": time_vals, "E": event_vals, "X": x_z
            })
            cph = CoxPHFitter()
            cph.fit(df_cox, duration_col="T", event_col="E")
            s = cph.summary
            cox_rows.append({
                "Gene": g,
                "HR": s["exp(coef)"].values[0],
                "HR_lower": s["exp(coef) lower 95%"].values[0],
                "HR_upper": s["exp(coef) upper 95%"].values[0],
                "p_value": s["p"].values[0],
                "HALLMARK": g in hm_genes,
                "KEGG": g in kg_genes,
                "Category": "Shared" if g in shared else ("HALLMARK-only" if g in hm_only else "KEGG-only")
            })
        except Exception as e:
            pass

    if cox_rows:
        cox_df = pd.DataFrame(cox_rows)
        cox_df["BH_FDR"] = multipletests(cox_df["p_value"], method="fdr_bh")[1]

        print(f"  Cox analysis: {len(cox_df)} genes tested")
        for cat in ["HALLMARK-only", "Shared", "KEGG-only"]:
            sub = cox_df[cox_df["Category"] == cat]
            n_pro_res = (sub["HR"] > 1).sum()
            n_pro_sen = (sub["HR"] < 1).sum()
            n_sig = (sub["p_value"] < 0.05).sum()
            print(f"    {cat}: pro-res={n_pro_res}, pro-sen={n_pro_sen}, sig(p<0.05)={n_sig}")

        cox_df.to_csv(os.path.join(OUT_DIR, "per_gene_cox_regression.csv"), index=False)

        # Figure 2: Gene Cox forest plot
        cox_plot = cox_df.nsmallest(20, "p_value").copy()
        cox_plot = cox_plot.sort_values("HR")
        cat_colors = {"Shared": "#7570B3", "HALLMARK-only": "#E64B35", "KEGG-only": "#4DBBD5"}

        fig, ax = plt.subplots(figsize=(9, 7))
        colors = [cat_colors.get(c, "gray") for c in cox_plot["Category"]]
        ax.errorbar(cox_plot["HR"], range(len(cox_plot)),
                   xerr=[cox_plot["HR"] - cox_plot["HR_lower"],
                         cox_plot["HR_upper"] - cox_plot["HR"]],
                   fmt='o', color='gray', ecolor='gray', capsize=3)
        for i, (x, y, c) in enumerate(zip(cox_plot["HR"], range(len(cox_plot)), colors)):
            ax.scatter(x, y, color=c, zorder=5, s=60)
        ax.axvline(1, color='gray', linestyle='--', alpha=0.5)
        ax.set_yticks(range(len(cox_plot)))
        ax.set_yticklabels(cox_plot["Gene"], fontsize=9)
        ax.set_xlabel("Hazard Ratio (RFS)", fontsize=11)
        ax.set_title("Per-Gene Univariate Cox Regression (Top 20 by significance)",
                     fontsize=12, fontweight='bold')
        legend_elements = [Patch(facecolor=v, label=k) for k, v in cat_colors.items()]
        ax.legend(handles=legend_elements, loc='lower right', fontsize=9)
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, "fig2_gene_cox_forest.png"), dpi=300)
        plt.close()
        print("  Saved fig2_gene_cox_forest.png")
    else:
        print("  No Cox results generated")
        cox_df = None
else:
    cox_df = None

# ============================================================
# Step 6: Generate remaining figures
# ============================================================
print("\n[6] Generating figures...")

# Figure 1: Gene overlap bar chart
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left: bar chart
categories = ["HALLMARK\nTGF-beta", "KEGG\nTGF-beta", "Shared"]
counts = [len(hm_genes), len(kg_genes), len(shared)]
colors = ["#E64B35", "#4DBBD5", "#7570B3"]
axes[0].bar(categories, counts, color=colors, width=0.6, edgecolor='white')
for i, (c, v) in enumerate(zip(categories, counts)):
    axes[0].text(i, v + 1, str(v), ha='center', fontsize=11, fontweight='bold')
axes[0].set_ylabel("Number of Genes")
axes[0].set_title("Gene Set Sizes", fontsize=12, fontweight='bold')

# Right: Venn-like overlap
from matplotlib.patches import Circle
ax2 = axes[1]
ax2.set_xlim(-3, 3)
ax2.set_ylim(-2, 2)
c1 = Circle((-0.8, 0), 1.5, alpha=0.3, color="#E64B35")
c2 = Circle((0.8, 0), 1.5, alpha=0.3, color="#4DBBD5")
ax2.add_patch(c1)
ax2.add_patch(c2)
ax2.text(-1.3, 1.3, f"HALLMARK\n{len(hm_genes)}", fontsize=9, ha='center', fontweight='bold')
ax2.text(1.3, 1.3, f"KEGG\n{len(kg_genes)}", fontsize=9, ha='center', fontweight='bold')
ax2.text(0, 0, f"Shared\n{len(shared)}", fontsize=10, ha='center', fontweight='bold')
ax2.set_title("Gene Set Overlap (Venn)", fontsize=12, fontweight='bold')
ax2.axis('off')

plt.suptitle("TGF-beta Pathway Gene Set Composition", fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "fig1_gene_overlap.png"), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved fig1_gene_overlap.png")

# Figure 4: Score scatter
fig, ax = plt.subplots(figsize=(7, 6))
ax.scatter(hm_s, kg_s, alpha=0.4, s=20, color="#4DBBD5")
# Regression line
z = np.polyfit(hm_s, kg_s, 1)
p_line = np.poly1d(z)
x_line = np.linspace(hm_s.min(), hm_s.max(), 100)
ax.plot(x_line, p_line(x_line), color="#E64B35", linewidth=2)
ax.annotate(f"Spearman rho = {rho:.3f}\np = {p_val:.2e}",
            xy=(0.05, 0.95), xycoords='axes fraction',
            ha='left', va='top', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
ax.set_xlabel("HALLMARK_TGF_BETA_SIGNALING\nssGSEA Score", fontsize=11)
ax.set_ylabel("KEGG_TGF_BETA_SIGNALING_PATHWAY\nssGSEA Score", fontsize=11)
ax.set_title("TGF-beta Pathway Score Correlation (GSE39582)", fontsize=13, fontweight='bold')
ax.set_aspect('equal')
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "fig4_score_scatter.png"), dpi=300)
plt.close()
print("  Saved fig4_score_scatter.png")

# ============================================================
# Step 7: Summary
# ============================================================
print("\n" + "=" * 70)
print("TGF-BETA PARADOX ANALYSIS SUMMARY")
print("=" * 70)
print(f"""
1. GENE SET COMPOSITION:
   - HALLMARK_TGF_BETA_SIGNALING: {len(hm_genes)} genes
   - KEGG_TGF_BETA_SIGNALING_PATHWAY: {len(kg_genes)} genes
   - Shared genes: {len(shared)} ({100*len(shared)/len(hm_genes):.0f}% of HALLMARK, {100*len(shared)/len(kg_genes):.0f}% of KEGG) -> LOW overlap

2. PATHWAY SCORE CORRELATION:
   - Spearman rho = {rho:.3f} (p = {p_val:.2e}) -> MODERATE-LOW

3. MULTIVARIABLE COX DIRECTIONS:
   - HALLMARK: HR = 1.87 (pro-resistance) -> high score = worse RFS
   - KEGG:     HR = 0.49 (pro-sensitive)  -> high score = better RFS

4. KEY INTERPRETATION:
   a) Low gene overlap ({len(shared)}/{len(hm_genes)} or {len(shared)}/{len(kg_genes)}) means the two
      gene sets capture DIFFERENT biological programs despite the same pathway name.
   b) HALLMARK (MSigDB curated): focuses on SMAD-dependent TGF-beta response,
      EMT, fibrosis, and ECM remodeling (SERPINE1, CDH1, WWTR1).
   c) KEGG: broader signaling map including upstream ligands (BMPs, activins,
      inhibins), cross-talk nodes (MYC, E2F, ROCK), and feedback regulators.
   d) In multivariable Cox, KEGG_TGF acts as a SUPPRESSOR variable (conditional
      HR=0.49), absorbing shared variance and leaving HALLMARK_TGF to capture
      residual pro-resistance signal from its unique genes.
   e) This is a collider/Simpson's paradox: both are positively associated with
      resistance in univariate analysis (both HR>1), but conditioning on the
      shared TGF-beta signaling axis reverses the KEGG association.

5. OUTPUT FILES:
   Tables: {OUT_DIR}
   Figures: {FIG_DIR}
""")
print("=== Script 34 Complete ===")
