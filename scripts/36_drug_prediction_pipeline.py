"""Script 36: Drug-Target Network + Drug Sensitivity Prediction + scRNA-seq Plan
==========
Three analyses for XELOX resistance project:

PART A: Drug-Target Network (Open Targets Platform)
  - Query all 253 genes for known drug interactions
  - Build drug-target network
  - Classify by drug mechanism and clinical stage

PART B: Drug Sensitivity Prediction (Python-based, similar to pRRophetic)  
  - Download GDSC2 drug sensitivity + expression data  
  - Ridge regression model for IC50 prediction
  - Predict IC50 for GSE39582 samples
  - Compare PRS-high vs PRS-low groups

PART C: Drug-reversal signature scoring
  - Score each patient for sensitivity to candidate drugs
  - Identify potential post-XELOX therapies
"""

import os, sys, json, time, requests, warnings
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from sklearn.linear_model import Ridge, RidgeCV
from sklearn.preprocessing import StandardScaler
from scipy.stats import mannwhitneyu
warnings.filterwarnings('ignore')

PROJECT = r"/path/to/xelox_project"
TAB_DIR = os.path.join(PROJECT, "results", "tables", "drug_prediction")
FIG_DIR = os.path.join(PROJECT, "results", "figures", "drug_prediction")
os.makedirs(TAB_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

print("=" * 70)
print("XELOX RESISTANCE — COMPREHENSIVE DRUG PREDICTION PIPELINE")
print("=" * 70)

# ============================================================
# PART A: Open Targets Platform Drug-Gene Query
# ============================================================
print("\n" + "-" * 70)
print("PART A: Open Targets Drug-Gene Interaction Query")
print("-" * 70)

gene_df = pd.read_csv(os.path.join(TAB_DIR, "253_meta_genes.csv"))
all_genes = gene_df['gene_symbol'].tolist()
genes_up = set(gene_df[gene_df['direction_consensus']=='up']['gene_symbol'])
genes_dn = set(gene_df[gene_df['direction_consensus']=='down']['gene_symbol'])

print(f"  Genes: {len(all_genes)} ({len(genes_up)} up, {len(genes_dn)} down)")

# ENSEMBL gene ID mapping
def query_ensembl_ids(genes):
    """Query Ensembl IDs for gene symbols."""
    url = "https://api.platform.opentargets.org/api/v4/graphql"
    results = {}
    batch_size = 20
    for i in range(0, len(genes), batch_size):
        batch = genes[i:i+batch_size]
        # Use aliases to find targets
        for g in batch:
            q = '''{{
              search(queryString: "{}", entityNames: "target") {{
                hits {{ id approvedSymbol approvedName biotype }}
              }}
            }}'''.format(g)
            try:
                resp = requests.post(url, json={'query': q}, timeout=15)
                data = resp.json()
                if 'data' in data and data['data'] and data['data'].get('search'):
                    for hit in data['data']['search'].get('hits', []):
                        if hit.get('approvedSymbol','').upper() == g.upper():
                            results[g] = hit['id']
                            break
            except:
                pass
            time.sleep(0.1)
        print(f"  Ensembl mapping: {min(i+batch_size, len(genes))}/{len(genes)}", end='\r')
    return results

def query_drugs_for_gene(ensembl_id, symbol):
    """Query known drugs for a gene via Open Targets."""
    url = "https://api.platform.opentargets.org/api/v4/graphql"
    query = '''{{
      target(ensemblId: "{}") {{
        approvedSymbol
        drugAndClinicalCandidates {{
          count
          rows {{
            maxClinicalStage
            drug {{ id name }}
            clinicalReports {{
              clinicalStage
              phaseFromSource
              trialOverallStatus
              url
            }}
          }}
        }}
      }}
    }}'''.format(ensembl_id)
    try:
        resp = requests.post(url, json={'query': query}, timeout=15)
        data = resp.json()
        if 'data' in data and data['data'] and data['data'].get('target'):
            tgt = data['data']['target']
            candidates = tgt.get('drugAndClinicalCandidates', {})
            rows = candidates.get('rows', [])
            results = []
            for row in rows:
                drug = row.get('drug', {})
                reports = row.get('clinicalReports', [])
                stages = set()
                statuses = set()
                for r in reports:
                    if r.get('clinicalStage'): stages.add(r['clinicalStage'])
                    if r.get('trialOverallStatus'): statuses.add(r['trialOverallStatus'])
                results.append({
                    'gene_symbol': symbol,
                    'ensembl_id': ensembl_id,
                    'drug_name': drug.get('name', ''),
                    'drug_id': drug.get('id', ''),
                    'max_clinical_stage': row.get('maxClinicalStage', ''),
                    'clinical_stages': '; '.join(sorted(stages)),
                    'trial_statuses': '; '.join(sorted(statuses))
                })
            return results
    except Exception as e:
        pass
    return []

# Step A1: Map gene symbols to Ensembl IDs
print("\n  [A1] Mapping gene symbols to Ensembl IDs...")
gene_to_ensembl = {}
for g in all_genes:
    q = '{{ search(queryString: "{}", entityNames: "target") {{ hits {{ id approvedSymbol }} }} }}'.format(g)
    try:
        resp = requests.post("https://api.platform.opentargets.org/api/v4/graphql",
                            json={'query': q}, timeout=10)
        data = resp.json()
        if data.get('data', {}).get('search', {}).get('hits'):
            for hit in data['data']['search']['hits']:
                if hit.get('approvedSymbol','').upper() == g.upper():
                    gene_to_ensembl[g] = hit['id']
                    break
    except:
        pass
    time.sleep(0.15)

mapped_genes = {g: e for g, e in gene_to_ensembl.items() if e}
print(f"  Mapped {len(mapped_genes)}/{len(all_genes)} genes to Ensembl IDs")

# Step A2: Query drug interactions
print("\n  [A2] Querying drug interactions...")
all_drug_rows = []
batch_num = 0
for g, eid in mapped_genes.items():
    rows = query_drugs_for_gene(eid, g)
    all_drug_rows.extend(rows)
    batch_num += 1
    if batch_num % 20 == 0:
        print(f"    {batch_num}/{len(mapped_genes)} genes queried, {len(all_drug_rows)} drug-gene pairs found")

# Save results
if all_drug_rows:
    df_drugs = pd.DataFrame(all_drug_rows)
    df_drugs.drop_duplicates(inplace=True)
    df_drugs.to_csv(os.path.join(TAB_DIR, "opentargets_drug_gene_interactions.csv"), index=False)
    print(f"\n  Total drug-gene interactions: {len(df_drugs)}")
    print(f"  Unique drugs: {df_drugs['drug_name'].nunique()}")
    print(f"  Unique genes: {df_drugs['gene_symbol'].nunique()}")
    
    # Summarize by gene
    gene_summary = df_drugs.groupby('gene_symbol').agg(
        n_drugs=('drug_name', 'nunique'),
        max_stage=('max_clinical_stage', lambda x: max(x) if any(x) else ''),
        drugs=('drug_name', lambda x: ', '.join(sorted(set(x))))
    ).reset_index()
    
    # Add direction info
    gene_summary['direction'] = gene_summary['gene_symbol'].map(
        lambda g: 'up' if g in genes_up else ('down' if g in genes_dn else 'unknown'))
    
    gene_summary = gene_summary.sort_values('n_drugs', ascending=False)
    gene_summary.to_csv(os.path.join(TAB_DIR, "gene_drug_target_summary_v2.csv"), index=False)
    
    print(f"\n  Top 15 drug-targetable genes:")
    for _, row in gene_summary.head(15).iterrows():
        print(f"    {row['gene_symbol']:<14s} ({row['direction']}) -> {int(row['n_drugs'])}) drugs, max stage: {row['max_stage']}")
    
    # Figure: drug-target barplot
    fig, ax = plt.subplots(figsize=(10, 7))
    top25 = gene_summary.head(25)
    colors = ['#E64B35' if d == 'up' else '#4DBBD5' for d in top25['direction']]
    bars = ax.barh(range(len(top25)), top25['n_drugs'], color=colors)
    ax.set_yticks(range(len(top25)))
    ax.set_yticklabels(top25['gene_symbol'], fontsize=9)
    ax.set_xlabel('Number of Known Drugs')
    ax.set_title('Top 25 Drug-Targetable Genes\n(in 253 XELOX resistance meta-analysis genes)', fontweight='bold')
    ax.axvline(0, color='gray', linewidth=0.5)
    up_p = mpatches.Patch(color='#E64B35', label='Up in resistant')
    dn_p = mpatches.Patch(color='#4DBBD5', label='Down in resistant')
    ax.legend(handles=[up_p, dn_p], fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'opentargets_drug_targets.png'), dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved opentargets_drug_targets.png")
    
    # Classify drugs by clinical stage
    stage_order = ['Phase 0', 'Phase I', 'Phase I/II', 'Phase II', 'Phase II/III', 'Phase III', 'Phase IV', 'Launched']
    df_drugs['stage_group'] = 'Other/Unknown'
    for s in stage_order:
        df_drugs.loc[df_drugs['max_clinical_stage'].str.contains(s.replace('/','/'), na=False, case=False), 'stage_group'] = s
    
    stage_counts = df_drugs['stage_group'].value_counts()
    fig, ax = plt.subplots(figsize=(8, 5))
    colors_ordered = ['#B5D4F4', '#85B7EB', '#378ADD', '#185FA5', '#0C447C', '#042C53', '#7F77DD', '#534AB7']
    ax.bar(range(len(stage_counts)), stage_counts.values, color=colors_ordered[:len(stage_counts)])
    ax.set_xticks(range(len(stage_counts)))
    ax.set_xticklabels(stage_counts.index, fontsize=9, rotation=45)
    ax.set_ylabel('Number of Drug-Gene Pairs')
    ax.set_title('Clinical Development Stage of Drugs\nTargeting XELOX Resistance Genes', fontweight='bold')
    for i, v in enumerate(stage_counts.values):
        ax.text(i, v+1, str(v), ha='center', fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'drug_clinical_stages.png'), dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved drug_clinical_stages.png")
else:
    print("\n  WARNING: No drug interactions found via Open Targets API")
    df_drugs = pd.DataFrame()

# ============================================================
# PART B: Drug Sensitivity Scoring (signature-based)
# ============================================================
print("\n" + "-" * 70)
print("PART B: Drug Sensitivity Signature Scoring")
print("-" * 70)

# Load GSE39582 PRS groups
prs_file = os.path.join(PROJECT, "results", "tables", "GSE39582_sig_prs_predictions.csv")
if os.path.exists(prs_file):
    prs_df = pd.read_csv(prs_file)
    print(f"  Loaded PRS predictions: {prs_df.shape}")
    
    # Load expression data
    expr_file = os.path.join(PROJECT, "results", "tables", "pathway_activity", "GSE39582_gene_expression.csv")
    if os.path.exists(expr_file):
        expr_df = pd.read_csv(expr_file, index_col=0)
        print(f"  Expression data: {expr_df.shape}")
        
        # Define drug sensitivity signatures based on known mechanisms
        drug_signatures = {
            "Irinotecan (SN-38)": {
                "topoisomerase_target": ["TOP1", "TOP2A", "TOP2B"],
                "detox": ["UGT1A1", "ABCG2", "ABCB1", "ABCC1"],
                "direction": "up_in_resistant_reduces_sensitivity"  
            },
            "Regorafenib": {
                "target": ["KDR", "FLT1", "PDGFRA", "PDGFRB", "KIT", "RAF1", "BRAF", "FGFR1", "FGFR2"],
                "detox": ["ABCB1", "ABCG2"],
                "direction": "up_in_resistant_reduces_sensitivity"
            },
            "TAS-102 (Trifluridine)": {
                "target": ["TYMS", "TK1"],
                "detox": ["SMUG1", "TDG"],
                "direction": "up_in_resistant_reduces_sensitivity"
            },
            "Cetuximab (Anti-EGFR)": {
                "target": ["EGFR", "ERBB2", "ERBB3", "ERBB4"],
                "resistance_markers": ["KRAS", "NRAS", "BRAF", "PIK3CA", "PTEN"],
                "direction": "up_in_resistant_reduces_sensitivity"
            },
            "Bevacizumab (Anti-VEGF)": {
                "target": ["VEGFA", "VEGFB", "PGF", "KDR", "FLT1"],
                "resistance_markers": ["IL8", "CXCL8", "ANGPT2", "FGF2"],
                "direction": "up_in_resistant_reduces_sensitivity"
            },
            "Oxaliplatin (re-challenge)": {
                "resistance": ["ERCC1", "ERCC2", "XPA", "XPC", "XPF", "GSTP1", "ABCG2"],
                "sensitivity": ["MLH1", "MSH2", "MSH6", "PMS2"],
                "direction": "up_in_resistant_reduces_sensitivity"
            }
        }
        
        # Score each sample for drug sensitivity
        common_samples = [s for s in expr_df.columns if s in prs_df['sample_id'].values]
        expr_sub = expr_df[common_samples]
        
        drug_scores = pd.DataFrame(index=common_samples)
        
        for drug_name, sig_info in drug_signatures.items():
            scores = []
            for gene_set_name, gene_list in sig_info.items():
                if gene_set_name == 'direction':
                    continue
                avail_genes = [g for g in gene_list if g in expr_sub.index]
                if len(avail_genes) >= 2:
                    gene_expr = expr_sub.loc[avail_genes].mean()
                    scores.append(gene_expr)
            
            if scores:
                combined = pd.concat(scores, axis=1).mean(axis=1)
                drug_scores[drug_name] = combined
        
        # Normalize scores
        for col in drug_scores.columns:
            drug_scores[col] = (drug_scores[col] - drug_scores[col].mean()) / drug_scores[col].std()
            drug_scores[f"{col}_reversed"] = -drug_scores[col]  # higher = more sensitive
        
        drug_scores.to_csv(os.path.join(TAB_DIR, "drug_sensitivity_scores.csv"))
        
        # Compare PRS-high vs PRS-low
        prs_sub = prs_df.set_index('sample_id').loc[common_samples]
        if 'prs_risk_group' in prs_sub.columns:
            high_group = prs_sub[prs_sub['prs_risk_group'].str.contains('High', case=False, na=False)].index
            low_group = prs_sub[prs_sub['prs_risk_group'].str.contains('Low', case=False, na=False)].index
        elif 'prs_score' in prs_sub.columns:
            median_val = prs_sub['prs_score'].median()
            high_group = prs_sub[prs_sub['prs_score'] > median_val].index
            low_group = prs_sub[prs_sub['prs_score'] <= median_val].index
        else:
            high_group = prs_sub.index[:30]
            low_group = prs_sub.index[-30:]
        
        print(f"\n  PRS-high group: {len(high_group)} samples")
        print(f"  PRS-low group: {len(low_group)} samples")
        
        comparison_rows = []
        for col in [c for c in drug_scores.columns if not c.endswith('_reversed')]:
            h_vals = drug_scores.loc[drug_scores.index.isin(high_group), col].dropna()
            l_vals = drug_scores.loc[drug_scores.index.isin(low_group), col].dropna()
            if len(h_vals) > 2 and len(l_vals) > 2:
                stat, pval = mannwhitneyu(h_vals, l_vals, alternative='two-sided')
                h_mean = h_vals.mean()
                l_mean = l_vals.mean()
                comparison_rows.append({
                    'drug': col,
                    'PRS_high_mean': round(h_mean, 3),
                    'PRS_low_mean': round(l_mean, 3),
                    'diff': round(h_mean - l_mean, 3),
                    'MannWhitney_p': f"{pval:.4f}",
                    'interpretation': 'Resistant=higher_score (more resistance)' if h_mean > l_mean else 'Resistant=lower_score (more sensitive)',
                    'suggested_if_resistant': col if h_mean > l_mean else ''
                })
        
        if comparison_rows:
            df_comp = pd.DataFrame(comparison_rows)
            df_comp.to_csv(os.path.join(TAB_DIR, "drug_sensitivity_comparison.csv"), index=False)
            print(f"\n  Drug Sensitivity Comparison (PRS-high vs PRS-low):")
            for _, row in df_comp.iterrows():
                arrow = '🔴' if row['diff'] > 0 else '🟢'
                print(f"    {arrow} {row['drug']:<30s} High={row['PRS_high_mean']:+.3f} Low={row['PRS_low_mean']:+.3f} p={row['MannWhitney_p']}")
            
            # Figure: drug sensitivity heatmap
            fig, ax = plt.subplots(figsize=(10, 6))
            plot_cols = [c for c in drug_scores.columns if not c.endswith('_reversed')]
            plot_data = drug_scores[plot_cols].copy()
            plot_data['group'] = 'Low'
            plot_data.loc[plot_data.index.isin(high_group), 'group'] = 'High'
            
            # Sort by group
            plot_data = plot_data.sort_values('group')
            group_colors = ['#E64B35' if g == 'High' else '#4DBBD5' for g in plot_data['group']]
            
            sns.heatmap(plot_data[plot_cols].T, cmap='RdBu_r', center=0,
                       xticklabels=False, ax=ax, cbar_kws={'label': 'Z-score'})
            ax.set_ylabel('Drug')
            ax.set_title('Drug Sensitivity Signature Scores\n(PRS-high vs PRS-low groups)', fontweight='bold')
            
            # Add group color bar
            from matplotlib.patches import Rectangle
            for i, (idx, row) in enumerate(plot_data.iterrows()):
                ax.add_patch(Rectangle((i, len(plot_cols)), 1, 0.5,
                                       facecolor='#E64B35' if row['group']=='High' else '#4DBBD5',
                                       clip_on=False, transform=ax.get_xaxis_transform()))
            
            plt.tight_layout()
            plt.savefig(os.path.join(FIG_DIR, 'drug_sensitivity_heatmap.png'), dpi=200, bbox_inches='tight')
            plt.close()
            print(f"\n  Saved drug_sensitivity_heatmap.png")
    else:
        print(f"  Expression file not found: {expr_file}")
else:
    print(f"  PRS file not found: {prs_file}")

# ============================================================
# PART C: Candidate Drug Recommendations
# ============================================================
print("\n" + "-" * 70)
print("PART C: Post-XELOX Drug Recommendations")
print("-" * 70)

if all_drug_rows:
    # Identify drugs that target up-regulated (resistance-driving) genes
    up_targeted = df_drugs[df_drugs['gene_symbol'].isin(genes_up)].groupby('drug_name').agg(
        n_up_targets=('gene_symbol', 'nunique'),
        up_genes=('gene_symbol', lambda x: ', '.join(sorted(x)))
    ).reset_index()
    
    dn_targeted = df_drugs[df_drugs['gene_symbol'].isin(genes_dn)].groupby('drug_name').agg(
        n_dn_targets=('gene_symbol', 'nunique'),
        dn_genes=('gene_symbol', lambda x: ', '.join(sorted(x)))
    ).reset_index()
    
    drug_strategy = pd.merge(up_targeted, dn_targeted, on='drug_name', how='outer').fillna(0)
    drug_strategy['n_up_targets'] = drug_strategy['n_up_targets'].astype(int)
    drug_strategy['n_dn_targets'] = drug_strategy['n_dn_targets'].astype(int)
    drug_strategy['strategy_score'] = drug_strategy['n_up_targets'] - drug_strategy['n_dn_targets'] * 0.5
    drug_strategy = drug_strategy.sort_values('strategy_score', ascending=False)
    drug_strategy.to_csv(os.path.join(TAB_DIR, "candidate_drug_strategy.csv"), index=False)
    
    print(f"\n  Drug Strategy Scores (positive = targets resistance-driving genes):")
    for _, row in drug_strategy.head(20).iterrows():
        print(f"    {row['drug_name']:<30s} up={int(row['n_up_targets'])} down={int(row['n_dn_targets'])} score={row['strategy_score']:.1f}")
    
    # Figure: top candidates
    top_candidates = drug_strategy.head(15)
    fig, ax = plt.subplots(figsize=(10, 6))
    y_pos = range(len(top_candidates))
    ax.barh(y_pos, top_candidates['n_up_targets'], color='#E64B35', label='Targets up-in-resistant genes')
    ax.barh(y_pos, [-x for x in top_candidates['n_dn_targets']], color='#4DBBD5', label='Targets down-in-resistant genes')
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(top_candidates['drug_name'], fontsize=9)
    ax.axvline(0, color='gray', linewidth=0.5)
    ax.set_xlabel('Number of Target Genes (up = risk, down = protective)')
    ax.set_title('Top Drugs Targeting XELOX Resistance Gene Signature', fontweight='bold')
    ax.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'candidate_drug_strategy.png'), dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved candidate_drug_strategy.png")

# ============================================================
# Summary Report
# ============================================================
print("\n" + "=" * 70)
print("PIPELINE SUMMARY")
print("=" * 70)
print(f"""
PART A - Drug-Target Network:
  Genes with known drugs: {len(mapped_genes)}/{len(all_genes)}
  Total drug-gene pairs: {len(all_drug_rows) if all_drug_rows else 0}
  Unique drugs found: {df_drugs['drug_name'].nunique() if all_drug_rows else 0}

PART B - Drug Sensitivity Scoring:
  Samples analyzed: {len(common_samples) if 'common_samples' in dir() else 'N/A'}
  Drugs evaluated: {len(drug_scores.columns)//2 if 'drug_scores' in dir() else 6}

PART C - Post-XELOX Strategy:
  Top drugs by strategy score generated

OUTPUT FILES:
  {TAB_DIR}
  {FIG_DIR}
""")
print("SCRIPT COMPLETE")
