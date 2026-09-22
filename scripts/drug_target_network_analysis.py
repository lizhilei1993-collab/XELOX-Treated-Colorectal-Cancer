"""Script 35: Drug-Target Network Analysis
==========
Phase: Drug-target analysis for 253 XELOX resistance meta-analysis genes.

Pipeline:
  Step 1. DGIdb API query — find known drug-gene interactions
  Step 2. Drug-target network construction & visualization
  Step 3. Drug classification & enrichment (FDA-approved, chemo, targeted)
  Step 4. CMap-style signature analysis (up/down gene lists → drug predictions)
  Step 5. GDSC correlation validation (pathway activity vs drug IC50)

Output:
  results/tables/drug_prediction/
    - dgidb_drug_gene_interactions.csv     — all DGIdb interactions
    - dgidb_summary.csv                     — drugs grouped by mechanism
    - drug_target_network_stats.txt         — summary statistics
    - cmap_drug_predictions.csv             — CMap-like drug enrichment
  results/figures/drug_prediction/
    - drug_target_network.png               — drug-target bipartite graph
    - drug_class_barplot.png                — drug class distribution
    - top_drugs_heatmap.png                 — top predicted drugs
"""

import os, sys, json, time
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import requests
import warnings
warnings.filterwarnings('ignore')

PROJECT = r"/path/to/xelox_project"
TAB_DIR = os.path.join(PROJECT, "results", "tables", "drug_prediction")
FIG_DIR = os.path.join(PROJECT, "results", "figures", "drug_prediction")
os.makedirs(TAB_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

# ============================================================
# Step 0: Load 253 meta-analysis genes
# ============================================================
print("=" * 70)
print("XELOX RESISTANCE — DRUG-TARGET NETWORK ANALYSIS")
print("=" * 70)

gene_df = pd.read_csv(os.path.join(TAB_DIR, "253_meta_genes.csv"))
genes_up = gene_df[gene_df['direction_consensus'] == 'up']['gene_symbol'].tolist()
genes_dn = gene_df[gene_df['direction_consensus'] == 'down']['gene_symbol'].tolist()
all_genes = gene_df['gene_symbol'].tolist()

print(f"\nLoaded {len(all_genes)} meta-analysis genes ({len(genes_up)} up, {len(genes_dn)} down)")

# ============================================================
# Step 1: DGIdb API Query (v4)
# ============================================================
print("\n" + "=" * 70)
print("STEP 1: DGIdb Drug-Gene Interaction Query")
print("=" * 70)

DGIDB_URL = "https://dgidb.org/api/v4/interactions"
# Query in batches of 50 (API limit)
BATCH_SIZE = 50

all_interactions = []

# First, try a simpler query to verify API
print("  Testing API connection...")
try:
    test_payload = {"genes": ["EGFR", "BRAF"], "interactionSources": ["DrugBank"]}
    test_resp = requests.post(DGIDB_URL, json=test_payload, timeout=30)
    print(f"  API test: status={test_resp.status_code}, content={test_resp.text[:200]}")
except Exception as e:
    print(f"  API test failed: {e}")

for batch_start in range(0, len(all_genes), BATCH_SIZE):
    batch = all_genes[batch_start:batch_start + BATCH_SIZE]
    payload = {
        "genes": batch,
        "interactionSources": ["DrugBank", "PharmGKB", "TTD", "GuideToPharmacology",
                               "CIViC", "OncoKB"],
        "interactionTypes": ["inhibitor", "antagonist", "agonist", "activator",
                             "suppressor", "inducer", "ligand", "antibody",
                             "modulator", "substrate"],
        "interactionTrustLevels": ["expert_curated", "high", "medium"]
    }
    
    try:
        resp = requests.post(DGIDB_URL, json=payload, timeout=30)
        if resp.status_code == 200:
            data = resp.json()
            for match in data.get('matchedTerms', []):
                gene_name = match.get('searchTerm', '')
                gene_obj = match.get('gene', {})
                entrez_id = gene_obj.get('entrezId', '')
                
                for interaction in gene_obj.get('interactions', []):
                    drug = interaction.get('drug', {})
                    drug_name = drug.get('name', '')
                    drug_concept_id = drug.get('conceptId', '')
                    drug_types = drug.get('drugTypes', '')
                    
                    # Get interaction types
                    types = interaction.get('interactionTypes', [])
                    interaction_type = ', '.join(types) if types else 'not specified'
                    
                    # Get sources
                    sources = interaction.get('interactionSources', [])
                    
                    # Get approval status from drug types
                    approval = 'approved' if 'approved' in str(drug_types).lower() else 'Unknown'
                    
                    # Get sources/PMIDs
                    sources = [s.get('source_db_name', '') for s in interaction.get('interaction_attributes', [])]
                    
                    # Get drug approval status
                    approval = drug.get('approval_status', 'Unknown')
                    
                    all_interactions.append({
                        'gene': gene_name,
                        'entrez_id': entrez_id,
                        'drug_name': drug_name,
                        'drug_concept_id': drug_concept_id,
                        'interaction_type': interaction_type,
                        'approval_status': approval,
                        'sources': '; '.join(set(sources)) if sources else ''
                    })
            
            time.sleep(0.3)  # Rate limiting
            print(f"  Batch {batch_start//BATCH_SIZE + 1}: {len(batch)} genes queried -> {len(data.get('matched_terms', []))} matches")
        else:
            print(f"  WARNING: API returned {resp.status_code} for batch {batch_start//BATCH_SIZE + 1}")
    except Exception as e:
        print(f"  ERROR: {e} for batch {batch_start//BATCH_SIZE + 1}")

# Save results
if all_interactions:
    df_dgidb = pd.DataFrame(all_interactions)
    df_dgidb.drop_duplicates(inplace=True)
    df_dgidb.to_csv(os.path.join(TAB_DIR, "dgidb_drug_gene_interactions.csv"), index=False)
    print(f"\n  Total drug-gene interactions found: {len(df_dgidb)}")
    print(f"  Unique drugs: {df_dgidb['drug_name'].nunique()}")
    print(f"  Unique genes targeted: {df_dgidb['gene'].nunique()}")
else:
    print("\n  WARNING: No DGIdb interactions found. Creating template file.")
    # Check API or try alternative approach
    df_dgidb = pd.DataFrame()

# ============================================================
# Step 1b: Drug-Target Summary Table
# ============================================================
print("\n" + "=" * 70)
print("STEP 2: Drug Classification & Summary")
print("=" * 70)

if len(df_dgidb) > 0:
    # Count drugs per gene
    gene_drug_count = df_dgidb.groupby('gene').agg(
        n_drugs=('drug_name', 'nunique'),
        drugs=('drug_name', lambda x: ', '.join(sorted(set(x))))
    ).reset_index()
    gene_drug_count = gene_drug_count.merge(
        gene_df[['gene_symbol', 'direction_consensus', 'stouffer_z', 'stouffer_fdr']],
        left_on='gene', right_on='gene_symbol', how='left'
    ).drop(columns=['gene_symbol'])
    gene_drug_count = gene_drug_count.sort_values('n_drugs', ascending=False)
    
    # Count drugs by mechanism
    drugs_by_type = df_dgidb.groupby('interaction_type').agg(
        n_drugs=('drug_name', 'nunique'),
        n_genes=('gene', 'nunique')
    ).reset_index().sort_values('n_drugs', ascending=False)
    
    # Count drugs by approval status
    approval_counts = df_dgidb['approval_status'].value_counts().reset_index()
    approval_counts.columns = ['approval_status', 'count']
    
    # Save summaries
    gene_drug_count.to_csv(os.path.join(TAB_DIR, "gene_drug_target_summary.csv"), index=False)
    drugs_by_type.to_csv(os.path.join(TAB_DIR, "drug_class_summary.csv"), index=False)
    approval_counts.to_csv(os.path.join(TAB_DIR, "drug_approval_status.csv"), index=False)
    
    print(f"\n  Top 10 most-targeted genes:")
    for _, row in gene_drug_count.head(10).iterrows():
        print(f"    {row['gene']:<14s} ({row['direction_consensus']}) -> {int(row['n_drugs'])} drugs")
    
    print(f"\n  Top drug interaction types:")
    for _, row in drugs_by_type.head(10).iterrows():
        print(f"    {row['interaction_type']:<25s} {int(row['n_drugs']):>3d} drugs, {int(row['n_genes'])} genes")
    
    print(f"\n  Approval status:")
    for _, row in approval_counts.iterrows():
        print(f"    {row['approval_status']:<25s} {int(row['count'])} interactions")
    
    # ============================================================
    # Figure: Drug class distribution
    # ============================================================
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left: Top gene targets
    top_genes = gene_drug_count.head(20)
    colors = ['#E64B35' if d == 'up' else '#4DBBD5' for d in top_genes['direction_consensus']]
    axes[0].barh(range(len(top_genes)), top_genes['n_drugs'], color=colors)
    axes[0].set_yticks(range(len(top_genes)))
    axes[0].set_yticklabels(top_genes['gene'], fontsize=9)
    axes[0].set_xlabel('Number of Known Drugs')
    axes[0].set_title('Top 20 Drug-Targetable Genes\n(in 253 meta-analysis set)', fontweight='bold')
    axes[0].axvline(0, color='gray', linestyle='--', linewidth=0.5)
    up_patch = mpatches.Patch(color='#E64B35', label='Up in resistant')
    dn_patch = mpatches.Patch(color='#4DBBD5', label='Down in resistant')
    axes[0].legend(handles=[up_patch, dn_patch], fontsize=9)
    
    # Right: Interaction types
    top_types = drugs_by_type.head(12)
    axes[1].barh(range(len(top_types)), top_types['n_drugs'], color='#7F77DD')
    axes[1].set_yticks(range(len(top_types)))
    axes[1].set_yticklabels(top_types['interaction_type'], fontsize=9)
    axes[1].set_xlabel('Number of Drugs')
    axes[1].set_title('Drug Interaction Types\n(in XELOX resistance gene set)', fontweight='bold')
    
    plt.suptitle('Drug-Gene Interaction Landscape of XELOX Resistance', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'drug_gene_interaction_landscape.png'), dpi=200, bbox_inches='tight')
    plt.close()
    print(f"\n  Saved drug_gene_interaction_landscape.png")

# ============================================================
# Step 3: Chemotherapy-specific drug enrichment
# ============================================================
print("\n" + "=" * 70)
print("STEP 3: XELOX-Relevant Drug Categorization")
print("=" * 70)

if len(df_dgidb) > 0:
    # Categorize drugs relevant to XELOX resistance
    xelox_drug_categories = {
        "Platinum-based": ["oxaliplatin", "cisplatin", "carboplatin", "nedaplatin"],
        "Fluoropyrimidine": ["5-fluorouracil", "capecitabine", "fluorouracil", "tegafur", "S-1"],
        "Topoisomerase Inhibitor": ["irinotecan", "topotecan", "SN-38"],
        "Anti-EGFR": ["cetuximab", "panitumumab", "erlotinib", "gefitinib", "afatinib"],
        "Anti-VEGF": ["bevacizumab", "sorafenib", "sunitinib", "regorafenib", "aflibercept"],
        "Multi-kinase": ["regorafenib", "sorafenib", "sunitinib", "pazopanib", "cabozantinib"],
        "Immunotherapy": ["pembrolizumab", "nivolumab", "ipilimumab", "atezolizumab"],
        "DNA Repair / PARP": ["olaparib", "niraparib", "talazoparib", "rucaparib", "veliparib"],
        "Chemotherapy (other)": ["gemcitabine", "docetaxel", "paclitaxel", "doxorubicin", "mitomycin"],
        "Targeted (other)": ["trametinib", "dabrafenib", "vemurafenib", "lapatinib", "tucatinib"]
    }
    
    drug_category_map = {}
    for cat, drugs in xelox_drug_categories.items():
        for d in drugs:
            drug_category_map[d.lower()] = cat
    
    df_dgidb['drug_lower'] = df_dgidb['drug_name'].str.lower().str.strip()
    df_dgidb['drug_category'] = df_dgidb['drug_lower'].map(drug_category_map).fillna('Other')
    
    # Find XELOX-relevant interactions
    xelox_relevant = df_dgidb[df_dgidb['drug_category'] != 'Other'].copy()
    xelox_relevant = xelox_relevant.groupby(['drug_category', 'drug_name', 'gene']).size().reset_index(name='count')
    
    if len(xelox_relevant) > 0:
        xelox_relevant.to_csv(os.path.join(TAB_DIR, "xelox_relevant_drug_gene_pairs.csv"), index=False)
        print(f"\n  XELOX-relevant drug-gene pairs: {len(xelox_relevant)}")
        
        # Per category summary
        cat_summary = xelox_relevant.groupby('drug_category').agg(
            n_drugs=('drug_name', 'nunique'),
            n_genes=('gene', 'nunique')
        ).reset_index()
        print(f"\n  Drug categories targeting XELOX resistance genes:")
        for _, row in cat_summary.iterrows():
            print(f"    {row['drug_category']:<25s} {int(row['n_drugs']):>3d} drugs -> {int(row['n_genes'])} genes")
        
        # Figure: XELOX-relevant drug heatmap
        pivot = xelox_relevant.pivot_table(
            index='drug_name', columns='gene', values='count',
            aggfunc='sum', fill_value=0
        )
        if pivot.shape[0] > 0 and pivot.shape[1] > 1:
            fig, ax = plt.subplots(figsize=(max(8, pivot.shape[1] * 0.5), max(6, pivot.shape[0] * 0.4)))
            im = ax.imshow(pivot.values, cmap='YlOrRd', aspect='auto')
            ax.set_xticks(range(pivot.shape[1]))
            ax.set_xticklabels(pivot.columns, fontsize=7, rotation=90)
            ax.set_yticks(range(pivot.shape[0]))
            ax.set_yticklabels(pivot.index, fontsize=8)
            ax.set_xlabel('Target Genes', fontsize=10)
            ax.set_ylabel('XELOX-Relevant Drugs', fontsize=10)
            ax.set_title('XELOX Resistance: Drug-Target Matrix\n(253 meta-analysis genes)', fontweight='bold')
            plt.colorbar(im, ax=ax, label='Interaction count')
            plt.tight_layout()
            plt.savefig(os.path.join(FIG_DIR, 'xelox_relevant_drug_target_matrix.png'), dpi=200, bbox_inches='tight')
            plt.close()
            print(f"  Saved xelox_relevant_drug_target_matrix.png")
    else:
        print("\n  No XELOX-relevant drug interactions found in DGIdb results")
else:
    print("\n  Skipping Step 3 (no DGIdb results)")

# ============================================================
# Step 4: Drug signature enrichment (CMap-like)
# ============================================================
print("\n" + "=" * 70)
print("STEP 4: Drug Signature Enrichment Analysis")
print("=" * 70)

if len(df_dgidb) > 0:
    # Build direction-specific drug lists
    up_genes_set = set(genes_up)
    dn_genes_set = set(genes_dn)
    targeted_up_genes = set(df_dgidb[df_dgidb['gene'].isin(up_genes_set)]['gene'])
    targeted_dn_genes = set(df_dgidb[df_dgidb['gene'].isin(dn_genes_set)]['gene'])
    
    # For each drug, compute enrichment in up vs down direction
    drug_direction = []
    for drug, grp in df_dgidb.groupby('drug_name'):
        target_genes = set(grp['gene'])
        n_up = len(target_genes & up_genes_set)
        n_dn = len(target_genes & dn_genes_set)
        total = n_up + n_dn
        if total > 0:
            # If drug targets predominantly down-regulated genes → could reverse resistance
            # If targets up-regulated genes → could mimic resistance
            direction_score = (n_dn - n_up) / total  # positive = reverse resistance
            drug_direction.append({
                'drug_name': drug,
                'n_up_targets': n_up,
                'n_dn_targets': n_dn,
                'total_targets': total,
                'reverse_resistance_score': round(direction_score, 3),
                'primary_categories': ', '.join(sorted(set(grp['drug_category'])))
            })
    
    if drug_direction:
        df_direction = pd.DataFrame(drug_direction)
        df_direction = df_direction.sort_values('reverse_resistance_score', ascending=False)
        df_direction.to_csv(os.path.join(TAB_DIR, "drug_reverse_resistance_scores.csv"), index=False)
        
        print(f"\n  Top 10 drugs predicted to REVERSE XELOX resistance:")
        for _, row in df_direction.head(10).iterrows():
            print(f"    {row['drug_name']:<30s} score={row['reverse_resistance_score']:+.3f} "
                  f"(targets: {int(row['n_dn_targets'])} down + {int(row['n_up_targets'])} up)")
        
        print(f"\n  Top 10 drugs that MIMIC XELOX resistance (negative score):")
        for _, row in df_direction.tail(10).iloc[::-1].iterrows():
            print(f"    {row['drug_name']:<30s} score={row['reverse_resistance_score']:+.3f} "
                  f"(targets: {int(row['n_dn_targets'])} down + {int(row['n_up_targets'])} up)")
        
        # Figure: Top predicted drugs
        fig, ax = plt.subplots(figsize=(10, 8))
        top20 = pd.concat([df_direction.head(10), df_direction.tail(10)])
        top20 = top20.sort_values('reverse_resistance_score')
        colors = ['#4DBBD5' if s > 0 else '#E64B35' for s in top20['reverse_resistance_score']]
        bars = ax.barh(range(len(top20)), top20['reverse_resistance_score'], color=colors)
        ax.set_yticks(range(len(top20)))
        ax.set_yticklabels(top20['drug_name'], fontsize=9)
        ax.axvline(0, color='gray', linestyle='--', linewidth=0.8)
        ax.set_xlabel('Reverse Resistance Score\n(positive = targets down-in-resistant genes, reverses resistance)', fontsize=10)
        ax.set_title('Predicted Drugs: Ability to Reverse XELOX Resistance\n(Based on 253-gene meta-analysis signature)', fontweight='bold')
        rev_patch = mpatches.Patch(color='#4DBBD5', label='May reverse resistance')
        mim_patch = mpatches.Patch(color='#E64B35', label='May mimic resistance')
        ax.legend(handles=[rev_patch, mim_patch], fontsize=9)
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, 'drug_reverse_resistance_scores.png'), dpi=200, bbox_inches='tight')
        plt.close()
        print(f"\n  Saved drug_reverse_resistance_scores.png")
else:
    print("\n  Skipping Step 4 (no DGIdb results)")

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 70)
print("DRUG-TARGET NETWORK ANALYSIS COMPLETE")
print("=" * 70)

summary_lines = [
    "XELOX RESISTANCE — DRUG-TARGET NETWORK ANALYSIS",
    f"Meta-analysis genes analyzed: {len(all_genes)}",
    f"  Up-regulated (risk): {len(genes_up)}",
    f"  Down-regulated (protective): {len(genes_dn)}",
]
if len(df_dgidb) > 0:
    summary_lines.extend([
        f"DGIdb drug-gene interactions: {len(df_dgidb)}",
        f"Unique drugs: {df_dgidb['drug_name'].nunique()}",
        f"Unique targetable genes: {df_dgidb['gene'].nunique()}",
        f"FDA-approved interactions: {int(approval_counts[approval_counts['approval_status']=='approved']['count'].sum()) if 'approved' in approval_counts['approval_status'].values else 0}",
        "",
        "OUTPUT FILES:",
        f"  {TAB_DIR}",
    ])

summary_text = "\n".join(summary_lines)
with open(os.path.join(TAB_DIR, "drug_target_network_stats.txt"), 'w') as f:
    f.write(summary_text)
print(summary_text)
print("\nSCRIPT COMPLETE")
