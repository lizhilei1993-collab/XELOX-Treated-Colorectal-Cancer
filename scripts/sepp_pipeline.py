#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SEPP流水线脚本 - 独立运行
SEA靶点预测 → OpenTargets关联分析 → PHAWES风险注释 → DeepSeek AI预测

用法: python sepp_pipeline.py [--skip-sea] [--skip-ot] [--skip-phewas] [--skip-ai] [--sea-only]
"""

import os, sys, json, time, asyncio, concurrent.futures
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional

# 添加SEPP到Python路径
SEPP_DIR = Path(__file__).parent.parent / "SideEffectPredictorPro_1.2生产版 - 副本"
sys.path.insert(0, str(SEPP_DIR))

# 输出目录
OUTPUT_DIR = Path(__file__).parent.parent / "results" / "sepp_pipeline"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# XELOX化合物定义
# ============================================================
XELOX_COMPOUNDS = {
    "oxaliplatin": {
        "name": "奥沙利铂",
        "smiles": "C1CC(C1)(C(=O)O)C2CC(C2)(C(=O)O).[Pt]",
        "role": "DNA交联剂",
        "category": "chemotherapy"
    },
    "capecitabine": {
        "name": "卡培他滨",
        "smiles": "CCCCCOC(=O)NC1=NC(=O)N(C=C1F)[C@@H]2O[C@H]([C@@H]([C@H]2O)O)CO",
        "role": "5-FU前药",
        "category": "chemotherapy"
    },
    "5-fluorouracil": {
        "name": "5-氟尿嘧啶",
        "smiles": "C1=C(C(=O)NC(=O)N1)F",
        "role": "活性代谢物",
        "category": "chemotherapy"
    },
    "leucovorin": {
        "name": "亚叶酸钙",
        "smiles": "C1C(N(C2=C(N1)N=C(NC2=O)N)CNC3=CC=C(C=C3)C(=O)NC(CCC(=O)O)C(=O)O)CNC4=CC=C(C=C4)C(=O)NC(CCC(=O)O)C(=O)O",
        "role": "5-FU增效剂",
        "category": "chemotherapy"
    },
}


# ============================================================
# 步骤1: SEA靶点预测 (同步调用, crawler内部已处理async)
# ============================================================
def run_sea_prediction(compound_key: str, compound_info: Dict) -> Optional[Dict]:
    """对单个化合物运行SEA靶点预测 (同步)"""
    print(f"\n{'='*60}")
    print(f"[步骤1-SEA] 预测 {compound_info['name']} ({compound_key})")
    print(f"  SMILES: {compound_info['smiles'][:60]}...")
    
    try:
        from modules.sea_crawlee_crawler import SEACrawleeCrawler
        
        data_dir = str(SEPP_DIR / "data" / "sea_predictions")
        crawler = SEACrawleeCrawler(data_dir=data_dir)
        
        result = crawler.crawl_sea_predictions(
            smiles=compound_info['smiles'],
            compound_name=compound_info['name']
        )
        
        # 增强结果
        result['compound_key'] = compound_key
        result['compound_role'] = compound_info['role']
        result['pipeline_timestamp'] = datetime.now().isoformat()
        
        # 过滤---只保留人类基因
        human_targets = []
        for t in result.get('targets', []):
            gene_name = t.get('gene_name', '').strip()
            uniprot = t.get('uniprot_id', '')
            # 保留有基因名且是人类靶点的
            if gene_name and gene_name not in ['', '-']:
                human_targets.append(t)
        result['targets'] = human_targets
        result['total_human_targets'] = len(human_targets)
        
        # 保存
        output_file = OUTPUT_DIR / f"01_SEA_{compound_key}.json"
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(result, f, indent=2, ensure_ascii=False)
        print(f"  [OK] 保存到: {output_file.name}")
        print(f"  [结果] {len(human_targets)} 个人类靶点 (总{len(result.get('targets', [])) + len([t for t in result.get('targets', []) if not t.get('gene_name','').strip()])}个)")
        
        # 列出前10个靶点
        for i, t in enumerate(human_targets[:10]):
            print(f"    {i+1}. {t.get('gene_name', '?')} | {t.get('target_name', '?')} | score={t.get('score', 0):.2f} p={t.get('p_value', 1):.2e}")
        
        return result
    
    except Exception as e:
        import traceback
        print(f"  [FAIL] {e}")
        traceback.print_exc()
        return None


def run_all_sea_predictions() -> Dict:
    """运行所有XELOX化合物的SEA预测"""
    all_results = {}
    
    for key, info in XELOX_COMPOUNDS.items():
        result = run_sea_prediction(key, info)
        if result:
            all_results[key] = result
    
    # 汇总所有靶点
    all_targets = set()
    compound_target_map = {}
    for key, result in all_results.items():
        targets = [t.get('gene_name', '') for t in result.get('targets', []) if t.get('gene_name', '').strip()]
        compound_target_map[key] = sorted(set(targets))
        all_targets.update(targets)
    
    summary = {
        'compounds_analyzed': len(all_results),
        'total_unique_targets': len(all_targets),
        'compound_targets': compound_target_map,
        'all_targets_list': sorted(list(all_targets)),
        'timestamp': datetime.now().isoformat()
    }
    
    summary_file = OUTPUT_DIR / "01_SEA_summary.json"
    with open(summary_file, 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(f"\n{'='*60}")
    print(f"[SEA总结] {len(all_targets)} 个独特靶点基因")
    for key, targets in compound_target_map.items():
        print(f"  {XELOX_COMPOUNDS[key]['name']}: {len(targets)} 个靶点")
    print(f"总结文件: {summary_file}")
    
    return summary


# ============================================================
# 步骤2: OpenTargets关联分析
# ============================================================
def run_opentargets_analysis(all_targets: List[str]) -> Dict:
    """对靶点列表运行OpenTargets基因-疾病关联查询"""
    print(f"\n{'='*60}")
    print(f"[步骤2-OpenTargets] 查询 {len(all_targets)} 个靶点的疾病关联")
    
    import requests
    
    gene_associations = {}
    success_count = 0
    
    for i, gene in enumerate(all_targets):
        if not gene or len(gene) < 2:
            continue
        
        query = f'''
        {{
          search(queryString: "{gene}") {{
            hits {{
              id
              name
              entity
              description
            }}
          }}
        }}
        '''
        
        try:
            r = requests.post(
                'https://api.platform.opentargets.org/api/v4/graphql',
                json={'query': query},
                timeout=20,
                headers={'Content-Type': 'application/json'}
            )
            
            if r.status_code == 200:
                data = r.json()
                hits = data.get('data', {}).get('search', {}).get('hits', [])
                gene_associations[gene] = {
                    'hits_count': len(hits),
                    'hits': [{
                        'id': h.get('id'),
                        'name': h.get('name'),
                        'entity': h.get('entity'),
                        'description': h.get('description')
                    } for h in hits[:20]]  # 限制Top 20
                }
                success_count += 1
                print(f"  [{i+1}/{len(all_targets)}] {gene}: {len(hits)} associations")
            else:
                gene_associations[gene] = {'hits_count': 0, 'error': f'HTTP {r.status_code}'}
                print(f"  [{i+1}/{len(all_targets)}] {gene}: ERROR {r.status_code}")
        
        except Exception as e:
            gene_associations[gene] = {'hits_count': 0, 'error': str(e)}
            print(f"  [{i+1}/{len(all_targets)}] {gene}: FAILED")
        
        time.sleep(0.3)  # 速率限制
    
    output_file = OUTPUT_DIR / "02_OpenTargets_associations.json"
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(gene_associations, f, indent=2, ensure_ascii=False)
    
    print(f"\n[OpenTargets总结] {success_count}/{len(all_targets)} 成功, 保存到 {output_file.name}")
    return gene_associations


# ============================================================
# 步骤3: PHAWES风险注释 (需要单独的event loop)
# ============================================================
async def _phewas_single_gene(gene: str, i: int, total: int) -> Dict:
    """单个基因的PHAWES查询"""
    from modules.playwright_azphewas_crawler import AzphewasPlaywrightCrawler
    
    try:
        async with AzphewasPlaywrightCrawler(headless=True) as crawler:
            result = await crawler.crawl_gene_data(gene)
            binary_count = len(result.get('binary_traits', []))
            cont_count = len(result.get('continuous_traits', []))
            print(f"  [{i}/{total}] {gene}: binary={binary_count}, continuous={cont_count}")
            return {'gene': gene, 'result': result}
    except Exception as e:
        print(f"  [{i}/{total}] {gene}: FAILED - {e}")
        return {'gene': gene, 'error': str(e)}


async def _phewas_all_async(all_targets: List[str]) -> Dict:
    """异步运行所有基因的PHAWES"""
    results = {}
    for i, gene in enumerate(all_targets):
        if not gene or len(gene) < 2:
            continue
        r = await _phewas_single_gene(gene, i+1, len(all_targets))
        if 'error' in r:
            results[gene] = {'error': r['error']}
        else:
            results[gene] = r['result']
    return results


def run_phewas_annotation(all_targets: List[str]) -> Dict:
    """对靶点列表运行PHAWES全基因组风险注释"""
    print(f"\n{'='*60}")
    print(f"[步骤3-PHAWES] 注释 {len(all_targets)} 个靶点的基因-表型关联")
    
    phewas_results = asyncio.run(_phewas_all_async(all_targets))
    
    success_count = sum(1 for v in phewas_results.values() if 'error' not in v)
    
    output_file = OUTPUT_DIR / "03_PHAWES_annotations.json"
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(phewas_results, f, indent=2, ensure_ascii=False)
    
    print(f"\n[PHAWES总结] {success_count}/{len(all_targets)} 成功, 保存到 {output_file.name}")
    return phewas_results


# ============================================================
# 步骤4: DeepSeek AI智能预测
# ============================================================
def run_deepseek_prediction(sea_summary: Dict, ot_results: Dict, phewas_results: Dict) -> Dict:
    """整合前三步数据, 调用DeepSeek API生成综合ADR预测"""
    print(f"\n{'='*60}")
    print(f"[步骤4-DeepSeek] AI综合分析")
    
    import requests
    
    api_key = os.environ.get("DEEPSEEK_API_KEY", "")
    if not api_key:
        print("  [SKIP] DeepSeek API Key未配置")
        return {'status': 'skipped', 'reason': 'No API key'}
    
    target_list = sea_summary.get('all_targets_list', [])
    compound_targets = sea_summary.get('compound_targets', {})
    
    system_prompt = """你是一位资深的临床药理学和肿瘤药物安全性专家。

请基于提供的XELOX化疗方案（奥沙利铂+卡培他滨+5-FU+亚叶酸钙）靶点预测数据，进行综合分析。

分析要求:
1. 识别哪些SEA预测靶点已知与XELOX耐药相关（查阅你的知识库）
2. 识别哪些靶点已知与XELOX典型ADR（神经毒性、骨髓抑制、手足综合征、消化道毒性）相关
3. 分析"耐药"和"毒性"是否有共享的分子通路和驱动基因
4. 给出最有临床转化价值的TOP5靶点
5. 评估整体风险等级

输出要求：严格的JSON格式，包含以下字段：
{
  "overall_risk_assessment": "string",
  "resistance_associated_targets": [{"gene": "string", "mechanism": "string"}],
  "adr_associated_targets": [{"gene": "string", "adr_type": "string", "mechanism": "string"}],
  "shared_resistance_toxicity_pathways": ["string"],
  "top5_translational_targets": [{"gene": "string", "rationale": "string"}],
  "clinical_recommendations": ["string"]
}"""
    
    user_prompt = f"""## XELOX化疗方案靶点分析数据

### 化合物-靶点对应关系（SEA预测）
{json.dumps(compound_targets, ensure_ascii=False, indent=2)}

### 全部SEA预测靶点列表（{len(target_list)}个独特基因）
{', '.join(target_list[:50])}...
(共{len(target_list)}个基因)

### OpenTargets疾病关联
{sum(1 for v in ot_results.values() if 'error' not in v)} 个基因在OpenTargets中有已知疾病关联

请基于以上数据和你的知识库进行综合分析。"""
    
    try:
        r = requests.post(
            'https://api.deepseek.com/v1/chat/completions',
            headers={
                'Content-Type': 'application/json',
                'Authorization': f'Bearer {api_key}'
            },
            json={
                'model': 'deepseek-chat',
                'messages': [
                    {'role': 'system', 'content': system_prompt},
                    {'role': 'user', 'content': user_prompt}
                ],
                'max_tokens': 4096,
                'temperature': 0
            },
            timeout=120
        )
        
        if r.status_code == 200:
            result = r.json()
            ai_response = result['choices'][0]['message']['content']
            
            output = {
                'timestamp': datetime.now().isoformat(),
                'model': 'deepseek-chat',
                'analysis': ai_response,
                'usage': result.get('usage', {})
            }
            
            output_file = OUTPUT_DIR / "04_DeepSeek_AI_analysis.json"
            with open(output_file, 'w', encoding='utf-8') as f:
                json.dump(output, f, indent=2, ensure_ascii=False)
            
            print(f"  [OK] 分析完成, 字数: {len(ai_response)}, 保存到 {output_file.name}")
            return output
        else:
            print(f"  [FAIL] HTTP {r.status_code}")
            return {'status': 'error', 'code': r.status_code, 'response': r.text[:200]}
    
    except Exception as e:
        print(f"  [FAIL] {e}")
        return {'status': 'error', 'error': str(e)}


# ============================================================
# 主流程
# ============================================================
def main():
    import argparse
    parser = argparse.ArgumentParser(description='SEPP Pipeline for XELOX')
    parser.add_argument('--skip-sea', action='store_true')
    parser.add_argument('--skip-ot', action='store_true')
    parser.add_argument('--skip-phewas', action='store_true')
    parser.add_argument('--skip-ai', action='store_true')
    parser.add_argument('--sea-only', action='store_true')
    parser.add_argument('--targets-only', type=int, help='Limit targets for OT/Phewas to top N')
    args = parser.parse_args()
    
    print("=" * 60)
    print("SEPP流水线 - XELOX化疗方案不良反应预测")
    print(f"时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"输出目录: {OUTPUT_DIR}")
    print("=" * 60)
    
    all_targets = []
    sea_summary = None
    
    # ---- 步骤1: SEA ----
    if not args.skip_sea:
        sea_summary = run_all_sea_predictions()
        all_targets = sea_summary.get('all_targets_list', [])
    else:
        summary_file = OUTPUT_DIR / "01_SEA_summary.json"
        if summary_file.exists():
            with open(summary_file, 'r', encoding='utf-8') as f:
                sea_summary = json.load(f)
            all_targets = sea_summary.get('all_targets_list', [])
            print(f"[加载] 从缓存加载 {len(all_targets)} 个靶点")
    
    if args.sea_only:
        print("\n[完成] 仅运行SEA预测")
        return
    
    if not all_targets:
        print("[ERROR] 无靶点数据, 终止")
        return
    
    # 可选: 限制靶点数量
    if args.targets_only and args.targets_only < len(all_targets):
        all_targets = all_targets[:args.targets_only]
        print(f"[限制] 仅分析前 {args.targets_only} 个靶点")
    
    # ---- 步骤2: OpenTargets ----
    ot_results = {}
    if not args.skip_ot:
        ot_results = run_opentargets_analysis(all_targets)
    
    # ---- 步骤3: PHAWES ----
    phewas_results = {}
    if not args.skip_phewas:
        phewas_results = run_phewas_annotation(all_targets)
    
    # ---- 步骤4: DeepSeek ----
    if not args.skip_ai and sea_summary:
        run_deepseek_prediction(sea_summary, ot_results, phewas_results)
    
    print(f"\n{'=' * 60}")
    print(f"SEPP流水线完成!")
    print(f"输出目录: {OUTPUT_DIR}")
    print(f"文件列表:")
    for f in sorted(OUTPUT_DIR.glob("*.json")):
        size_kb = f.stat().st_size / 1024
        print(f"  {f.name} ({size_kb:.1f} KB)")
    print("=" * 60)


if __name__ == '__main__':
    main()
