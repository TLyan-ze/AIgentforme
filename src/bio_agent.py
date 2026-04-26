#!/usr/bin/env python3
"""
BioAgent — 生物信息学分析智能体 v2.0
====================================
多任务生信分析 + LLM 报告生成平台

任务:
  1. rnaseq   — RNA-seq 差异表达 + GO/KEGG 富集 (从表达矩阵到报告)
  2. variant  — 变异检测报告 (VCF → 统计 + 可视化 + 报告)
  3. scrna    — 单细胞分析报告 (标记基因 → 细胞类型解读)
  4. enrich   — 独立功能富集分析 (基因列表 → GO/KEGG/Reactome)

作者: 颜泽钦
"""

import os
import sys
import argparse
import logging
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)

# 确保 src/ 在 import 路径中
sys.path.insert(0, str(Path(__file__).parent))


def _load_env():
    """加载 .env 文件。"""
    env_file = Path(__file__).parent.parent / ".env"
    if env_file.exists():
        with open(env_file) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#") and "=" in line:
                    k, v = line.split("=", 1)
                    os.environ.setdefault(k.strip(), v.strip())


def cmd_rnaseq(args):
    """Task 1: RNA-seq 差异表达 + 富集分析。"""
    from task_rnaseq import run_pipeline
    from generate_data import generate_rnaseq_data

    if args.sample_data:
        data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
        counts_path, meta_path = generate_rnaseq_data(
            data_dir, species=args.species, seed=42)
    else:
        counts_path = args.counts
        meta_path = args.metadata
        if not counts_path or not meta_path:
            print("[错误] 请指定 --counts 和 --metadata，或使用 --sample-data")
            sys.exit(1)

    result = run_pipeline(
        counts_path=counts_path,
        metadata_path=meta_path,
        output_dir=args.output,
        project_name=args.name,
        species=args.species,
        fc_threshold=args.fc,
        pval_threshold=args.pval,
        condition_col=args.condition,
    )

    _print_result("RNA-seq DEG", result)


def cmd_variant(args):
    """Task 2: 变异检测报告。"""
    from task_variant import run_variant_pipeline
    from generate_data import generate_vcf_data

    if args.sample_data:
        data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
        vcf_path = generate_vcf_data(data_dir)
    else:
        vcf_path = args.vcf
        if not vcf_path:
            print("[错误] 请指定 --vcf 或使用 --sample-data")
            sys.exit(1)

    result = run_variant_pipeline(vcf_path, args.output, args.name)
    _print_result("Variant Calling", result)


def cmd_scrna(args):
    """Task 3: 单细胞分析报告。"""
    from task_scrna import run_scrna_pipeline
    from generate_data import generate_scrna_data

    if args.sample_data:
        data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
        marker_path, summary_path = generate_scrna_data(data_dir, species=args.species)
    else:
        marker_path = args.markers
        summary_path = args.summary
        if not marker_path or not summary_path:
            print("[错误] 请指定 --markers 和 --summary，或使用 --sample-data")
            sys.exit(1)

    result = run_scrna_pipeline(marker_path, summary_path, args.output,
                                args.name, args.species)
    _print_result("scRNA-seq", result)


def cmd_enrich(args):
    """Task 4: 独立功能富集分析。"""
    from task_enrichment import run_enrichment_pipeline
    from generate_data import generate_gene_list

    if args.sample_data:
        data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
        gene_path = generate_gene_list(data_dir)
    else:
        gene_path = args.genes
        if not gene_path:
            print("[错误] 请指定 --genes 或使用 --sample-data")
            sys.exit(1)

    result = run_enrichment_pipeline(gene_path, args.output, args.name, args.species)
    _print_result("Enrichment", result)


def _print_result(task_name: str, result: dict):
    """打印结果摘要。"""
    print(f"\n{'='*60}")
    print(f"[{task_name}] 分析完成!")
    print(f"  报告: {result.get('report', 'N/A')}")
    if "figures" in result:
        print(f"  图表: {len(result['figures'])} 张")
    if "stats" in result:
        s = result["stats"]
        if "sig_genes" in s:
            print(f"  显著差异基因: {s['sig_genes']} (上调 {s['up_genes']}, 下调 {s['down_genes']})")
        elif "total_variants" in s:
            print(f"  变异数: {s['total_variants']} (SNV {s['snv_count']}, InDel {s['indel_count']})")
        elif "n_cell_types" in s:
            print(f"  细胞类型: {s['n_cell_types']} | 总细胞: {s['total_cells']}")
    print(f"  耗时: {result.get('elapsed', 0):.1f} 秒")
    print(f"{'='*60}")


def main():
    parser = argparse.ArgumentParser(
        description="BioAgent — 生物信息学分析智能体 v2.0",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # RNA-seq 分析 (示例数据)
  python bio_agent.py rnaseq --sample-data --species human

  # RNA-seq 分析 (真实数据)
  python bio_agent.py rnaseq -c counts.csv -m metadata.csv --species mouse

  # 变异检测
  python bio_agent.py variant --sample-data

  # 单细胞分析
  python bio_agent.py scrna --sample-data --species human

  # 功能富集
  python bio_agent.py enrich --sample-data
        """,
    )
    subparsers = parser.add_subparsers(dest="task", help="选择任务")

    # ── Task 1: RNA-seq ──
    p1 = subparsers.add_parser("rnaseq", help="RNA-seq 差异表达 + GO/KEGG 富集")
    p1.add_argument("-c", "--counts", help="表达矩阵 CSV (基因 × 样本)")
    p1.add_argument("-m", "--metadata", help="样本信息 CSV (sample_id, condition)")
    p1.add_argument("-s", "--species", choices=["human", "mouse"], default="human")
    p1.add_argument("-o", "--output", default="output")
    p1.add_argument("-n", "--name", default="RNA-seq 差异表达分析")
    p1.add_argument("--fc", type=float, default=1.0, help="log2FC 阈值")
    p1.add_argument("--pval", type=float, default=0.05, help="padj 阈值")
    p1.add_argument("--condition", default="condition", help="条件列名")
    p1.add_argument("--sample-data", action="store_true", help="使用示例数据")

    # ── Task 2: Variant ──
    p2 = subparsers.add_parser("variant", help="变异检测报告")
    p2.add_argument("--vcf", help="VCF 文件路径")
    p2.add_argument("-o", "--output", default="output")
    p2.add_argument("-n", "--name", default="变异检测分析")
    p2.add_argument("--sample-data", action="store_true")

    # ── Task 3: scRNA ──
    p3 = subparsers.add_parser("scrna", help="单细胞分析报告")
    p3.add_argument("--markers", help="标记基因 CSV")
    p3.add_argument("--summary", help="细胞类型汇总 CSV")
    p3.add_argument("-s", "--species", choices=["human", "mouse"], default="human")
    p3.add_argument("-o", "--output", default="output")
    p3.add_argument("-n", "--name", default="单细胞分析")
    p3.add_argument("--sample-data", action="store_true")

    # ── Task 4: Enrichment ──
    p4 = subparsers.add_parser("enrich", help="独立功能富集分析")
    p4.add_argument("--genes", help="基因列表文件 (每行一个基因)")
    p4.add_argument("-s", "--species", choices=["human", "mouse"], default="human")
    p4.add_argument("-o", "--output", default="output")
    p4.add_argument("-n", "--name", default="功能富集分析")
    p4.add_argument("--sample-data", action="store_true")

    args = parser.parse_args()

    if not args.task:
        parser.print_help()
        sys.exit(0)

    _load_env()

    dispatch = {
        "rnaseq": cmd_rnaseq,
        "variant": cmd_variant,
        "scrna": cmd_scrna,
        "enrich": cmd_enrich,
    }
    dispatch[args.task](args)


if __name__ == "__main__":
    main()
