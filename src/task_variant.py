#!/usr/bin/env python3
"""
Task 2: 变异检测 (Variant Calling) 报告智能体
从 VCF 文件 → 变异分类统计 → 功能影响注释 → LLM 报告
"""

import os
import json
import logging
from datetime import datetime
from collections import Counter

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams["font.sans-serif"] = ["Arial Unicode MS", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

logger = logging.getLogger("BioAgent.Variant")


def parse_vcf(vcf_path: str) -> pd.DataFrame:
    """解析 VCF 文件，提取变异信息和 INFO 字段。"""
    logger.info(f"解析 VCF: {vcf_path}")
    records = []

    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue

            chrom, pos, vid, ref, alt, qual, filt, info = parts[:8]
            info_dict = {}
            for item in info.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info_dict[k] = v
                else:
                    info_dict[item] = True

            records.append({
                "chrom": chrom, "pos": int(pos), "id": vid,
                "ref": ref, "alt": alt, "qual": float(qual) if qual != "." else 0,
                "filter": filt,
                "dp": int(info_dict.get("DP", 0)),
                "af": float(info_dict.get("AF", 0)),
                "type": info_dict.get("TYPE", "unknown"),
                "impact": info_dict.get("IMPACT", "unknown"),
                "gene": info_dict.get("GENE", ""),
            })

    df = pd.DataFrame(records)
    logger.info(f"解析完成: {len(df)} 个变异位点")
    return df


def variant_statistics(vcf_df: pd.DataFrame) -> dict:
    """计算变异统计指标。"""
    type_counts = Counter(vcf_df["type"])
    impact_counts = Counter(vcf_df["impact"])
    chrom_counts = Counter(vcf_df["chrom"])

    # 按基因统计
    gene_counts = vcf_df[vcf_df["gene"] != ""]["gene"].value_counts().head(20)

    return {
        "total_variants": len(vcf_df),
        "snv_count": type_counts.get("SNV", 0),
        "indel_count": type_counts.get("InDel", 0),
        "type_distribution": dict(type_counts),
        "impact_distribution": dict(impact_counts),
        "high_impact": int(impact_counts.get("HIGH", 0)),
        "moderate_impact": int(impact_counts.get("MODERATE", 0)),
        "mean_dp": round(float(vcf_df["dp"].mean()), 1),
        "mean_af": round(float(vcf_df["af"].mean()), 4),
        "qual_range": (round(float(vcf_df["qual"].min()), 1),
                       round(float(vcf_df["qual"].max()), 1)),
        "top_genes": gene_counts.to_dict(),
        "chrom_distribution": dict(chrom_counts),
    }


def plot_variant_stats(vcf_df: pd.DataFrame, stats: dict, output_dir: str) -> list:
    """绘制变异统计图表。"""
    paths = []

    # 1. 变异类型分布
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    type_counts = pd.Series(stats["type_distribution"])
    axes[0].bar(type_counts.index, type_counts.values, color=["#3C5488", "#E64B35"])
    axes[0].set_title("Variant Type Distribution")
    axes[0].set_ylabel("Count")

    impact_counts = pd.Series(stats["impact_distribution"])
    colors = {"HIGH": "#E64B35", "MODERATE": "#F39B7F", "LOW": "#91D1C2", "MODIFIER": "#8491B4B2"}
    axes[1].bar(impact_counts.index, impact_counts.values,
                color=[colors.get(x, "#3C5488") for x in impact_counts.index])
    axes[1].set_title("Functional Impact Distribution")
    axes[1].set_ylabel("Count")

    path = os.path.join(output_dir, "variant_types.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)

    # 2. 染色体分布
    fig, ax = plt.subplots(figsize=(12, 5))
    chrom_series = pd.Series(stats["chrom_distribution"]).sort_index()
    ax.bar(chrom_series.index, chrom_series.values, color="#3C5488", alpha=0.8)
    ax.set_title("Variants per Chromosome")
    ax.set_ylabel("Count")
    ax.tick_params(axis="x", rotation=45)

    path = os.path.join(output_dir, "variant_chromosomes.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)

    # 3. AF / DP 分布
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    axes[0].hist(vcf_df["af"], bins=50, color="#00A087", alpha=0.7, edgecolor="white")
    axes[0].set_title("Allele Frequency Distribution")
    axes[0].set_xlabel("AF"); axes[0].set_ylabel("Count")

    axes[1].hist(vcf_df["dp"], bins=50, color="#E64B35", alpha=0.7, edgecolor="white")
    axes[1].set_title("Read Depth Distribution")
    axes[1].set_xlabel("DP"); axes[1].set_ylabel("Count")

    path = os.path.join(output_dir, "variant_af_dp.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)

    return paths


def generate_variant_report(stats: dict, figures: list, output_dir: str,
                            project_name: str) -> str:
    """生成变异检测报告。"""
    from task_rnaseq import _call_llm

    sys_prompt = "你是资深基因组学分析师，请根据变异检测结果撰写中文报告。"

    top_genes_str = "\n".join(f"  {g}: {c} 个变异" for g, c in list(stats["top_genes"].items())[:10])
    user_prompt = f"""变异检测统计结果：

总变异数: {stats['total_variants']}
SNV: {stats['snv_count']}, InDel: {stats['indel_count']}
高影响: {stats['high_impact']}, 中等影响: {stats['moderate_impact']}
平均深度: {stats['mean_dp']}, 平均 AF: {stats['mean_af']}

Top 变异基因:
{top_genes_str}

请简要解读这些变异检测结果。"""

    llm_text = _call_llm(sys_prompt, user_prompt)
    if not llm_text:
        llm_text = (f"共检测到 {stats['total_variants']} 个变异位点，"
                    f"其中 SNV {stats['snv_count']} 个，InDel {stats['indel_count']} 个。"
                    f"高功能影响变异 {stats['high_impact']} 个，建议重点关注。")

    fig_md = "\n".join(f"![{os.path.basename(p)}](../{p})" for p in figures)
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    report = f"""# {project_name} — 变异检测分析报告

> **生成时间**: {now} | **工具**: BioAgent v2.0

## 变异统计

| 指标 | 值 |
|------|------|
| 总变异数 | {stats['total_variants']} |
| SNV | {stats['snv_count']} |
| InDel | {stats['indel_count']} |
| HIGH impact | {stats['high_impact']} |
| MODERATE impact | {stats['moderate_impact']} |
| 平均深度 | {stats['mean_dp']} |
| 平均 AF | {stats['mean_af']} |

## 可视化

{fig_md}

## 分析解读

{llm_text}

## Top 变异基因

| 基因 | 变异数 |
|------|--------|
"""
    for g, c in stats["top_genes"].items():
        report += f"| {g} | {c} |\n"

    report += f"\n---\n*BioAgent 自动生成 — {now}*\n"

    path = os.path.join(output_dir, "variant_report.md")
    with open(path, "w", encoding="utf-8") as f:
        f.write(report)
    return path


def run_variant_pipeline(vcf_path: str, output_dir: str = "output",
                         project_name: str = "变异检测分析") -> dict:
    """运行变异检测报告 Pipeline。"""
    start = datetime.now()
    os.makedirs(output_dir, exist_ok=True)

    vcf_df = parse_vcf(vcf_path)
    stats = variant_statistics(vcf_df)
    figures = plot_variant_stats(vcf_df, stats, output_dir)
    report_path = generate_variant_report(stats, figures, output_dir, project_name)

    with open(os.path.join(output_dir, "variant_stats.json"), "w") as f:
        json.dump(stats, f, ensure_ascii=False, indent=2, default=str)

    elapsed = (datetime.now() - start).total_seconds()
    logger.info(f"变异检测 Pipeline 完成 | {elapsed:.1f}s")
    return {"report": report_path, "figures": figures, "stats": stats, "elapsed": elapsed}
