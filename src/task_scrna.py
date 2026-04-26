#!/usr/bin/env python3
"""
Task 3: 单细胞 RNA-seq 分析报告智能体
从标记基因结果 + 细胞类型 → 统计分析 → 可视化 → LLM 报告
"""

import os
import json
import logging
from datetime import datetime

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams["font.sans-serif"] = ["Arial Unicode MS", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

logger = logging.getLogger("BioAgent.scRNA")


def load_scrna_data(marker_path: str, summary_path: str):
    """加载单细胞标记基因结果和细胞类型汇总。"""
    logger.info(f"加载标记基因: {marker_path}")
    markers = pd.read_csv(marker_path)
    logger.info(f"加载细胞汇总: {summary_path}")
    summary = pd.read_csv(summary_path)
    return markers, summary


def scrna_statistics(markers: pd.DataFrame, summary: pd.DataFrame) -> dict:
    """计算单细胞分析统计。"""
    cell_types = summary["cell_type"].tolist()
    n_cells_total = int(summary["n_cells"].sum())

    # 每个 cluster 的标记基因数
    marker_counts = markers.groupby("cluster").size().to_dict()

    # 高质量标记基因 (avg_log2FC > 1)
    high_quality = markers[markers["avg_log2FC"] > 1.0]
    hq_per_cluster = high_quality.groupby("cluster").size().to_dict()

    # 共享标记基因（出现在多个 cluster）
    gene_cluster_map = markers.groupby("gene")["cluster"].apply(set)
    shared_genes = {g: list(cs) for g, cs in gene_cluster_map.items() if len(cs) > 1}

    return {
        "n_cell_types": len(cell_types),
        "total_cells": n_cells_total,
        "cell_types": cell_types,
        "cells_per_type": dict(zip(summary["cell_type"], summary["n_cells"])),
        "total_markers": len(markers),
        "markers_per_cluster": marker_counts,
        "high_quality_markers": len(high_quality),
        "hq_per_cluster": hq_per_cluster,
        "shared_marker_genes": len(shared_genes),
        "top_shared": dict(list(shared_genes.items())[:10]),
    }


def plot_scrna_stats(markers: pd.DataFrame, summary: pd.DataFrame,
                     stats: dict, output_dir: str) -> list:
    """绘制单细胞分析图表。"""
    paths = []

    # 1. 细胞类型分布
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.barh(summary["cell_type"], summary["n_cells"], color="#3C5488", alpha=0.8)
    ax.set_xlabel("Cell Count")
    ax.set_title("Cell Type Distribution")
    ax.invert_yaxis()
    path = os.path.join(output_dir, "cell_type_distribution.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)

    # 2. 每个 cluster 的标记基因数
    fig, ax = plt.subplots(figsize=(10, 6))
    clusters = list(stats["markers_per_cluster"].keys())
    counts = list(stats["markers_per_cluster"].values())
    ax.bar(clusters, counts, color="#00A087", alpha=0.8)
    ax.set_ylabel("Marker Gene Count")
    ax.set_title("Marker Genes per Cell Type")
    ax.tick_params(axis="x", rotation=45)
    path = os.path.join(output_dir, "markers_per_cluster.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)

    # 3. Top 标记基因 log2FC 热图
    top_list = []
    for ct, grp in markers.groupby("cluster"):
        top_list.append(grp.nlargest(5, "avg_log2FC"))
    top_per_cluster = pd.concat(top_list, ignore_index=True)
    pivot = top_per_cluster.pivot_table(index="gene", columns="cluster",
                                         values="avg_log2FC", fill_value=0)

    fig, ax = plt.subplots(figsize=(10, max(8, len(pivot) * 0.4)))
    sns.heatmap(pivot, cmap="Reds", annot=True, fmt=".1f", ax=ax, linewidths=0.3)
    ax.set_title("Top Marker Genes (log2FC)")
    path = os.path.join(output_dir, "marker_heatmap.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)

    return paths


def generate_scrna_report(stats: dict, figures: list, output_dir: str,
                          project_name: str, species: str) -> str:
    """生成单细胞分析报告。"""
    from task_rnaseq import _call_llm

    sys_prompt = "你是资深单细胞分析专家，请根据结果撰写中文报告。"

    cells_str = "\n".join(f"  {ct}: {n}" for ct, n in stats["cells_per_type"].items())
    user_prompt = f"""单细胞分析结果（物种: {species}）：

细胞类型数: {stats['n_cell_types']}
总细胞数: {stats['total_cells']}
标记基因总数: {stats['total_markers']}
高质量标记基因: {stats['high_quality_markers']}

各细胞类型细胞数:
{cells_str}

请撰写单细胞分析报告，包括细胞类型注释讨论和标记基因解读。"""

    llm_text = _call_llm(sys_prompt, user_prompt)
    if not llm_text:
        llm_text = (f"共鉴定 {stats['n_cell_types']} 种细胞类型，"
                    f"总计 {stats['total_cells']} 个细胞。"
                    f"发现 {stats['total_markers']} 个标记基因。")

    fig_md = "\n".join(f"![{os.path.basename(p)}](../{p})" for p in figures)
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    report = f"""# {project_name} — 单细胞 RNA-seq 分析报告

> **生成时间**: {now} | **物种**: {species} | **工具**: BioAgent v2.0

## 分析统计

| 指标 | 值 |
|------|------|
| 细胞类型数 | {stats['n_cell_types']} |
| 总细胞数 | {stats['total_cells']} |
| 标记基因总数 | {stats['total_markers']} |
| 高质量标记基因 | {stats['high_quality_markers']} |

## 细胞类型详情

| 细胞类型 | 细胞数 | 标记基因数 |
|----------|--------|-----------|
"""
    for ct in stats["cell_types"]:
        n_cells = stats["cells_per_type"].get(ct, 0)
        n_markers = stats["markers_per_cluster"].get(ct, 0)
        report += f"| {ct} | {n_cells} | {n_markers} |\n"

    report += f"""
## 可视化

{fig_md}

## 分析解读

{llm_text}

---
*BioAgent 自动生成 — {now}*
"""

    path = os.path.join(output_dir, "scrna_report.md")
    with open(path, "w", encoding="utf-8") as f:
        f.write(report)
    return path


def run_scrna_pipeline(marker_path: str, summary_path: str,
                       output_dir: str = "output",
                       project_name: str = "单细胞分析",
                       species: str = "human") -> dict:
    """运行单细胞分析报告 Pipeline。"""
    start = datetime.now()
    os.makedirs(output_dir, exist_ok=True)

    markers, summary = load_scrna_data(marker_path, summary_path)
    stats = scrna_statistics(markers, summary)
    figures = plot_scrna_stats(markers, summary, stats, output_dir)
    report_path = generate_scrna_report(stats, figures, output_dir, project_name, species)

    with open(os.path.join(output_dir, "scrna_stats.json"), "w") as f:
        json.dump(stats, f, ensure_ascii=False, indent=2, default=str)

    elapsed = (datetime.now() - start).total_seconds()
    logger.info(f"单细胞 Pipeline 完成 | {elapsed:.1f}s")
    return {"report": report_path, "figures": figures, "stats": stats, "elapsed": elapsed}
