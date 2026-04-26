#!/usr/bin/env python3
"""
Task 4: 独立功能富集分析智能体
基因列表 → GO/KEGG/Reactome 富集 → 可视化 → LLM 报告
"""

import os
import json
import logging
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["Arial Unicode MS", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

logger = logging.getLogger("BioAgent.Enrichment")


def load_gene_list(path: str) -> list:
    """加载基因列表（每行一个基因名）。"""
    with open(path) as f:
        genes = [line.strip() for line in f if line.strip()]
    logger.info(f"加载基因列表: {len(genes)} 基因")
    return genes


def run_enrichment(gene_list: list, species: str = "human") -> dict:
    """执行 GO/KEGG/Reactome 富集分析。"""
    logger.info(f"执行富集分析 (物种: {species}, 基因数: {len(gene_list)})...")

    try:
        import gseapy as gp

        if species == "human":
            gene_sets = ["GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022"]
        else:
            gene_sets = ["GO_Biological_Process_2023", "KEGG_2021_Mouse", "Reactome_2022"]

        results = {}
        for gs in gene_sets:
            try:
                enr = gp.enrichr(gene_list=gene_list, gene_sets=[gs],
                                 organism=species, outdir=None, no_plot=True)
                if enr.results is not None and not enr.results.empty:
                    results[gs] = enr.results
                    logger.info(f"  {gs}: {len(enr.results)} 条通路")
                else:
                    results[gs] = pd.DataFrame()
            except Exception as e:
                logger.warning(f"  {gs} 失败: {e}")
                results[gs] = pd.DataFrame()

        return results

    except ImportError:
        logger.warning("gseapy 未安装，跳过富集分析")
        return {}


def plot_enrichment_results(enr_results: dict, output_dir: str) -> list:
    """绘制富集分析图表。"""
    paths = []
    for gs_name, df in enr_results.items():
        if df is None or df.empty:
            continue

        term_col = "Term" if "Term" in df.columns else df.columns[0]
        pval_col = "Adjusted P-value" if "Adjusted P-value" in df.columns else "P-value"
        if pval_col not in df.columns:
            pval_col = df.columns[1]

        top = df.head(20).copy()
        top["-log10(padj)"] = -np.log10(top[pval_col].clip(lower=1e-50))

        fig, ax = plt.subplots(figsize=(10, 8))
        ax.barh(range(len(top)), top["-log10(padj)"], color="#3C5488", alpha=0.8)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top[term_col].str[:60], fontsize=8)
        ax.set_xlabel("-log10(adjusted p-value)")
        label = gs_name.replace("_2023", "").replace("_2021", "").replace("_2022", "")
        ax.set_title(f"{label} Enrichment (Top 20)", fontsize=14)
        ax.invert_yaxis()

        path = os.path.join(output_dir, f"enrichment_{gs_name.split('_')[0].lower()}.png")
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        paths.append(path)
        logger.info(f"富集图 [{gs_name}]: {path}")

    return paths


def generate_enrichment_report(gene_list: list, enr_results: dict, figures: list,
                               output_dir: str, project_name: str,
                               species: str) -> str:
    """生成富集分析报告。"""
    from task_rnaseq import _call_llm

    # 构建通路文本
    pathway_text = ""
    for gs_name, df in enr_results.items():
        if df is None or df.empty:
            continue
        term_col = "Term" if "Term" in df.columns else df.columns[0]
        pval_col = "Adjusted P-value" if "Adjusted P-value" in df.columns else "P-value"
        if pval_col not in df.columns:
            pval_col = df.columns[1]
        pathway_text += f"\n### {gs_name}\n"
        for _, r in df.head(8).iterrows():
            pathway_text += f"- {r[term_col][:70]} (padj={r[pval_col]:.2e})\n"

    sys_prompt = "你是资深生物信息学分析师，请根据功能富集分析结果撰写中文报告。"
    user_prompt = f"""基因功能富集分析结果（物种: {species}）：
输入基因数: {len(gene_list)}

{pathway_text}

请撰写富集分析报告，讨论最显著的生物学通路及其意义。"""

    llm_text = _call_llm(sys_prompt, user_prompt)
    if not llm_text:
        llm_text = f"对 {len(gene_list)} 个基因进行了功能富集分析，结果见下方通路列表。"

    fig_md = "\n".join(f"![{os.path.basename(p)}](../{p})" for p in figures)
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    report = f"""# {project_name} — 功能富集分析报告

> **生成时间**: {now} | **物种**: {species} | **工具**: BioAgent v2.0

## 分析概况

| 指标 | 值 |
|------|------|
| 输入基因数 | {len(gene_list)} |
| 物种 | {species} |
| 富集数据库 | GO BP / KEGG / Reactome |

## 可视化

{fig_md}

## 分析解读

{llm_text}

---

"""
    # 通路表格
    for gs_name, df in enr_results.items():
        if df is None or df.empty:
            continue
        term_col = "Term" if "Term" in df.columns else df.columns[0]
        pval_col = "Adjusted P-value" if "Adjusted P-value" in df.columns else "P-value"
        if pval_col not in df.columns:
            pval_col = df.columns[1]
        report += f"## {gs_name} (Top 20)\n\n| 通路 | padj |\n|------|------|\n"
        for _, r in df.head(20).iterrows():
            report += f"| {r[term_col][:65]} | {r[pval_col]:.2e} |\n"
        report += "\n"

    report += f"\n---\n*BioAgent 自动生成 — {now}*\n"

    path = os.path.join(output_dir, "enrichment_report.md")
    with open(path, "w", encoding="utf-8") as f:
        f.write(report)
    return path


def run_enrichment_pipeline(gene_list_path: str, output_dir: str = "output",
                            project_name: str = "功能富集分析",
                            species: str = "human") -> dict:
    """运行独立富集分析 Pipeline。"""
    start = datetime.now()
    os.makedirs(output_dir, exist_ok=True)

    gene_list = load_gene_list(gene_list_path)
    enr_results = run_enrichment(gene_list, species)
    figures = plot_enrichment_results(enr_results, output_dir)
    report_path = generate_enrichment_report(gene_list, enr_results, figures,
                                              output_dir, project_name, species)

    elapsed = (datetime.now() - start).total_seconds()
    logger.info(f"富集分析 Pipeline 完成 | {elapsed:.1f}s")
    return {"report": report_path, "figures": figures, "elapsed": elapsed}
