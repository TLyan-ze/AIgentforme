#!/usr/bin/env python3
"""
BioReportAgent — 生物信息学差异表达分析报告智能体
===================================================

功能：读取差异表达分析结果 (CSV)，自动完成统计分析、可视化、
     并调用 LLM 生成结构化中文分析报告。

作者：颜泽钦
单位：某三甲医院 实验研究部
"""

import os
import sys
import json
import argparse
import logging
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# ── 日志配置 ──────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("BioReportAgent")

# ── matplotlib 中文字体配置 ───────────────────────────────
plt.rcParams["font.sans-serif"] = ["Arial Unicode MS", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


# ============================================================
#  1. 数据加载模块
# ============================================================

def load_deg_results(filepath: str) -> pd.DataFrame:
    """加载差异表达分析结果 CSV 文件。

    期望列：gene_name 或 gene_id, log2FoldChange, pvalue, padj
    """
    logger.info(f"加载数据: {filepath}")
    df = pd.read_csv(filepath)
    required_cols = {"log2FoldChange", "pvalue", "padj"}
    if not required_cols.issubset(set(df.columns)):
        raise ValueError(f"数据缺少必要列，需要: {required_cols}，实际有: {set(df.columns)}")
    logger.info(f"数据加载完成，共 {len(df)} 个基因")
    return df


# ============================================================
#  2. 统计分析模块
# ============================================================

def run_statistical_analysis(df: pd.DataFrame,
                             fc_threshold: float = 1.0,
                             pval_threshold: float = 0.05) -> dict:
    """对差异表达结果进行统计分析。"""
    logger.info("执行统计分析...")

    total_genes = len(df)

    # 显著差异基因
    sig_mask = (df["padj"] < pval_threshold) & (abs(df["log2FoldChange"]) > fc_threshold)
    sig_genes = df[sig_mask]
    n_sig = len(sig_genes)

    # 上调 / 下调
    up_mask = (df["padj"] < pval_threshold) & (df["log2FoldChange"] > fc_threshold)
    down_mask = (df["padj"] < pval_threshold) & (df["log2FoldChange"] < -fc_threshold)
    n_up = up_mask.sum()
    n_down = down_mask.sum()

    up_genes = df[up_mask].sort_values("padj")
    down_genes = df[down_mask].sort_values("padj")

    # Top 差异基因
    top_up = up_genes.head(20)
    top_down = down_genes.head(20)

    # 统计摘要
    stats = {
        "total_genes": int(total_genes),
        "significant_genes": int(n_sig),
        "significant_ratio": round(n_sig / total_genes * 100, 2),
        "upregulated": int(n_up),
        "downregulated": int(n_down),
        "fc_threshold": fc_threshold,
        "pval_threshold": pval_threshold,
        "log2fc_range": (
            round(float(df["log2FoldChange"].min()), 4),
            round(float(df["log2FoldChange"].max()), 4),
        ),
        "log2fc_mean": round(float(df["log2FoldChange"].mean()), 4),
        "log2fc_median": round(float(df["log2FoldChange"].median()), 4),
        "top_up_genes": _extract_gene_list(top_up),
        "top_down_genes": _extract_gene_list(top_down),
    }

    logger.info(f"统计完成: 总基因 {total_genes}, 显著差异 {n_sig} "
                f"(上调 {n_up}, 下调 {n_down})")
    return stats


def _extract_gene_list(subdf: pd.DataFrame) -> list:
    """从子数据框提取基因信息列表。"""
    results = []
    name_col = "gene_name" if "gene_name" in subdf.columns else "gene_id"
    for _, row in subdf.iterrows():
        name = row.get("gene_name", "")
        if pd.isna(name) or name == "":
            name = row.get("gene_id", f"unknown_{_}")
        results.append({
            "gene": name,
            "log2FC": round(float(row["log2FoldChange"]), 4),
            "padj": f"{row['padj']:.2e}",
        })
    return results


# ============================================================
#  3. 可视化模块
# ============================================================

def plot_volcano(df: pd.DataFrame, stats: dict, output_dir: str) -> str:
    """绘制火山图。"""
    logger.info("绘制火山图...")
    fc_thr = stats["fc_threshold"]
    pv_thr = stats["pval_threshold"]

    fig, ax = plt.subplots(figsize=(10, 8))

    # 分类着色
    df["neg_log10_padj"] = -np.log10(df["padj"].clip(lower=1e-50))
    ns = df[(df["padj"] >= pv_thr) | (abs(df["log2FoldChange"]) <= fc_thr)]
    up = df[(df["padj"] < pv_thr) & (df["log2FoldChange"] > fc_thr)]
    down = df[(df["padj"] < pv_thr) & (df["log2FoldChange"] < -fc_thr)]

    ax.scatter(ns["log2FoldChange"], ns["neg_log10_padj"],
               c="grey", alpha=0.3, s=5, label="Not significant")
    ax.scatter(up["log2FoldChange"], up["neg_log10_padj"],
               c="#E64B35", alpha=0.6, s=10, label=f"Up ({len(up)})")
    ax.scatter(down["log2FoldChange"], down["neg_log10_padj"],
               c="#4DBBD5", alpha=0.6, s=10, label=f"Down ({len(down)})")

    # 阈值线
    ax.axhline(-np.log10(pv_thr), color="black", linestyle="--", linewidth=0.8)
    ax.axvline(fc_thr, color="black", linestyle="--", linewidth=0.8)
    ax.axvline(-fc_thr, color="black", linestyle="--", linewidth=0.8)

    # 标注 top 基因
    name_col = "gene_name" if "gene_name" in df.columns else "gene_id"
    top_genes = pd.concat([
        up.nsmallest(5, "padj"),
        down_genes := down.nsmallest(5, "padj"),
    ])
    for _, row in top_genes.iterrows():
        name = row.get("gene_name", "")
        if pd.isna(name) or name == "":
            name = row.get("gene_id", "")
        if name:
            ax.annotate(name, (row["log2FoldChange"], row["neg_log10_padj"]),
                        fontsize=7, alpha=0.8,
                        xytext=(5, 5), textcoords="offset points")

    ax.set_xlabel("log2(Fold Change)", fontsize=12)
    ax.set_ylabel("-log10(adjusted p-value)", fontsize=12)
    ax.set_title("Volcano Plot — Differential Expression Analysis", fontsize=14)
    ax.legend(fontsize=10)

    path = os.path.join(output_dir, "volcano_plot.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"火山图已保存: {path}")
    return path


def plot_ma(df: pd.DataFrame, stats: dict, output_dir: str) -> str:
    """绘制 MA 图。"""
    logger.info("绘制 MA 图...")
    fc_thr = stats["fc_threshold"]
    pv_thr = stats["pval_threshold"]
    base_col = "baseMean" if "baseMean" in df.columns else None

    fig, ax = plt.subplots(figsize=(10, 8))

    if base_col:
        x = np.log2(df[base_col].clip(lower=0.1))
        ax.set_xlabel("log2(baseMean)", fontsize=12)
    else:
        x = np.arange(len(df))
        ax.set_xlabel("Gene index", fontsize=12)

    y = df["log2FoldChange"]
    sig = df["padj"] < pv_thr
    ns = ~sig

    ax.scatter(x[ns], y[ns], c="grey", alpha=0.2, s=3)
    up = sig & (y > fc_thr)
    down = sig & (y < -fc_thr)
    ax.scatter(x[up], y[up], c="#E64B35", alpha=0.5, s=8, label=f"Up ({up.sum()})")
    ax.scatter(x[down], y[down], c="#4DBBD5", alpha=0.5, s=8, label=f"Down ({down.sum()})")

    ax.axhline(0, color="black", linewidth=0.5)
    ax.axhline(fc_thr, color="grey", linestyle="--", linewidth=0.5)
    ax.axhline(-fc_thr, color="grey", linestyle="--", linewidth=0.5)
    ax.set_ylabel("log2(Fold Change)", fontsize=12)
    ax.set_title("MA Plot — Differential Expression Analysis", fontsize=14)
    ax.legend(fontsize=10)

    path = os.path.join(output_dir, "ma_plot.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"MA 图已保存: {path}")
    return path


def plot_distribution(df: pd.DataFrame, output_dir: str) -> str:
    """绘制 log2FC 分布直方图。"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # log2FC 分布
    axes[0].hist(df["log2FoldChange"], bins=80, color="#3C5488", alpha=0.7, edgecolor="white")
    axes[0].axvline(0, color="red", linestyle="--")
    axes[0].set_xlabel("log2(Fold Change)")
    axes[0].set_ylabel("Frequency")
    axes[0].set_title("Distribution of log2(Fold Change)")

    # p 值分布
    axes[1].hist(df["pvalue"], bins=50, color="#00A087", alpha=0.7, edgecolor="white")
    axes[1].set_xlabel("p-value")
    axes[1].set_ylabel("Frequency")
    axes[1].set_title("Distribution of p-values")

    path = os.path.join(output_dir, "distribution_plots.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"分布图已保存: {path}")
    return path


# ============================================================
#  4. LLM 报告生成模块
# ============================================================

def call_llm(system_prompt: str, user_prompt: str) -> str:
    """调用 LLM API 生成文本。支持 OpenAI 兼容接口。"""
    try:
        from openai import OpenAI
    except ImportError:
        logger.warning("openai 库未安装，使用离线模板生成报告")
        return _offline_report_template(system_prompt, user_prompt)

    api_key = os.environ.get("LLM_API_KEY", "")
    base_url = os.environ.get("LLM_BASE_URL", "https://api.openai.com/v1")
    model = os.environ.get("LLM_MODEL", "gpt-4o-mini")

    if not api_key or api_key == "your_api_key_here":
        logger.warning("未配置 LLM_API_KEY，使用离线模板生成报告")
        return _offline_report_template(system_prompt, user_prompt)

    logger.info(f"调用 LLM: {model} @ {base_url}")
    client = OpenAI(api_key=api_key, base_url=base_url)

    response = client.chat.completions.create(
        model=model,
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ],
        temperature=0.3,
        max_tokens=4000,
    )
    return response.choices[0].message.content


def _offline_report_template(system_prompt: str, user_prompt: str) -> str:
    """离线模式：当无法连接 LLM 时，使用模板生成报告框架。"""
    return """## 差异表达分析报告（离线模板模式）

> [离线模式] 未配置 LLM API，以下为模板生成内容。
> 配置 .env 文件中的 LLM_API_KEY 后可获得 AI 增强的分析报告。

### 1. 分析概述

本次分析对两组样本进行了差异表达基因（DEG）分析，采用 DESeq2 方法进行统计检验。

### 2. 结果摘要

根据设定的阈值（|log2FC| > 1.0, padj < 0.05），筛选出显著差异表达的基因。
其中上调基因和下调基因的分布见火山图和 MA 图。

### 3. 关键发现

- 上调基因中包含多个与肿瘤发生发展相关的基因
- 下调基因中包含肿瘤抑制基因
- 部分差异基因可能作为潜在的生物标志物

### 4. 建议

- 对显著差异基因进行 GO/KEGG 富集分析
- 结合临床表型数据进一步验证
- 考虑使用 qRT-PCR 对关键基因进行实验验证

---
*本报告由 BioReportAgent 离线模板生成*
"""


def generate_report(stats: dict, figure_paths: list, output_dir: str,
                    project_name: str = "差异表达分析") -> str:
    """生成完整的分析报告。"""

    # 构建 LLM 提示
    system_prompt = (
        "你是一位资深的生物信息学分析师，擅长转录组学和差异表达分析。"
        "请根据提供的统计结果，用中文撰写一份专业、结构化的差异表达分析报告。"
        "报告应包含：分析概述、结果解读、关键基因讨论、生物学意义、后续建议。"
        "语言风格：学术正式，条理清晰。"
    )

    top_up_str = "\n".join(
        [f"  {g['gene']}: log2FC={g['log2FC']}, padj={g['padj']}" for g in stats["top_up_genes"][:10]]
    )
    top_down_str = "\n".join(
        [f"  {g['gene']}: log2FC={g['log2FC']}, padj={g['padj']}" for g in stats["top_down_genes"][:10]]
    )

    user_prompt = f"""请根据以下差异表达分析统计数据撰写分析报告：

## 项目名称
{project_name}

## 分析参数
- Fold Change 阈值: {stats['fc_threshold']}
- 校正 p 值阈值: {stats['pval_threshold']}

## 统计结果
- 总基因数: {stats['total_genes']}
- 显著差异基因数: {stats['significant_genes']} ({stats['significant_ratio']}%)
- 上调基因: {stats['upregulated']}
- 下调基因: {stats['downregulated']}
- log2FC 范围: [{stats['log2fc_range'][0]}, {stats['log2fc_range'][1]}]
- log2FC 均值: {stats['log2fc_mean']}
- log2FC 中位数: {stats['log2fc_median']}

## Top 上调基因 (Top 10)
{top_up_str}

## Top 下调基因 (Top 10)
{top_down_str}

请撰写一份完整的中文分析报告，包含以下部分：
1. 分析概述（实验背景和方法简介）
2. 结果解读（统计数据的意义说明）
3. 关键基因讨论（Top 上调和下调基因的功能和意义）
4. 生物学意义（可能的通路和调控机制）
5. 后续分析建议

注意：如果基因名看起来是 ENSG 编号而非标准基因名，请基于基因编号进行分析。
"""

    # 调用 LLM 生成报告正文
    llm_content = call_llm(system_prompt, user_prompt)

    # 组装完整报告
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    figures_md = "\n".join([f"![{Path(p).stem}](../{p})" for p in figure_paths])

    report = f"""# {project_name} — 差异表达分析报告

> **生成时间**: {now}
> **生成工具**: BioReportAgent v1.0
> **分析平台**: 实验研究部 生物信息分析平台

---

## 统计概览

| 指标 | 数值 |
|------|------|
| 总基因数 | {stats['total_genes']} |
| 显著差异基因数 | {stats['significant_genes']} ({stats['significant_ratio']}%) |
| 上调基因数 | {stats['upregulated']} |
| 下调基因数 | {stats['downregulated']} |
| log2FC 范围 | [{stats['log2fc_range'][0]}, {stats['log2fc_range'][1]}] |
| 阈值 | \\|log2FC\\| > {stats['fc_threshold']}, padj < {stats['pval_threshold']} |

---

## 可视化结果

{figures_md}

---

## LLM 分析报告

{llm_content}

---

## Top 差异基因列表

### 上调基因 (Top 20)

| 基因 | log2FC | padj |
|------|--------|------|
"""
    for g in stats["top_up_genes"][:20]:
        report += f"| {g['gene']} | {g['log2FC']} | {g['padj']} |\n"

    report += f"""
### 下调基因 (Top 20)

| 基因 | log2FC | padj |
|------|--------|------|
"""
    for g in stats["top_down_genes"][:20]:
        report += f"| {g['gene']} | {g['log2FC']} | {g['padj']} |\n"

    report += f"""
---

*本报告由 BioReportAgent 自动生成，仅供科研参考，不构成临床诊断依据。*
*实验研究部 — 颜泽钦*
"""

    # 保存报告
    report_path = os.path.join(output_dir, "analysis_report.md")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(report)
    logger.info(f"分析报告已保存: {report_path}")
    return report_path


# ============================================================
#  5. 主流程 (Agent Pipeline)
# ============================================================

def run_agent(input_file: str, output_dir: str = "output",
              project_name: str = "差异表达分析",
              fc_threshold: float = 1.0,
              pval_threshold: float = 0.05,
              use_llm: bool = True) -> dict:
    """BioReportAgent 主流程。

    Parameters
    ----------
    input_file : str
        差异表达结果 CSV 文件路径
    output_dir : str
        输出目录
    project_name : str
        项目名称
    fc_threshold : float
        log2FC 阈值
    pval_threshold : float
        校正 p 值阈值
    use_llm : bool
        是否使用 LLM 生成报告

    Returns
    -------
    dict
        包含所有输出路径和统计结果的字典
    """
    start_time = datetime.now()
    logger.info("=" * 60)
    logger.info("BioReportAgent 启动")
    logger.info(f"输入文件: {input_file}")
    logger.info(f"输出目录: {output_dir}")
    logger.info("=" * 60)

    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)

    # Step 1: 数据加载
    df = load_deg_results(input_file)

    # Step 2: 统计分析
    stats = run_statistical_analysis(df, fc_threshold, pval_threshold)

    # Step 3: 可视化
    figure_paths = []
    try:
        figure_paths.append(plot_volcano(df, stats, output_dir))
    except Exception as e:
        logger.error(f"火山图生成失败: {e}")
    try:
        figure_paths.append(plot_ma(df, stats, output_dir))
    except Exception as e:
        logger.error(f"MA 图生成失败: {e}")
    try:
        figure_paths.append(plot_distribution(df, output_dir))
    except Exception as e:
        logger.error(f"分布图生成失败: {e}")

    # Step 4: LLM 报告生成
    if use_llm:
        report_path = generate_report(stats, figure_paths, output_dir, project_name)
    else:
        report_path = generate_report(stats, figure_paths, output_dir, project_name)
        logger.info("LLM 报告生成完成")

    # Step 5: 保存统计数据 (JSON)
    stats_path = os.path.join(output_dir, "statistics.json")
    with open(stats_path, "w", encoding="utf-8") as f:
        json.dump(stats, f, ensure_ascii=False, indent=2)

    elapsed = (datetime.now() - start_time).total_seconds()
    logger.info("=" * 60)
    logger.info(f"BioReportAgent 完成，耗时 {elapsed:.2f} 秒")
    logger.info(f"报告: {report_path}")
    logger.info(f"统计: {stats_path}")
    logger.info("=" * 60)

    return {
        "report_path": report_path,
        "stats_path": stats_path,
        "figure_paths": figure_paths,
        "statistics": stats,
        "elapsed_seconds": elapsed,
    }


# ============================================================
#  6. CLI 入口
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="BioReportAgent — 生物信息学差异表达分析报告智能体",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 使用示例数据运行
  python bio_report_agent.py --sample-data

  # 指定输入文件运行
  python bio_report_agent.py -i data/deg_results.csv -o output/

  # 自定义阈值
  python bio_report_agent.py -i data/deg_results.csv --fc 1.5 --pval 0.01
        """,
    )
    parser.add_argument("-i", "--input", type=str, help="差异表达结果 CSV 文件路径")
    parser.add_argument("-o", "--output", type=str, default="output", help="输出目录 (默认: output)")
    parser.add_argument("-n", "--name", type=str, default="差异表达分析", help="项目名称")
    parser.add_argument("--fc", type=float, default=1.0, help="log2FC 阈值 (默认: 1.0)")
    parser.add_argument("--pval", type=float, default=0.05, help="校正 p 值阈值 (默认: 0.05)")
    parser.add_argument("--sample-data", action="store_true", help="使用内置示例数据运行演示")
    parser.add_argument("--no-llm", action="store_true", help="不调用 LLM，使用离线模板")

    args = parser.parse_args()

    # 加载 .env
    env_file = Path(__file__).parent.parent / ".env"
    if env_file.exists():
        with open(env_file) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#") and "=" in line:
                    key, val = line.split("=", 1)
                    os.environ.setdefault(key.strip(), val.strip())

    if args.sample_data:
        # 生成示例数据
        sys.path.insert(0, str(Path(__file__).parent))
        from generate_sample_data import generate_deg_results
        sample_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
        os.makedirs(sample_dir, exist_ok=True)
        sample_path = os.path.join(sample_dir, "sample_deg_results.csv")
        df = generate_deg_results()
        df.to_csv(sample_path, index=False)
        logger.info(f"已生成示例数据: {sample_path}")
        input_file = sample_path
    elif args.input:
        input_file = args.input
    else:
        parser.print_help()
        print("\n[错误] 请指定 --input 或 --sample-data")
        sys.exit(1)

    result = run_agent(
        input_file=input_file,
        output_dir=args.output,
        project_name=args.name,
        fc_threshold=args.fc,
        pval_threshold=args.pval,
        use_llm=not args.no_llm,
    )

    # 打印结果摘要
    print(f"\n{'='*60}")
    print(f"分析完成!")
    print(f"  报告文件: {result['report_path']}")
    print(f"  统计数据: {result['stats_path']}")
    print(f"  可视化图表: {len(result['figure_paths'])} 张")
    print(f"  显著差异基因: {result['statistics']['significant_genes']} "
          f"(上调 {result['statistics']['upregulated']}, "
          f"下调 {result['statistics']['downregulated']})")
    print(f"  耗时: {result['elapsed_seconds']:.2f} 秒")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
