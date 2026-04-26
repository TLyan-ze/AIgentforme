#!/usr/bin/env python3
"""
Task 1: RNA-seq 差异表达分析 + 功能富集 Pipeline
=================================================
从表达矩阵开始 → 数据质控 → DESeq2 → GO/KEGG 富集 → 可视化 → LLM 报告
支持物种: human / mouse
"""

import os
import json
import logging
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams["font.sans-serif"] = ["Arial Unicode MS", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

logger = logging.getLogger("BioAgent.RNAseq")


# ============================================================
#  1. 数据加载与质控
# ============================================================

def load_data(counts_path: str, metadata_path: str):
    """加载表达矩阵和样本信息。"""
    logger.info(f"加载表达矩阵: {counts_path}")
    raw_counts = pd.read_csv(counts_path, index_col=0)

    logger.info(f"加载样本信息: {metadata_path}")
    metadata = pd.read_csv(metadata_path, index_col=0)

    # 分离 gene_name 列（如果存在）
    gene_names = None
    if "gene_name" in raw_counts.columns:
        gene_names = raw_counts["gene_name"]
        counts = raw_counts.drop(columns=["gene_name"])
    else:
        counts = raw_counts

    # 确保数值类型
    counts = counts.astype(float)
    logger.info(f"数据维度: {counts.shape[0]} 基因 × {counts.shape[1]} 样本")
    return counts, metadata, gene_names


def quality_control(counts: pd.DataFrame, min_counts: int = 10,
                    min_samples: int = 2) -> dict:
    """数据质控：过滤低表达基因，计算样本统计。"""
    logger.info("执行数据质控...")

    total_before = len(counts)

    # 过滤低表达基因
    mask = (counts >= min_counts).sum(axis=1) >= min_samples
    filtered_idx = mask[mask].index
    genes_filtered = total_before - len(filtered_idx)

    qc_stats = {
        "total_genes_raw": total_before,
        "genes_after_filter": int(len(filtered_idx)),
        "genes_filtered": int(genes_filtered),
        "filter_criteria": f"counts>={min_counts} in >={min_samples} samples",
        "total_counts": int(counts.sum().sum()),
        "mean_counts_per_sample": counts.sum().mean().round(0).tolist(),
        "sample_names": counts.columns.tolist(),
    }

    logger.info(f"质控: {total_before} → {len(filtered_idx)} 基因 "
                f"(过滤 {genes_filtered} 低表达)")
    return qc_stats, filtered_idx


# ============================================================
#  2. 差异表达分析 (DESeq2 / scipy fallback)
# ============================================================

def _run_pydeseq2(counts: pd.DataFrame, metadata: pd.DataFrame,
                  condition_col: str = "condition") -> pd.DataFrame:
    """使用 pydeseq2 进行差异表达分析。"""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    # 确保 metadata 和 counts 的样本顺序一致
    common = counts.columns.intersection(metadata.index)
    counts = counts[common]
    metadata = metadata.loc[common]

    dds = DeseqDataSet(
        counts=counts.T,  # pydeseq2: samples × genes
        metadata=metadata,
        design_factors=condition_col,
    )
    dds.deseq2()

    groups = sorted(metadata[condition_col].unique())
    stat_res = DeseqStats(dds, contrast=[condition_col, groups[1], groups[0]])
    stat_res.summary()

    result = stat_res.results_df.copy()
    result.index.name = "gene_id"
    return result


def _run_simple_deg(counts: pd.DataFrame, metadata: pd.DataFrame,
                    condition_col: str = "condition") -> pd.DataFrame:
    """使用 scipy t-test 作为 DESeq2 不可用时的降级方案。"""
    from scipy import stats

    logger.warning("pydeseq2 未安装，使用 scipy t-test 进行简化差异分析")

    groups = sorted(metadata[condition_col].unique())
    common = counts.columns.intersection(metadata.index)
    counts = counts[common]
    metadata = metadata.loc[common]

    grp0 = metadata[metadata[condition_col] == groups[0]].index
    grp1 = metadata[metadata[condition_col] == groups[1]].index

    results = []
    for gene in counts.index:
        v0 = np.log2(counts.loc[gene, grp0].values + 1)
        v1 = np.log2(counts.loc[gene, grp1].values + 1)
        log2fc = float(np.mean(v1) - np.mean(v0))
        base_mean = float(counts.loc[gene].mean())

        if np.std(v0) == 0 and np.std(v1) == 0:
            pval = 1.0
        else:
            _, pval = stats.ttest_ind(v0, v1, equal_var=False)
            if np.isnan(pval):
                pval = 1.0

        results.append({
            "baseMean": base_mean,
            "log2FoldChange": log2fc,
            "lfcSE": 0.0,
            "stat": 0.0,
            "pvalue": pval,
        })

    df = pd.DataFrame(results, index=counts.index)

    # BH 多重检验校正
    pvals = df["pvalue"].values.copy()
    n = len(pvals)
    sorted_idx = np.argsort(pvals)
    adjusted = np.zeros(n)
    for i in range(n):
        adjusted[sorted_idx[i]] = pvals[sorted_idx[i]] * n / (i + 1)
    for i in range(n - 2, -1, -1):
        j = sorted_idx[i]
        jp1 = sorted_idx[i + 1]
        adjusted[j] = min(adjusted[j], adjusted[jp1])
    adjusted = np.clip(adjusted, 0, 1)
    df["padj"] = adjusted

    df.index.name = "gene_id"
    return df


def differential_expression(counts: pd.DataFrame, metadata: pd.DataFrame,
                            gene_names: pd.Series = None,
                            condition_col: str = "condition") -> pd.DataFrame:
    """运行差异表达分析。"""
    logger.info("运行差异表达分析...")

    try:
        result = _run_pydeseq2(counts, metadata, condition_col)
        method = "DESeq2 (pydeseq2)"
    except ImportError:
        result = _run_simple_deg(counts, metadata, condition_col)
        method = "t-test (scipy fallback)"
    except Exception as e:
        logger.warning(f"pydeseq2 运行失败 ({e})，切换到 t-test")
        result = _run_simple_deg(counts, metadata, condition_col)
        method = "t-test (scipy fallback)"

    # 添加基因名
    if gene_names is not None:
        result["gene_name"] = gene_names.reindex(result.index).fillna("")
    else:
        result["gene_name"] = ""

    # 计算显著差异统计
    n_sig = ((result["padj"] < 0.05) & (abs(result["log2FoldChange"]) > 1.0)).sum()
    n_up = ((result["padj"] < 0.05) & (result["log2FoldChange"] > 1.0)).sum()
    n_down = ((result["padj"] < 0.05) & (result["log2FoldChange"] < -1.0)).sum()
    logger.info(f"差异分析完成 [{method}]: 显著 {n_sig} (上调 {n_up}, 下调 {n_down})")

    result.attrs["method"] = method
    return result


# ============================================================
#  3. 功能富集分析 (GO / KEGG)
# ============================================================

def enrichment_analysis(deg_result: pd.DataFrame, species: str = "human",
                        fc_threshold: float = 1.0, pval_threshold: float = 0.05) -> dict:
    """对差异基因进行 GO 和 KEGG 富集分析。"""
    logger.info(f"执行功能富集分析 (物种: {species})...")

    # 获取显著差异基因列表
    sig = deg_result[
        (deg_result["padj"] < pval_threshold) &
        (abs(deg_result["log2FoldChange"]) > fc_threshold)
    ]
    gene_list = sig["gene_name"].replace("", np.nan).dropna().unique().tolist()

    if not gene_list:
        logger.warning("无有效基因名可用于富集分析，尝试使用 gene_id")
        gene_list = sig.index.tolist()

    if len(gene_list) < 5:
        logger.warning(f"差异基因过少 ({len(gene_list)})，跳过富集分析")
        return {"go": pd.DataFrame(), "kegg": pd.DataFrame(), "gene_count": len(gene_list)}

    try:
        import gseapy as gp

        # 选择基因集
        if species == "human":
            go_set = "GO_Biological_Process_2023"
            kegg_set = "KEGG_2021_Human"
        else:
            go_set = "GO_Biological_Process_2023"
            kegg_set = "KEGG_2021_Mouse"

        go_enr = gp.enrichr(gene_list=gene_list, gene_sets=[go_set],
                            organism=species, outdir=None, no_plot=True)
        kegg_enr = gp.enrichr(gene_list=gene_list, gene_sets=[kegg_set],
                               organism=species, outdir=None, no_plot=True)

        go_df = go_enr.results.head(30) if go_enr.results is not None else pd.DataFrame()
        kegg_df = kegg_enr.results.head(20) if kegg_enr.results is not None else pd.DataFrame()

        logger.info(f"富集分析完成: GO {len(go_df)} 条, KEGG {len(kegg_df)} 条")
        return {"go": go_df, "kegg": kegg_df, "gene_count": len(gene_list)}

    except ImportError:
        logger.warning("gseapy 未安装，跳过在线富集分析")
        return {"go": pd.DataFrame(), "kegg": pd.DataFrame(),
                "gene_count": len(gene_list), "fallback": True}
    except Exception as e:
        logger.warning(f"富集分析失败 ({e})，使用离线模式")
        return {"go": pd.DataFrame(), "kegg": pd.DataFrame(),
                "gene_count": len(gene_list), "error": str(e)}


# ============================================================
#  4. 可视化
# ============================================================

def plot_volcano(deg: pd.DataFrame, output_dir: str, fc_thr: float = 1.0,
                 pv_thr: float = 0.05) -> str:
    """火山图。"""
    fig, ax = plt.subplots(figsize=(10, 8))
    deg = deg.copy()
    deg["neg_log10_padj"] = -np.log10(deg["padj"].clip(lower=1e-50))

    ns = deg[(deg["padj"] >= pv_thr) | (abs(deg["log2FoldChange"]) <= fc_thr)]
    up = deg[(deg["padj"] < pv_thr) & (deg["log2FoldChange"] > fc_thr)]
    down = deg[(deg["padj"] < pv_thr) & (deg["log2FoldChange"] < -fc_thr)]

    ax.scatter(ns["log2FoldChange"], ns["neg_log10_padj"], c="grey", alpha=0.3, s=5)
    ax.scatter(up["log2FoldChange"], up["neg_log10_padj"], c="#E64B35", alpha=0.6, s=10, label=f"Up ({len(up)})")
    ax.scatter(down["log2FoldChange"], down["neg_log10_padj"], c="#4DBBD5", alpha=0.6, s=10, label=f"Down ({len(down)})")

    ax.axhline(-np.log10(pv_thr), color="k", ls="--", lw=0.8)
    ax.axvline(fc_thr, color="k", ls="--", lw=0.8)
    ax.axvline(-fc_thr, color="k", ls="--", lw=0.8)

    # 标注 top 基因
    top = pd.concat([up.nsmallest(5, "padj"), down.nsmallest(5, "padj")])
    for _, row in top.iterrows():
        name = row.get("gene_name", "") or row.name
        if isinstance(name, str) and name:
            ax.annotate(name, (row["log2FoldChange"], row["neg_log10_padj"]),
                        fontsize=7, alpha=0.8, xytext=(5, 5), textcoords="offset points")

    ax.set_xlabel("log2(Fold Change)", fontsize=12)
    ax.set_ylabel("-log10(adjusted p-value)", fontsize=12)
    ax.set_title("Volcano Plot", fontsize=14)
    ax.legend(fontsize=10)

    path = os.path.join(output_dir, "volcano_plot.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"火山图: {path}")
    return path


def plot_sample_correlation(counts: pd.DataFrame, output_dir: str) -> str:
    """样本相关性热图。"""
    log_counts = np.log2(counts + 1)
    corr = log_counts.corr()

    fig, ax = plt.subplots(figsize=(8, 7))
    sns.heatmap(corr, annot=True, fmt=".2f", cmap="RdBu_r", vmin=0.8, vmax=1,
                square=True, ax=ax, linewidths=0.5)
    ax.set_title("Sample Correlation Heatmap", fontsize=14)

    path = os.path.join(output_dir, "sample_correlation.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"相关性热图: {path}")
    return path


def plot_pca(counts: pd.DataFrame, metadata: pd.DataFrame,
             output_dir: str, condition_col: str = "condition") -> str:
    """PCA 降维可视化。"""
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    log_counts = np.log2(counts.T + 1)  # samples × genes
    scaler = StandardScaler()
    scaled = scaler.fit_transform(log_counts)

    pca = PCA(n_components=2)
    pcs = pca.fit_transform(scaled)

    fig, ax = plt.subplots(figsize=(8, 7))
    groups = metadata[condition_col].unique()
    colors = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488"]

    for i, grp in enumerate(groups):
        mask = metadata[condition_col] == grp
        ax.scatter(pcs[mask, 0], pcs[mask, 1],
                   c=colors[i % len(colors)], s=80, label=grp, edgecolors="white", lw=1)

    for j, sample in enumerate(metadata.index):
        ax.annotate(sample, (pcs[j, 0], pcs[j, 1]), fontsize=8,
                    xytext=(5, 5), textcoords="offset points")

    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
    ax.set_title("PCA Plot", fontsize=14)
    ax.legend()

    path = os.path.join(output_dir, "pca_plot.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"PCA 图: {path}")
    return path


def plot_top_genes_heatmap(counts: pd.DataFrame, deg: pd.DataFrame,
                           metadata: pd.DataFrame, output_dir: str,
                           n_top: int = 30, condition_col: str = "condition") -> str:
    """Top 差异基因热图。"""
    top_genes = deg.reindex(deg["padj"].abs().nsmallest(n_top).index)
    top_ids = top_genes.index.tolist()

    log_c = np.log2(counts.loc[top_ids] + 1)
    # Z-score normalization
    z = log_c.sub(log_c.mean(axis=1), axis=0).div(log_c.std(axis=1), axis=0)
    z = z.clip(-3, 3)

    # 按条件排序
    sorted_meta = metadata.sort_values(condition_col)
    z = z[sorted_meta.index]

    # 行标签用基因名
    labels = []
    for gid in top_ids:
        name = top_genes.loc[gid, "gene_name"] if "gene_name" in top_genes.columns else ""
        labels.append(name if name else gid)

    fig, ax = plt.subplots(figsize=(10, 10))
    sns.heatmap(z, cmap="RdBu_r", center=0, yticklabels=labels,
                xticklabels=sorted_meta.index, ax=ax, linewidths=0.3)
    ax.set_title(f"Top {n_top} DEGs Heatmap", fontsize=14)
    ax.set_ylabel("")

    path = os.path.join(output_dir, "top_genes_heatmap.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"热图: {path}")
    return path


def plot_enrichment(enr_results: dict, output_dir: str) -> list:
    """绘制 GO 和 KEGG 富集气泡图。"""
    paths = []
    for key, label in [("go", "GO Biological Process"), ("kegg", "KEGG Pathway")]:
        df = enr_results.get(key)
        if df is None or df.empty:
            continue

        # 兼容 gseapy 返回的列名
        term_col = "Term" if "Term" in df.columns else df.columns[0]
        pval_col = "Adjusted P-value" if "Adjusted P-value" in df.columns else "P-value"
        if pval_col not in df.columns:
            pval_col = df.columns[1]

        top = df.head(15).copy()
        top["-log10(padj)"] = -np.log10(top[pval_col].clip(lower=1e-50))

        fig, ax = plt.subplots(figsize=(10, 8))
        ax.barh(range(len(top)), top["-log10(padj)"], color="#3C5488", alpha=0.8)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top[term_col].str[:60], fontsize=8)
        ax.set_xlabel("-log10(adjusted p-value)")
        ax.set_title(f"{label} Enrichment (Top 15)", fontsize=14)
        ax.invert_yaxis()

        path = os.path.join(output_dir, f"enrichment_{key}.png")
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        paths.append(path)
        logger.info(f"富集图 [{label}]: {path}")

    return paths


# ============================================================
#  5. LLM 报告生成
# ============================================================

def _call_llm(system_prompt: str, user_prompt: str) -> str:
    """调用 LLM API。"""
    try:
        from openai import OpenAI
    except ImportError:
        return None

    import os
    api_key = os.environ.get("LLM_API_KEY", "")
    base_url = os.environ.get("LLM_BASE_URL", "https://api.openai.com/v1")
    model = os.environ.get("LLM_MODEL", "gpt-4o-mini")

    if not api_key or api_key == "your_api_key_here":
        return None

    client = OpenAI(api_key=api_key, base_url=base_url)
    resp = client.chat.completions.create(
        model=model,
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ],
        temperature=0.3, max_tokens=4000,
    )
    return resp.choices[0].message.content


def _build_stats_summary(deg: pd.DataFrame, fc_thr: float, pv_thr: float,
                         species: str, method: str) -> dict:
    """构建统计摘要字典。"""
    sig = deg[(deg["padj"] < pv_thr) & (abs(deg["log2FoldChange"]) > fc_thr)]
    up = sig[sig["log2FoldChange"] > 0]
    down = sig[sig["log2FoldChange"] < 0]

    def _gene_list(sub, n=10):
        return [
            {"gene": row.get("gene_name", "") or idx, "log2FC": round(float(row["log2FoldChange"]), 4),
             "padj": f"{row['padj']:.2e}"}
            for idx, row in sub.head(n).iterrows()
        ]

    return {
        "species": species, "method": method,
        "total_genes": len(deg),
        "sig_genes": len(sig), "sig_ratio": round(len(sig) / len(deg) * 100, 2),
        "up_genes": len(up), "down_genes": len(down),
        "fc_threshold": fc_thr, "pval_threshold": pv_thr,
        "log2fc_range": (round(float(deg["log2FoldChange"].min()), 4),
                         round(float(deg["log2FoldChange"].max()), 4)),
        "top_up": _gene_list(up.sort_values("padj")),
        "top_down": _gene_list(down.sort_values("padj")),
    }


def generate_report(stats: dict, enr_results: dict, figures: list,
                    qc_stats: dict, output_dir: str,
                    project_name: str, species: str) -> str:
    """生成完整 Markdown 分析报告。"""

    # 构建 LLM 提示
    sys_prompt = (
        "你是一位资深生物信息学分析师。根据差异表达和功能富集分析结果，"
        "用中文撰写一份专业、结构化的分析报告。风格：学术正式，条理清晰。"
    )

    top_up_str = "\n".join(f"  {g['gene']}: log2FC={g['log2FC']}" for g in stats["top_up"][:8])
    top_down_str = "\n".join(f"  {g['gene']}: log2FC={g['log2FC']}" for g in stats["top_down"][:8])

    enr_text = ""
    for key, label in [("go", "GO"), ("kegg", "KEGG")]:
        df = enr_results.get(key)
        if df is not None and not df.empty:
            term_col = "Term" if "Term" in df.columns else df.columns[0]
            pval_col = "Adjusted P-value" if "Adjusted P-value" in df.columns else "P-value"
            if pval_col not in df.columns:
                pval_col = df.columns[1]
            enr_text += f"\n### {label} Top 通路:\n"
            for _, r in df.head(10).iterrows():
                enr_text += f"- {r[term_col][:80]} (padj={r[pval_col]:.2e})\n"

    user_prompt = f"""请根据以下分析结果撰写差异表达分析报告：

## 项目: {project_name}
## 物种: {species}
## 分析方法: {stats['method']}

## 统计结果
- 总基因数: {stats['total_genes']}
- 显著差异基因: {stats['sig_genes']} ({stats['sig_ratio']}%)
- 上调: {stats['up_genes']}, 下调: {stats['down_genes']}
- log2FC 范围: [{stats['log2fc_range'][0]}, {stats['log2fc_range'][1]}]

## Top 上调基因
{top_up_str}

## Top 下调基因
{top_down_str}

{enr_text}

请撰写包含以下部分的报告:
1. 分析概述（背景、方法）
2. 差异表达结果解读
3. 关键基因功能讨论
4. 富集通路分析（如有）
5. 生物学意义与后续建议
"""

    llm_text = _call_llm(sys_prompt, user_prompt)
    if not llm_text:
        llm_text = (
            "## 分析报告\n\n"
            "> [离线模式] 配置 .env 中的 LLM_API_KEY 可获得 AI 增强报告。\n\n"
            f"### 概述\n本次使用 {stats['method']} 对 {species} RNA-seq 数据进行差异表达分析，"
            f"共鉴定 {stats['sig_genes']} 个显著差异基因（上调 {stats['up_genes']}，"
            f"下调 {stats['down_genes']}）。\n\n"
            "### 建议\n- 对显著差异基因进行 GO/KEGG 富集分析\n"
            "- 结合临床表型进一步验证\n- 使用 qRT-PCR 验证关键基因"
        )

    # 构建富集表格
    go_table = ""
    go_df = enr_results.get("go")
    if go_df is not None and not go_df.empty:
        term_col = "Term" if "Term" in go_df.columns else go_df.columns[0]
        pval_col = "Adjusted P-value" if "Adjusted P-value" in go_df.columns else "P-value"
        if pval_col not in go_df.columns:
            pval_col = go_df.columns[1]
        go_table = "### GO 富集结果 (Top 15)\n\n| 通路 | padj |\n|------|------|\n"
        for _, r in go_df.head(15).iterrows():
            go_table += f"| {r[term_col][:60]} | {r[pval_col]:.2e} |\n"

    kegg_table = ""
    kegg_df = enr_results.get("kegg")
    if kegg_df is not None and not kegg_df.empty:
        term_col = "Term" if "Term" in kegg_df.columns else kegg_df.columns[0]
        pval_col = "Adjusted P-value" if "Adjusted P-value" in kegg_df.columns else "P-value"
        if pval_col not in kegg_df.columns:
            pval_col = kegg_df.columns[1]
        kegg_table = "### KEGG 富集结果 (Top 15)\n\n| 通路 | padj |\n|------|------|\n"
        for _, r in kegg_df.head(15).iterrows():
            kegg_table += f"| {r[term_col][:60]} | {r[pval_col]:.2e} |\n"

    fig_md = "\n".join(f"![{Path(p).stem}](../{p})" for p in figures)

    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    report = f"""# {project_name} — RNA-seq 差异表达分析报告

> **生成时间**: {now}
> **物种**: {species}
> **分析方法**: {stats['method']}
> **生成工具**: BioAgent v2.0

---

## 数据质控

| 指标 | 值 |
|------|------|
| 原始基因数 | {qc_stats['total_genes_raw']} |
| 质控后基因数 | {qc_stats['genes_after_filter']} |
| 过滤基因数 | {qc_stats['genes_filtered']} |
| 过滤标准 | {qc_stats['filter_criteria']} |
| 总 reads 数 | {qc_stats['total_counts']:,} |

## 差异表达统计

| 指标 | 值 |
|------|------|
| 分析基因数 | {stats['total_genes']} |
| 显著差异基因 | {stats['sig_genes']} ({stats['sig_ratio']}%) |
| 上调基因 | {stats['up_genes']} |
| 下调基因 | {stats['down_genes']} |
| log2FC 范围 | [{stats['log2fc_range'][0]}, {stats['log2fc_range'][1]}] |
| 筛选阈值 | \\|log2FC\\| > {stats['fc_threshold']}, padj < {stats['pval_threshold']} |

## 可视化

{fig_md}

---

## LLM 分析报告

{llm_text}

---

{go_table}

{kegg_table}

---

## Top 差异基因

### 上调 (Top 15)

| 基因 | log2FC | padj |
|------|--------|------|
"""
    for g in stats["top_up"][:15]:
        report += f"| {g['gene']} | {g['log2FC']} | {g['padj']} |\n"

    report += "\n### 下调 (Top 15)\n\n| 基因 | log2FC | padj |\n|------|--------|------|\n"
    for g in stats["top_down"][:15]:
        report += f"| {g['gene']} | {g['log2FC']} | {g['padj']} |\n"

    report += f"\n---\n*BioAgent 自动生成 — {now}*\n"

    report_path = os.path.join(output_dir, "analysis_report.md")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(report)
    logger.info(f"报告: {report_path}")
    return report_path


# ============================================================
#  6. Pipeline 主入口
# ============================================================

def run_pipeline(counts_path: str, metadata_path: str,
                 output_dir: str = "output", project_name: str = "RNA-seq DEG",
                 species: str = "human", fc_threshold: float = 1.0,
                 pval_threshold: float = 0.05,
                 condition_col: str = "condition") -> dict:
    """运行完整的 RNA-seq 差异表达 + 富集分析 Pipeline。"""
    start = datetime.now()
    logger.info("=" * 60)
    logger.info(f"BioAgent RNA-seq Pipeline 启动")
    logger.info(f"物种: {species} | 输出: {output_dir}")
    logger.info("=" * 60)

    os.makedirs(output_dir, exist_ok=True)

    # 1. 加载数据
    counts, metadata, gene_names = load_data(counts_path, metadata_path)

    # 2. 质控
    qc_stats, keep_idx = quality_control(counts)
    counts_filtered = counts.loc[keep_idx]
    gene_names_filtered = gene_names.reindex(keep_idx) if gene_names is not None else None

    # 3. 差异表达
    deg_result = differential_expression(counts_filtered, metadata, gene_names_filtered, condition_col)
    method = deg_result.attrs.get("method", "unknown")

    # 4. 富集分析
    enr_results = enrichment_analysis(deg_result, species, fc_threshold, pval_threshold)

    # 5. 统计摘要
    stats = _build_stats_summary(deg_result, fc_threshold, pval_threshold, species, method)

    # 6. 可视化
    figures = []
    for fn, args in [
        (plot_volcano, (deg_result, output_dir, fc_threshold, pval_threshold)),
        (plot_sample_correlation, (counts_filtered, output_dir)),
        (plot_pca, (counts_filtered, metadata, output_dir, condition_col)),
        (plot_top_genes_heatmap, (counts_filtered, deg_result, metadata, output_dir)),
    ]:
        try:
            figures.append(fn(*args))
        except Exception as e:
            logger.error(f"{fn.__name__} 失败: {e}")

    try:
        figures.extend(plot_enrichment(enr_results, output_dir))
    except Exception as e:
        logger.error(f"富集图失败: {e}")

    # 7. 保存中间结果
    deg_result.to_csv(os.path.join(output_dir, "deg_results.csv"))
    with open(os.path.join(output_dir, "statistics.json"), "w", encoding="utf-8") as f:
        json.dump(stats, f, ensure_ascii=False, indent=2, default=str)

    # 8. 生成报告
    report_path = generate_report(stats, enr_results, figures, qc_stats,
                                  output_dir, project_name, species)

    elapsed = (datetime.now() - start).total_seconds()
    logger.info("=" * 60)
    logger.info(f"Pipeline 完成 | 耗时 {elapsed:.1f}s")
    logger.info(f"报告: {report_path}")
    logger.info("=" * 60)

    return {
        "report": report_path, "figures": figures, "stats": stats,
        "qc": qc_stats, "elapsed": elapsed,
    }
