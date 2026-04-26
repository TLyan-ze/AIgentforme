#!/usr/bin/env python3
"""生成模拟差异表达分析结果数据，用于演示 BioReportAgent。"""

import numpy as np
import pandas as pd


def generate_deg_results(n_genes: int = 5000, n_sig_up: int = 320, n_sig_down: int = 280,
                         random_seed: int = 42) -> pd.DataFrame:
    """生成模拟的 RNA-seq 差异表达分析结果。

    Parameters
    ----------
    n_genes : int
        总基因数
    n_sig_up : int
        显著上调基因数
    n_sig_down : int
        显著下调基因数
    random_seed : int
        随机种子

    Returns
    -------
    pd.DataFrame
        包含 log2FoldChange, pvalue, padj 列的差异表达结果
    """
    np.random.seed(random_seed)

    # 生成基因名 (ENSG + 数字)
    gene_ids = [f"ENSG{i:09d}" for i in range(1, n_genes + 1)]

    # 从真实人类基因中选一些有意义的名字用于展示
    real_genes = [
        "TP53", "BRCA1", "MYC", "EGFR", "KRAS", "PTEN", "CDH1", "RB1",
        "PIK3CA", "APC", "VEGFA", "FGFR1", "ERBB2", "NOTCH1", "WNT1",
        "BCL2", "MDM2", "CDK4", "CCND1", "JAK2", "STAT3", "MTOR",
        "AKT1", "MAPK1", "RAF1", "HRAS", "NRAS", "SOS1", "GRB2",
        "SHC1", "IRS1", "IGF1R", "PDGFRA", "KIT", "FLT3", "RET",
        "MET", "ALK", "ROS1", "NTRK1", "FGFR2", "FGFR3", "DDR2",
        "EPHA2", "EPHB1", "AXL", "MERTK", "TIE2", "VEGFR2", "PDGFRB",
        "ESR1", "PGR", "AR", "FOXA1", "GATA3", "CDK6", "CCNE1",
        "CDKN1A", "CDKN2A", "CDKN2B", "TP73", "TP63", "ATM", "ATR",
        "CHEK1", "CHEK2", "BRCA2", "PALB2", "RAD51", "XRCC1",
    ]

    # 基础 log2FC：大部分不显著
    log2fc = np.random.normal(0, 0.3, n_genes)

    # 设定显著上调基因
    up_indices = np.random.choice(n_genes, n_sig_up, replace=False)
    log2fc[up_indices] = np.random.uniform(1.5, 6.0, n_sig_up)

    # 设定显著下调基因（不与上调重叠）
    remaining = list(set(range(n_genes)) - set(up_indices))
    down_indices = np.random.choice(remaining, n_sig_down, replace=False)
    log2fc[down_indices] = np.random.uniform(-6.0, -1.5, n_sig_down)

    # p 值：显著基因有更低的 p 值
    pvalues = np.random.uniform(0, 1, n_genes)
    pvalues[up_indices] = np.random.exponential(0.001, n_sig_up)
    pvalues[down_indices] = np.random.exponential(0.001, n_sig_down)
    pvalues = np.clip(pvalues, 1e-50, 1.0)

    # padj (Benjamini-Hochberg 简化模拟)
    padj = pvalues * np.sort(np.random.uniform(1, 3, n_genes))
    padj = np.clip(padj, 1e-50, 1.0)
    padj[up_indices] = np.random.exponential(0.005, n_sig_up)
    padj[down_indices] = np.random.exponential(0.005, n_sig_down)
    padj = np.clip(padj, 1e-50, 1.0)

    # baseMean (表达量)
    basemean = np.random.lognormal(5, 2, n_genes)
    basemean = np.clip(basemean, 0.1, 1e6)

    df = pd.DataFrame({
        "gene_id": gene_ids,
        "gene_name": "",
        "baseMean": basemean.round(2),
        "log2FoldChange": log2fc.round(4),
        "lfcSE": np.random.uniform(0.05, 0.5, n_genes).round(4),
        "stat": (log2fc / np.random.uniform(0.05, 0.5, n_genes)).round(4),
        "pvalue": pvalues,
        "padj": padj,
    })

    # 给一些基因赋予真实基因名
    name_indices = np.random.choice(n_genes, min(len(real_genes), 100), replace=False)
    for i, idx in enumerate(name_indices[:len(real_genes)]):
        df.loc[idx, "gene_name"] = real_genes[i]

    # 按 padj 排序
    df = df.sort_values("padj").reset_index(drop=True)

    return df


if __name__ == "__main__":
    df = generate_deg_results()
    output_path = "data/sample_deg_results.csv"
    df.to_csv(output_path, index=False)
    print(f"已生成示例数据: {output_path}")
    print(f"总基因数: {len(df)}")
    print(f"显著上调 (|log2FC|>1 & padj<0.05): {len(df[(df['log2FoldChange'] > 1) & (df['padj'] < 0.05)])}")
    print(f"显著下调 (|log2FC|>1 & padj<0.05): {len(df[(df['log2FoldChange'] < -1) & (df['padj'] < 0.05)])}")
    print(f"\n前10行预览:")
    print(df.head(10).to_string(index=False))
