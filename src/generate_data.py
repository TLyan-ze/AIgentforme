#!/usr/bin/env python3
"""
为 BioAgent 所有任务生成模拟数据。
- RNA-seq: 表达矩阵 (counts) + 样本信息 (metadata)
- Variant Calling: 模拟 VCF
- Single Cell: 模拟标记基因结果 + 细胞类型
- Enrichment: 模拟基因列表
"""

import os
import numpy as np
import pandas as pd


# ── 真实基因名 ──────────────────────────────────────────────

HUMAN_GENES = [
    "TP53","BRCA1","MYC","EGFR","KRAS","PTEN","CDH1","RB1","PIK3CA","APC",
    "VEGFA","FGFR1","ERBB2","NOTCH1","WNT1","BCL2","MDM2","CDK4","CCND1","JAK2",
    "STAT3","MTOR","AKT1","MAPK1","RAF1","HRAS","NRAS","SOS1","GRB2","SHC1",
    "IRS1","IGF1R","PDGFRA","KIT","FLT3","RET","MET","ALK","ROS1","NTRK1",
    "FGFR2","FGFR3","DDR2","EPHA2","EPHB1","AXL","MERTK","ESR1","PGR","AR",
    "FOXA1","GATA3","CDK6","CCNE1","CDKN1A","CDKN2A","CDKN2B","TP73","TP63",
    "ATM","ATR","CHEK1","CHEK2","BRCA2","PALB2","RAD51","XRCC1","CD44","ALDH1A1",
    "PROM1","EPCAM","KRT5","KRT14","KRT8","KRT18","CDH2","VIM","FN1","SNAI1",
    "SNAI2","TWIST1","ZEB1","ZEB2","MMP2","MMP9","COL1A1","COL3A1","ACTA2",
    "PDGFRB","CXCL8","IL6","TNF","IFNG","CD274","PDCD1","CTLA4","FOXP3",
    "CD8A","CD4","CD19","MS4N1","CD68","ITGAM","LYZ","FCGR3A","NCAM1",
]

MOUSE_GENES = [
    "Trp53","Brca1","Myc","Egfr","Kras","Pten","Cdh1","Rb1","Pik3ca","Apc",
    "Vegfa","Fgfr1","Erbb2","Notch1","Wnt1","Bcl2","Mdm2","Cdk4","Ccnd1","Jak2",
    "Stat3","Mtor","Akt1","Mapk1","Raf1","Hras","Nras","Sos1","Grb2","Shc1",
    "Irs1","Igf1r","Pdgfra","Kit","Flt3","Ret","Met","Alk","Ros1","Ntrk1",
    "Fgfr2","Fgfr3","Ddr2","Epha2","Ephb1","Axl","Mertk","Esr1","Pgr","Ar",
    "Foxa1","Gata3","Cdk6","Ccne1","Cdkn1a","Cdkn2a","Cdkn2b","Trp73","Trp63",
    "Atm","Atr","Chek1","Chek2","Brca2","Palb2","Rad51","Xrcc1","Cd44","Aldh1a1",
    "Prom1","Epcam","Krt5","Krt14","Krt8","Krt18","Cdh2","Vim","Fn1","Snai1",
    "Snai2","Twist1","Zeb1","Zeb2","Mmp2","Mmp9","Col1a1","Col3a1","Acta2",
    "Pdgfrb","Cxcl8","Il6","Tnf","Ifng","Cd274","Pdcd1","Ctla4","Foxp3",
    "Cd8a","Cd4","Cd19","Ms4a1","Cd68","Itgam","Lyz","Fcgr3a","Ncam1",
]


def _neg_binom_sample(mu: float, alpha: float = 0.1) -> int:
    """从负二项分布采样 RNA-seq 计数。"""
    mu = max(mu, 0.01)
    r = 1.0 / alpha
    p = r / (r + mu)
    return np.random.negative_binomial(max(int(r), 1), min(p, 0.999))


def generate_rnaseq_data(output_dir: str, species: str = "human",
                         n_genes: int = 3000, n_control: int = 4,
                         n_treat: int = 4, n_deg_up: int = 150,
                         n_deg_down: int = 120, seed: int = 42):
    """生成 RNA-seq 表达矩阵和样本信息。

    Parameters
    ----------
    output_dir : str
        输出目录
    species : str
        'human' 或 'mouse'
    n_genes : int
        总基因数
    n_control, n_treat : int
        对照组/处理组样本数
    n_deg_up, n_deg_down : int
        模拟上调/下调基因数
    """
    np.random.seed(seed)
    os.makedirs(output_dir, exist_ok=True)

    gene_pool = HUMAN_GENES if species == "human" else MOUSE_GENES
    prefix = "ENSG" if species == "human" else "ENSMUSG"

    # 生成基因 ID
    gene_ids = [f"{prefix}{i:09d}" for i in range(1, n_genes + 1)]
    gene_names = [""] * n_genes

    # 为前 len(gene_pool) 个基因分配真实名称
    for i in range(min(len(gene_pool), n_genes)):
        gene_names[i] = gene_pool[i]

    # 差异表达基因索引
    deg_up_idx = set(range(0, n_deg_up))
    deg_down_idx = set(range(n_deg_up, n_deg_up + n_deg_down))

    # 样本名
    ctrl_samples = [f"Control_{i+1}" for i in range(n_control)]
    treat_samples = [f"Treat_{i+1}" for i in range(n_treat)]
    all_samples = ctrl_samples + treat_samples

    # 生成计数矩阵
    counts = np.zeros((n_genes, len(all_samples)), dtype=int)
    base_means = np.random.lognormal(3, 2, n_genes)

    for i in range(n_genes):
        for j, sample in enumerate(all_samples):
            mu = base_means[i]
            if j >= n_control:  # 处理组
                if i in deg_up_idx:
                    mu *= np.random.uniform(2.5, 8.0)
                elif i in deg_down_idx:
                    mu *= np.random.uniform(0.05, 0.4)
            counts[i, j] = _neg_binom_sample(mu)

    # 构建 DataFrame
    counts_df = pd.DataFrame(counts, index=gene_ids, columns=all_samples)
    counts_df.insert(0, "gene_name", gene_names)
    counts_df.index.name = "gene_id"

    # 样本信息
    metadata_df = pd.DataFrame({
        "sample_id": all_samples,
        "condition": ["control"] * n_control + ["treatment"] * n_treat,
        "batch": ["batch1"] * (n_control // 2 + n_treat // 2) +
                 ["batch2"] * (n_control - n_control // 2 + n_treat - n_treat // 2),
    })
    metadata_df.set_index("sample_id", inplace=True)

    counts_path = os.path.join(output_dir, f"counts_{species}.csv")
    meta_path = os.path.join(output_dir, f"metadata_{species}.csv")

    counts_df.to_csv(counts_path)
    metadata_df.to_csv(meta_path)

    print(f"[RNA-seq] 物种: {species}")
    print(f"  表达矩阵: {counts_path}  ({n_genes} 基因 × {len(all_samples)} 样本)")
    print(f"  样本信息: {meta_path}")
    print(f"  模拟上调: {n_deg_up}, 下调: {n_deg_down}")
    return counts_path, meta_path


def generate_vcf_data(output_dir: str, n_variants: int = 500, seed: int = 42):
    """生成模拟 VCF 文件。"""
    np.random.seed(seed)
    os.makedirs(output_dir, exist_ok=True)

    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    ref_bases = ["A", "T", "C", "G"]
    variant_types = ["SNV", "SNV", "SNV", "SNV", "InDel", "InDel"]  # SNV 更常见
    impacts = ["MODIFIER", "MODIFIER", "MODERATE", "MODERATE", "HIGH"]

    records = []
    for i in range(n_variants):
        chrom = np.random.choice(chroms)
        pos = np.random.randint(10000, 250000000)
        ref = np.random.choice(ref_bases)
        vtype = np.random.choice(variant_types)
        if vtype == "SNV":
            alt = np.random.choice([b for b in ref_bases if b != ref])
        else:
            alt = ref + "".join(np.random.choice(ref_bases, np.random.randint(1, 4)))
        qual = round(np.random.uniform(20, 200), 1)
        dp = np.random.randint(5, 200)
        af = round(np.random.uniform(0.01, 1.0), 4)
        impact = np.random.choice(impacts)
        gene = np.random.choice(HUMAN_GENES)
        records.append([chrom, pos, ".", ref, alt, qual, "PASS",
                         f"DP={dp};AF={af};TYPE={vtype};IMPACT={impact};GENE={gene}"])

    vcf_path = os.path.join(output_dir, "sample_variants.vcf")
    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
        f.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
        f.write('##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant Type">\n')
        f.write('##INFO=<ID=IMPACT,Number=1,Type=String,Description="Functional Impact">\n')
        f.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Gene Name">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for rec in records:
            f.write("\t".join(str(x) for x in rec) + "\n")

    print(f"[Variant] VCF: {vcf_path}  ({n_variants} 变异)")
    return vcf_path


def generate_scrna_data(output_dir: str, species: str = "human", seed: int = 42):
    """生成模拟单细胞标记基因结果和细胞类型。"""
    np.random.seed(seed)
    os.makedirs(output_dir, exist_ok=True)

    gene_pool = HUMAN_GENES if species == "human" else MOUSE_GENES

    cell_types = ["T cell", "B cell", "Macrophage", "NK cell", "Epithelial",
                  "Fibroblast", "Endothelial", "Dendritic cell"]
    markers_per_type = {}
    used_genes = set()

    for ct in cell_types:
        available = [g for g in gene_pool if g not in used_genes]
        n_markers = np.random.randint(8, 20)
        if len(available) >= n_markers:
            markers = np.random.choice(available, n_markers, replace=False).tolist()
        else:
            markers = available
        markers_per_type[ct] = markers
        used_genes.update(markers)

    # 生成 marker 基因结果
    rows = []
    for ct, markers in markers_per_type.items():
        for gene in markers:
            rows.append({
                "gene": gene,
                "cluster": ct,
                "avg_log2FC": round(np.random.uniform(0.8, 4.5), 4),
                "pct_in": round(np.random.uniform(30, 95), 1),
                "pct_out": round(np.random.uniform(0.5, 15), 1),
                "pvalue": f"{np.random.uniform(1e-50, 1e-5):.2e}",
                "padj": f"{np.random.uniform(1e-40, 1e-3):.2e}",
            })

    marker_df = pd.DataFrame(rows)
    marker_path = os.path.join(output_dir, f"scrna_markers_{species}.csv")
    marker_df.to_csv(marker_path, index=False)

    # 细胞类型汇总
    summary_df = pd.DataFrame([
        {"cell_type": ct, "n_cells": np.random.randint(200, 3000),
         "top_markers": ", ".join(markers_per_type[ct][:3])}
        for ct in cell_types
    ])
    summary_path = os.path.join(output_dir, f"scrna_summary_{species}.csv")
    summary_df.to_csv(summary_path, index=False)

    print(f"[scRNA] 物种: {species}")
    print(f"  标记基因: {marker_path}  ({len(marker_df)} 条记录)")
    print(f"  细胞汇总: {summary_path}  ({len(cell_types)} 种细胞类型)")
    return marker_path, summary_path


def generate_gene_list(output_dir: str, n_genes: int = 200, seed: int = 42):
    """生成模拟差异基因列表（用于独立富集分析）。"""
    np.random.seed(seed)
    os.makedirs(output_dir, exist_ok=True)

    genes = list(np.random.choice(HUMAN_GENES + MOUSE_GENES, n_genes, replace=False))
    path = os.path.join(output_dir, "deg_gene_list.txt")
    with open(path, "w") as f:
        f.write("\n".join(genes))

    print(f"[GeneList] 基因列表: {path}  ({len(genes)} 基因)")
    return path


if __name__ == "__main__":
    data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
    print("=" * 60)
    print("BioAgent — 生成所有示例数据")
    print("=" * 60)

    # RNA-seq (人 + 小鼠)
    generate_rnaseq_data(data_dir, species="human")
    generate_rnaseq_data(data_dir, species="mouse")

    # Variant Calling
    generate_vcf_data(data_dir)

    # Single Cell
    generate_scrna_data(data_dir, species="human")
    generate_scrna_data(data_dir, species="mouse")

    # Gene List (独立富集分析用)
    generate_gene_list(data_dir)

    print("\n所有示例数据生成完毕。")
