# RNA-seq 差异表达分析 — RNA-seq 差异表达分析报告

> **生成时间**: 2026-04-26 10:49:06
> **物种**: human
> **分析方法**: t-test (scipy fallback)
> **生成工具**: BioAgent v2.0

---

## 数据质控

| 指标 | 值 |
|------|------|
| 原始基因数 | 3000 |
| 质控后基因数 | 2162 |
| 过滤基因数 | 838 |
| 过滤标准 | counts>=10 in >=2 samples |
| 总 reads 数 | 4,219,501 |

## 差异表达统计

| 指标 | 值 |
|------|------|
| 分析基因数 | 2162 |
| 显著差异基因 | 31 (1.43%) |
| 上调基因 | 27 |
| 下调基因 | 4 |
| log2FC 范围 | [-3.5123, 3.6995] |
| 筛选阈值 | \|log2FC\| > 1.0, padj < 0.05 |

## 可视化

![volcano_plot](../examples/rnaseq_human/volcano_plot.png)
![sample_correlation](../examples/rnaseq_human/sample_correlation.png)
![pca_plot](../examples/rnaseq_human/pca_plot.png)
![top_genes_heatmap](../examples/rnaseq_human/top_genes_heatmap.png)

---

## LLM 分析报告

## 分析报告

> [离线模式] 配置 .env 中的 LLM_API_KEY 可获得 AI 增强报告。

### 概述
本次使用 t-test (scipy fallback) 对 human RNA-seq 数据进行差异表达分析，共鉴定 31 个显著差异基因（上调 27，下调 4）。

### 建议
- 对显著差异基因进行 GO/KEGG 富集分析
- 结合临床表型进一步验证
- 使用 qRT-PCR 验证关键基因

---





---

## Top 差异基因

### 上调 (Top 15)

| 基因 | log2FC | padj |
|------|--------|------|
| ENSG000000126 | 3.3751 | 3.71e-03 |
| ENSG000000142 | 3.0807 | 1.30e-02 |
| CDK6 | 2.4276 | 2.65e-02 |
| ALDH1A1 | 3.0286 | 2.65e-02 |
| ENSG000000117 | 2.9457 | 2.65e-02 |
| ENSG000000110 | 2.6093 | 2.65e-02 |
| CDKN2A | 3.2909 | 2.96e-02 |
| VIM | 2.4237 | 2.96e-02 |
| FOXP3 | 2.4389 | 2.96e-02 |
| MERTK | 2.2801 | 3.05e-02 |

### 下调 (Top 15)

| 基因 | log2FC | padj |
|------|--------|------|
| ENSG000000157 | -2.9176 | 3.54e-02 |
| ENSG000000216 | -1.88 | 3.54e-02 |
| ENSG000000223 | -3.2744 | 3.54e-02 |
| ENSG000000220 | -2.264 | 3.96e-02 |

---
*BioAgent 自动生成 — 2026-04-26 10:49:06*
