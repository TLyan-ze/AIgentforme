# BioAgent — 生物信息学分析智能体

> 基于 LLM 的多任务生物信息学分析 + 报告生成平台

## 简介

BioAgent 是面向生物信息学分析人员的智能分析平台，覆盖从测序公司交付数据到分析报告生成的完整工作流。支持 **4 类分析任务**，集成大语言模型（LLM）自动撰写结构化中文分析报告。

## 四大分析任务

| 任务 | 命令 | 输入 | 输出 |
|------|------|------|------|
| RNA-seq 差异表达 + 富集 | `rnaseq` | 表达矩阵 + 样本信息 | DESeq2 结果 + GO/KEGG + 报告 |
| 变异检测报告 | `variant` | VCF 文件 | 变异统计 + 功能影响 + 报告 |
| 单细胞分析报告 | `scrna` | 标记基因 + 细胞类型 | 细胞类型解读 + 报告 |
| 功能富集分析 | `enrich` | 基因列表 | GO/KEGG/Reactome + 报告 |

## 快速开始

### 安装

```bash
git clone https://github.com/TLyan-ze/AIgentforme.git
cd AIgentforme
pip install -r requirements.txt
```

### 演示运行

```bash
# RNA-seq 差异表达 + 富集 (支持 human/mouse)
python src/bio_agent.py rnaseq --sample-data --species human

# 变异检测
python src/bio_agent.py variant --sample-data

# 单细胞分析
python src/bio_agent.py scrna --sample-data --species human

# 功能富集
python src/bio_agent.py enrich --sample-data
```

### 使用真实数据

```bash
# RNA-seq: 传入表达矩阵和样本信息
python src/bio_agent.py rnaseq -c counts.csv -m metadata.csv --species mouse

# 变异检测: 传入 VCF 文件
python src/bio_agent.py variant --vcf sample.vcf

# 单细胞: 传入标记基因和细胞汇总
python src/bio_agent.py scrna --markers markers.csv --summary summary.csv

# 富集分析: 传入基因列表
python src/bio_agent.py enrich --genes gene_list.txt --species human
```

### 配置 LLM（可选）

```bash
cp .env.example .env
# 编辑 .env，填入 API Key
# 支持 OpenAI / DeepSeek / 智谱 GLM / 通义千问 等 OpenAI 兼容接口
```

## 项目结构

```
AIgentforme/
├── src/                          # 源代码目录
│   ├── bio_agent.py              # 主入口 CLI（4 个子命令：rnaseq/variant/scrna/enrich）
│   ├── generate_data.py          # 生成所有任务的模拟数据（表达矩阵、VCF、标记基因等）
│   ├── task_rnaseq.py            # Task 1: 完整 RNA-seq Pipeline（质控→DESeq2→富集→可视化→报告）
│   ├── task_variant.py           # Task 2: 变异检测报告（VCF 解析→统计→可视化→报告）
│   ├── task_scrna.py             # Task 3: 单细胞分析报告（标记基因→细胞类型统计→热图→报告）
│   ├── task_enrichment.py        # Task 4: 独立功能富集分析（GO/KEGG/Reactome→可视化→报告）
│   └── __init__.py
├── data/                         # 示例数据目录
│   ├── counts_human.csv          # 人类 RNA-seq 表达矩阵（3000 基因 × 8 样本）
│   ├── metadata_human.csv        # 人类样本信息（sample_id, condition, batch）
│   ├── counts_mouse.csv          # 小鼠 RNA-seq 表达矩阵
│   ├── metadata_mouse.csv        # 小鼠样本信息
│   ├── sample_variants.vcf       # 模拟 VCF 文件（500 个变异位点）
│   ├── scrna_markers_human.csv   # 人类单细胞标记基因结果
│   ├── scrna_summary_human.csv   # 人类单细胞细胞类型汇总
│   ├── scrna_markers_mouse.csv   # 小鼠单细胞标记基因结果
│   ├── scrna_summary_mouse.csv   # 小鼠单细胞细胞类型汇总
│   └── deg_gene_list.txt         # 差异基因列表（用于独立富集分析）
├── examples/                     # 示例输出目录
│   ├── rnaseq_human/             # Task 1 人类 RNA-seq 输出示例
│   │   ├── analysis_report.md    # LLM 生成的分析报告
│   │   ├── volcano_plot.png      # 火山图
│   │   ├── sample_correlation.png # 样本相关性热图
│   │   ├── pca_plot.png          # PCA 降维图
│   │   ├── top_genes_heatmap.png # Top 差异基因热图
│   │   ├── deg_results.csv       # 差异表达结果
│   │   └── statistics.json       # 统计数据
│   ├── variant/                  # Task 2 变异检测输出示例
│   ├── scrna_human/              # Task 3 单细胞输出示例
│   └── enrichment/               # Task 4 富集分析输出示例
├── .env.example                  # LLM API 配置模板
├── requirements.txt              # Python 依赖列表
├── PROJECT_REPORT.pdf            # 项目报告（PDF）
└── README.md                     # 本文件
```

### 重要文件说明

| 文件 | 说明 |
|------|------|
| `src/bio_agent.py` | CLI 主入口，通过 argparse 实现 4 个子命令分发，支持 `--sample-data` 快速演示 |
| `src/task_rnaseq.py` | 核心模块，实现完整 RNA-seq 分析 Pipeline：数据加载→质控→差异表达（pydeseq2/scipy 降级）→GO/KEGG 富集（gseapy）→可视化（5 张图表）→LLM 报告生成 |
| `src/task_variant.py` | VCF 文件解析与变异检测报告生成，包含 SNV/InDel 分类、染色体分布、AF/DP 统计 |
| `src/task_scrna.py` | 单细胞标记基因分析，支持 8 种细胞类型的统计与可视化 |
| `src/task_enrichment.py` | 独立功能富集分析，支持 Human/Mouse 双物种的 GO/KEGG/Reactome 数据库 |
| `src/generate_data.py` | 使用负二项分布生成真实的模拟 RNA-seq 数据，包含人类（TP53、BRCA1 等）和小鼠基因符号 |
| `requirements.txt` | 依赖：pandas, numpy, matplotlib, seaborn, scipy, scikit-learn, openai, gseapy |
| `.env.example` | LLM 配置模板，支持 OpenAI/DeepSeek/智谱 GLM/通义千问等 OpenAI 兼容接口 |

## RNA-seq Pipeline 完整流程

```
测序公司交付数据                BioAgent 自动化
┌──────────────┐
│ 表达矩阵 CSV  │───>  数据加载 & 质控
│ (基因 × 样本)  │     (过滤低表达基因)
├──────────────┤
│ 样本信息 CSV  │───>  DESeq2 差异分析
│ (condition)   │     (上调/下调基因筛选)
└──────────────┘
                       GO/KEGG 功能富集
                       (gseapy / Enrichr)

                       可视化
                       (火山图/PCA/热图/富集图)

                       LLM 报告生成
                       (结构化中文分析报告)
```

## 作者

**Zeqin Yan** — 医疗行业 生物信息工程师

## License

MIT
