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
src/
├── bio_agent.py           # 主入口 CLI（4 个子命令）
├── generate_data.py       # 生成所有任务的示例数据
├── task_rnaseq.py         # Task 1: RNA-seq DEG + Enrichment
├── task_variant.py        # Task 2: Variant Calling Report
├── task_scrna.py          # Task 3: Single Cell Report
└── task_enrichment.py     # Task 4: Enrichment Analysis
data/                      # 示例数据（表达矩阵、VCF、标记基因等）
examples/                  # 示例输出（报告 + 图表）
```

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
