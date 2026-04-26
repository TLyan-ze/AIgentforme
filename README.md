# BioReportAgent

**生物信息学差异表达分析报告智能体**

> 基于大语言模型（LLM）的自动化 RNA-seq 差异表达分析报告生成工具

## 简介

BioReportAgent 是面向生物信息学分析人员的智能报告生成工具。它能够自动读取 DESeq2/edgeR 的差异表达分析结果，完成统计分析、可视化绑图和报告撰写，将原本需要 1-2 小时的报告生成工作缩短至分钟级。

## 功能特点

- **自动化统计分析**：一键完成差异基因筛选、分类统计、Top 基因排序
- **专业可视化**：自动生成火山图、MA 图、分布图（PNG 格式）
- **LLM 智能报告**：调用大语言模型自动撰写结构化中文分析报告
- **离线降级**：未配置 API 时自动使用模板生成，确保基本可用
- **命令行接口**：参数灵活，易于集成到 Snakemake 等流程中

## 快速开始

### 安装

```bash
git clone https://github.com/TLyan-ze/AIgentforme.git
cd AIgentforme
pip install -r requirements.txt
```

### 演示运行

```bash
# 使用内置示例数据（离线模式）
python src/bio_report_agent.py --sample-data --no-llm

# 使用示例数据 + LLM 报告
python src/bio_report_agent.py --sample-data
```

### 使用真实数据

```bash
python src/bio_report_agent.py -i your_deg_results.csv -o output -n "项目名称"
```

### 配置 LLM（可选）

```bash
cp .env.example .env
# 编辑 .env，填入 API Key
# 支持 OpenAI / DeepSeek / 智谱 GLM / 通义千问 等 OpenAI 兼容接口
```

## 命令行参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-i, --input` | 差异表达结果 CSV 文件路径 | - |
| `-o, --output` | 输出目录 | `output` |
| `-n, --name` | 项目名称 | `差异表达分析` |
| `--fc` | log2FC 阈值 | `1.0` |
| `--pval` | 校正 p 值阈值 | `0.05` |
| `--sample-data` | 使用内置示例数据 | - |
| `--no-llm` | 不调用 LLM，使用离线模板 | - |

## 输入格式

CSV 文件应包含以下列（DESeq2 标准输出格式）：

```
gene_name,gene_id,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj
```

必要列：`log2FoldChange`, `pvalue`, `padj`
可选列：`gene_name`, `gene_id`, `baseMean`

## 输出内容

运行后在输出目录生成：

| 文件 | 说明 |
|------|------|
| `analysis_report.md` | Markdown 格式分析报告 |
| `volcano_plot.png` | 火山图（标注 Top 基因） |
| `ma_plot.png` | MA 图 |
| `distribution_plots.png` | log2FC 和 p 值分布图 |
| `statistics.json` | 结构化统计数据 |

## 项目结构

```
AIgentforme/
├── src/
│   ├── bio_report_agent.py      # 智能体主程序
│   └── generate_sample_data.py  # 示例数据生成器
├── data/                        # 示例数据
├── output/                      # 输出目录
├── PROJECT_REPORT.md            # 项目报告（详细）
├── requirements.txt             # Python 依赖
└── README.md                    # 本文件
```

## 作者

**颜泽钦** — 医疗行业 生物信息工程师

## License

MIT
