# BioAgent — 生物信息学分析智能体

## 项目报告

**作者**：颜泽钦
**单位**：医疗行业 生物信息分析部门
**日期**：2026年4月
**GitHub**：https://github.com/TLyan-ze/AIgentforme

---

## 一、单位与工作介绍

### 1.1 工作背景

本人就职于医疗行业科研部门，从事生物信息学分析与研究工作。所在部门配备高通量测序分析平台、单细胞测序分析平台、生物信息高性能计算集群等设施，承担基础与转化医学研究的技术支撑工作。

### 1.2 个人工作介绍

我是颜泽钦，现任医疗行业生物信息工程师。日常工作流程通常从测序公司交付原始数据开始，经过质控、比对、定量、差异分析、功能富集等一系列步骤，最终输出分析报告。主要工作内容包括：

- **RNA-seq 转录组分析**：从表达矩阵出发，使用 DESeq2 进行差异表达分析，GO/KEGG 通路富集分析，每月处理 5-10 个项目
- **基因组变异检测**：WES/WGS 数据的变异检测（SNV、InDel），功能影响注释
- **单细胞 RNA-seq 分析**：使用 Seurat/scanpy 进行细胞聚类、标记基因鉴定、细胞类型注释
- **流程搭建与维护**：基于 Snakemake 搭建自动化分析流程

---

## 二、有待改进的工作流程

### 2.1 当前流程：从测序数据到分析报告

测序公司交付数据后，典型的 RNA-seq 分析工作流如下：

```
测序公司交付                    人工操作                      最终交付
┌──────────┐    ┌──────────────┐    ┌───────────────┐    ┌──────────┐
│ 原始数据  │───>│ 质控 & 比对   │───>│ 差异表达分析   │───>│ 分析报告  │
│ (fastq)   │    │ (已有流程)   │    │ (DESeq2/R)    │    │ (Word)   │
└──────────┘    └──────────────┘    └───────────────┘    └──────────┘
                                         │
                                         ▼
                                   ┌───────────────┐    ┌──────────┐
                                   │ 功能富集分析   │───>│ 整合撰写  │
                                   │ (clusterPr.)   │    │ 分析报告  │
                                   └───────────────┘    └──────────┘
```

**当前痛点**（以下步骤可由 AI 自动化）：

| 步骤 | 当前方式 | 耗时/次 | 问题 |
|------|---------|--------|------|
| 数据质控统计 | 手动运行脚本 + 人工整理 | 15 min | 重复性高 |
| 差异表达分析 | 手动写 R 脚本 / 每项目重写 | 20 min | 代码不统一 |
| 功能富集分析 | 手动调用 clusterProfiler | 15 min | 重复劳动 |
| 可视化绑图 | 每次重写 matplotlib/R 代码 | 20 min | 格式不统一 |
| 撰写分析报告 | 人工阅读结果 + 写 Word | 40-60 min | 主观性强，质量参差 |
| **合计** | | **110-130 min** | |

### 2.2 工作流程中可使用 LLM/AI 替代的部分

| 步骤 | AI 替代方案 | 可替代程度 |
|------|-----------|-----------|
| 数据质控 → 统计汇总 | Python 自动化 + LLM 解读 | ★★★★★ 完全替代 |
| 差异表达分析 | DESeq2 (pydeseq2) 封装 | ★★★★★ 完全替代 |
| GO/KEGG 富集 | gseapy 自动化 | ★★★★★ 完全替代 |
| 可视化绑图 | 封装好的绑图函数 | ★★★★★ 完全替代 |
| 撰写分析报告 | LLM 基于结果自动生成 | ★★★★☆ 大部分替代 |
| 生物学解读 | LLM 结合知识库 | ★★★☆☆ 需人工审核 |

---

## 三、智能体设计

### 3.1 系统架构

```
                    BioAgent 系统架构 (v2.0)

┌─────────────────────────────────────────────────────────────┐
│                        CLI 入口层                            │
│   bio_agent.py rnaseq / variant / scrna / enrich            │
└────────────────────────┬────────────────────────────────────┘
                         │
         ┌───────────────┼───────────────┬──────────────┐
         ▼               ▼               ▼              ▼
   ┌───────────┐   ┌───────────┐  ┌───────────┐  ┌───────────┐
   │  Task 1   │   │  Task 2   │  │  Task 3   │  │  Task 4   │
   │  RNA-seq  │   │  Variant  │  │  scRNA    │  │  Enrich   │
   │  DEG+富集  │   │  变异检测  │  │  单细胞   │  │  功能富集  │
   └─────┬─────┘   └─────┬─────┘  └─────┬─────┘  └─────┬─────┘
         │               │               │              │
         ▼               ▼               ▼              ▼
   ┌─────────────────────────────────────────────────────────┐
   │                   共享分析引擎                            │
   │  DESeq2/pydeseq2 | gseapy | scipy | sklearn            │
   │  matplotlib | seaborn                                   │
   └────────────────────────┬────────────────────────────────┘
                            │
                            ▼
                   ┌──────────────────┐
                   │   LLM 报告生成    │
                   │  (OpenAI 兼容)   │
                   │  离线模板降级     │
                   └──────────────────┘
                            │
                            ▼
                   ┌──────────────────┐
                   │  Markdown 报告    │
                   │  PNG 可视化图表   │
                   │  JSON 统计数据    │
                   └──────────────────┘
```

### 3.2 四大任务模块

#### Task 1: RNA-seq 差异表达 + 功能富集（核心任务）

**完整 Pipeline**：表达矩阵 → 质控 → DESeq2 → GO/KEGG 富集 → 可视化 → LLM 报告

- **数据输入**：测序公司交付的基因表达矩阵 (CSV) + 样本信息 (CSV)
- **物种支持**：Human (GRCh38) / Mouse (GRCm38)
- **分析步骤**：
  1. 数据加载与质控（过滤低表达基因）
  2. DESeq2 差异表达分析（pydeseq2，scipy t-test 降级）
  3. GO/KEGG 功能富集分析（gseapy/Enrichr）
  4. 可视化：火山图、样本相关性热图、PCA、Top 基因热图、富集气泡图
  5. LLM 生成结构化中文分析报告

#### Task 2: 变异检测报告

- **数据输入**：VCF 文件
- **分析步骤**：解析 VCF → 变异分类统计 → 功能影响注释 → 染色体分布 → AF/DP 分布 → LLM 报告

#### Task 3: 单细胞 RNA-seq 分析报告

- **数据输入**：标记基因结果 (CSV) + 细胞类型汇总 (CSV)
- **分析步骤**：细胞类型统计 → 标记基因分析 → 共享标记基因 → 热图 → LLM 报告

#### Task 4: 独立功能富集分析

- **数据输入**：基因列表 (TXT，每行一个基因)
- **分析步骤**：GO/KEGG/Reactome 富集 → 可视化 → LLM 报告

### 3.3 技术栈

| 类别 | 技术 |
|------|------|
| 语言 | Python 3.8+ |
| 差异表达 | pydeseq2 (DESeq2) / scipy (降级) |
| 功能富集 | gseapy (Enrichr API) |
| 可视化 | matplotlib, seaborn |
| 降维 | scikit-learn (PCA) |
| LLM 接口 | OpenAI Python SDK (兼容多后端) |
| 数据处理 | pandas, numpy |

---

## 四、实现细节

### 4.1 代码结构

```
AIgentforme/
├── src/
│   ├── bio_agent.py           # 主入口 CLI（4 个子命令）
│   ├── generate_data.py       # 生成所有任务的模拟数据
│   ├── task_rnaseq.py         # Task 1: 完整 RNA-seq Pipeline
│   ├── task_variant.py        # Task 2: 变异检测报告
│   ├── task_scrna.py          # Task 3: 单细胞分析报告
│   ├── task_enrichment.py     # Task 4: 独立富集分析
│   └── __init__.py
├── data/                      # 示例数据
│   ├── counts_human.csv       # 人类表达矩阵 (3000 基因 × 8 样本)
│   ├── metadata_human.csv     # 人类样本信息
│   ├── counts_mouse.csv       # 小鼠表达矩阵
│   ├── metadata_mouse.csv     # 小鼠样本信息
│   ├── sample_variants.vcf    # 模拟 VCF
│   ├── scrna_markers_*.csv    # 单细胞标记基因
│   ├── scrna_summary_*.csv    # 单细胞细胞类型汇总
│   └── deg_gene_list.txt      # 差异基因列表
├── examples/                  # 示例输出
├── .env.example
├── requirements.txt
├── README.md
└── PROJECT_REPORT.md          # 本文件
```

### 4.2 Task 1 核心流程详解

```python
# 从表达矩阵到报告的完整流程
result = run_pipeline(
    counts_path="data/counts_human.csv",    # 测序公司交付的表达矩阵
    metadata_path="data/metadata_human.csv", # 样本分组信息
    species="human",                         # 支持 human / mouse
    fc_threshold=1.0,                        # log2FC 阈值
    pval_threshold=0.05,                     # padj 阈值
)
# 输出: analysis_report.md + 4张可视化图表 + statistics.json + deg_results.csv
```

**Step 1 — 数据加载与质控**：
- 读取表达矩阵（基因 × 样本）和样本信息（sample_id, condition）
- 过滤低表达基因（默认：counts >= 10 且在 >= 2 个样本中）

**Step 2 — 差异表达分析**：
- 优先使用 pydeseq2（DESeq2 的 Python 实现）
- 未安装时自动降级为 scipy t-test + BH 校正
- 输出：log2FoldChange, pvalue, padj

**Step 3 — 功能富集分析**：
- 使用 gseapy 调用 Enrichr API
- Human: GO_BP_2023 + KEGG_2021_Human
- Mouse: GO_BP_2023 + KEGG_2021_Mouse
- 未安装 gseapy 或网络不可用时跳过

**Step 4 — 可视化**：
- 火山图（标注 Top 基因）
- 样本相关性热图
- PCA 降维图
- Top 差异基因热图
- GO/KEGG 富集条形图

**Step 5 — LLM 报告生成**：
- System Prompt：设定"资深生物信息学分析师"角色
- User Prompt：注入统计数据、Top 基因、富集通路
- Temperature：0.3（确保专业性）
- 离线降级：未配置 API 时使用模板生成基本框架

### 4.3 LLM 集成与降级策略

```
LLM 可用？
  ├── 是 → 调用 OpenAI 兼容 API → 生成完整分析报告
  └── 否 → 使用离线模板 → 生成基本框架 + 统计数据
```

支持所有 OpenAI 兼容后端：OpenAI、DeepSeek、智谱 GLM、通义千问等。

---

## 五、项目展示

### 5.1 使用过程

#### 环境准备

```bash
git clone https://github.com/TLyan-ze/AIgentforme.git
cd AIgentforme
pip install -r requirements.txt
```

#### 运行各任务

```bash
# Task 1: RNA-seq 差异表达 + 富集 (人类)
python src/bio_agent.py rnaseq --sample-data --species human

# Task 1: RNA-seq (小鼠)
python src/bio_agent.py rnaseq --sample-data --species mouse

# Task 2: 变异检测
python src/bio_agent.py variant --sample-data

# Task 3: 单细胞
python src/bio_agent.py scrna --sample-data --species human

# Task 4: 独立富集
python src/bio_agent.py enrich --sample-data
```

#### 查看结果

```bash
# RNA-seq 输出
ls output/
# analysis_report.md    volcano_plot.png    sample_correlation.png
# pca_plot.png          top_genes_heatmap.png
# deg_results.csv       statistics.json
```

### 5.2 实验数据与结果

#### Task 1: RNA-seq 差异表达 + 富集

**输入**：模拟人类 RNA-seq 表达矩阵（3000 基因 × 8 样本：4 control + 4 treatment）

| 分析环节 | 结果 |
|---------|------|
| 质控 | 3000 → 2162 基因（过滤 838 低表达） |
| 差异分析 | 31 个显著差异基因（上调 27，下调 4） |
| 方法 | scipy t-test + BH 校正（pydeseq2 可选） |
| 可视化 | 火山图、相关性热图、PCA、Top 基因热图 |
| 报告 | analysis_report.md |

**可视化输出**：
- `volcano_plot.png`：火山图，标注 Top 差异基因
- `sample_correlation.png`：样本相关性热图，验证分组合理性
- `pca_plot.png`：PCA 降维图，展示样本聚类
- `top_genes_heatmap.png`：Top 30 差异基因表达热图

#### Task 2: 变异检测

**输入**：模拟 VCF 文件（500 个变异位点）

| 指标 | 值 |
|------|-----|
| 总变异数 | 500 |
| SNV | 332 |
| InDel | 168 |
| HIGH impact | 变异功能分级统计 |
| 可视化 | 变异类型分布、染色体分布、AF/DP 分布 |

#### Task 3: 单细胞分析

**输入**：模拟标记基因（106 条）+ 8 种细胞类型

| 指标 | 值 |
|------|-----|
| 细胞类型 | 8 种 |
| 总细胞数 | 15,991 |
| 标记基因 | 106 条 |
| 可视化 | 细胞分布、标记基因热图 |

### 5.3 使用效果评估

#### 效率对比

| 指标 | 人工流程 | BioAgent | 提升幅度 |
|------|---------|----------|---------|
| 数据质控 | 15 min | < 1s | ~900x |
| 差异分析 | 20 min | ~1s | ~1200x |
| 功能富集 | 15 min | ~2s (网络) | ~450x |
| 可视化 | 20 min | ~3s | ~400x |
| 报告撰写 | 40-60 min | 5-15s (LLM) | ~240x |
| **总计** | **110-130 min** | **< 30s** | **~250x** |

#### 测试结果

| 测试场景 | 状态 | 耗时 |
|---------|------|------|
| RNA-seq (human, 离线模式) | 通过 | ~1.3s |
| RNA-seq (mouse, 离线模式) | 通过 | ~1.3s |
| Variant Calling | 通过 | ~0.2s |
| scRNA-seq (human) | 通过 | ~0.2s |
| Enrichment (需网络) | 通过 (需 gseapy) | ~5-10s |

---

## 六、讨论与展望

### 6.1 优势

1. **覆盖完整工作流**：从测序公司交付数据到分析报告，一站式完成
2. **多任务支持**：差异表达、变异检测、单细胞、功能富集四大类分析
3. **物种支持**：Human / Mouse 双物种
4. **LLM 增强**：大语言模型自动撰写专业分析报告
5. **降级策略**：pydeseq2 → scipy，gseapy → 跳过，LLM → 离线模板，层层保障
6. **易于集成**：CLI 设计，可嵌入 Snakemake / Nextflow 等流程框架

### 6.2 局限性

1. pydeseq2 的 Python 实现与 R 版 DESeq2 存在细微差异
2. gseapy 富集分析依赖网络连接
3. LLM 可能产生幻觉，生成的生物学解读需人工审核
4. 当前单细胞任务基于预处理后的标记基因，未包含上游聚类分析

### 6.3 后续改进方向

1. **pydeseq2 集成优化**：确保与 R DESeq2 结果的一致性
2. **多格式输出**：支持 Word (.docx) 和 PDF
3. **交互式报告**：Plotly/Dash 生成交互式 HTML 报告
4. **RAG 增强**：构建本地文献知识库，提升 LLM 解读准确性
5. **流程编排**：与 Snakemake 深度集成，从 fastq 到报告的全自动流程

---

## 七、总结

本项目基于我在医疗行业的实际工作需求，设计并实现了 BioAgent——一个多任务生物信息学分析智能体。该智能体覆盖了 RNA-seq 差异表达分析（从表达矩阵到报告）、变异检测、单细胞分析、功能富集分析四大核心任务，将整个分析报告生成流程从 2 小时缩短至分钟级，有效提升了工作效率和分析质量的一致性。

---

*颜泽钦 — 医疗行业 生物信息工程师*
*2026年4月*
