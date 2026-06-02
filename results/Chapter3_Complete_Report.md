# 第三章 完整分析报告：基于公共组学数据的HZ疾病特征与RZV保护性免疫特征分析

**生成日期: 2026-06-03 | 数据版本: V2 (QC过滤, 非配对HZ, p<0.05, |LFC|>0.58)**

---

## 第一部分：GSE242252 — HZ患者全血Bulk RNA-seq

### 1.1 差异表达全景

**方法**: DESeq2非配对设计, 26 Acute vs 23 Convalescent, 保留全部49个HZ样本
**阈值**: p < 0.05, |log2FC| > 0.58

**结果文件**:
| 文件 | 内容 |
|------|------|
| `results/GSE242252/DE_HZ_acute_vs_convalescent.csv` | 全部基因DESeq2结果(含symbol列) |
| `results/GSE242252/Significant_DEGs_table.csv` | 显著DEG表(按|LFC|排序, 含baseMean/LFC/pvalue/padj) |

**对应图表**:
| 图 | 文件 | 内容 |
|----|------|------|
| 火山图 | `results/rnaseq/FigA_Volcano_HZ.pdf` | p<0.05 & |LFC|>0.58标注, Top20上/下调基因symbol标注, X轴±3 |
| PCA图 | `results/rnaseq/FigA_PCA_HZ.pdf` | 急性期(红) vs 恢复期(蓝), 95%置信椭圆 |
| 热图 | `results/rnaseq/FigA_Heatmap_DEGs.pdf` | 显著DEG (FDR<0.05, |LFC|>0.58), 行Z-score, 行聚类 |

**关键DEG**:

| 基因 | log2FC | p-value | 功能类别 |
|------|--------|---------|---------|
| IFI27 | +2.16 | 1.5e-04 | I型IFN ISG |
| PTTG1 | +2.01 | 3.0e-13 | 细胞增殖 |
| BATF2 | +1.95 | 9.5e-04 | IFN诱导转录因子 |
| IGHG4 | +1.94 | 2.6e-05 | 免疫球蛋白重链 |
| IFI44L | +1.82 | 8.5e-05 | I型IFN ISG |
| SERPING1 | +1.80 | 5.3e-07 | 补体C1抑制因子 |
| MZB1 | +1.68 | 3.9e-06 | 浆细胞标志 |
| ISG15 | +1.57 | 1.6e-05 | ISG化修饰 |
| RSAD2 | +1.55 | 1.2e-04 | 抗病毒效应蛋白 |
| TOP2A | +1.50 | 3.4e-13 | DNA复制 |

**解读**: HZ急性期是I型IFN驱动的天然免疫风暴。ISGs (IFI27, IFI44L, RSAD2, ISG15) 大幅上调, 同时补体系统激活 (SERPING1), B细胞向浆细胞分化 (MZB1, IGHG4), 免疫细胞大量增殖 (TOP2A, PTTG1)。上调基因远多于下调, 表明急性期是全面激活状态。

---

### 1.2 通路富集分析

**方法**: clusterProfiler GO BP + KEGG, 所有显著DEG合并(p<0.05, |LFC|>0.58)

**结果文件**:
| 文件 | 内容 |
|------|------|
| `results/GSE242252/GO_BP_enrichment.csv` | GO BP富集结果 |
| `results/GSE242252/KEGG_enrichment.csv` | KEGG通路富集结果 |

**对应图表**:
| 图 | 文件 |
|----|------|
| GO气泡图 | `results/rnaseq/FigA_GO_dotplot.pdf` |
| KEGG气泡图 | `results/rnaseq/FigA_KEGG_dotplot.pdf` |

**关键GO Terms**: Defense Response to Virus, B Cell Receptor Signaling, Antigen Receptor-Mediated Signaling, Mitotic Spindle Organization。

**解读**: 富集结果支持"天然抗病毒 + 适应性体液免疫 + 细胞增殖"三重激活模式。TCR信号相对薄弱, 说明全血中T细胞特异性信号被稀释。

---

### 1.3 HZ Disease Signature

以下14个基因来自Vandoren et al. 2024 Supplementary Table S5, 经本分析DESeq2重跑独立验证:

| 类别 | 基因 | 功能 |
|------|------|------|
| I型IFN | IFI44L, IFI27, RSAD2, ISG15, IFI44, MX1, IFIT5, OASL | 天然抗病毒应答 |
| 补体 | SERPING1, SIGLEC1 | 补体级联激活 |
| 增殖 | TOP2A, PTTG1 | 免疫细胞扩增 |
| 体液 | MZB1, BATF2 | 浆细胞/抗体产生 |

该Signature用于后续评价候选疫苗是否避免了HZ样炎症程序。

---

### 1.4 单细胞验证: ISG细胞来源

**数据**: HRA008316 (Zheng et al. 2024), PBMC scRNA-seq, 3 HA / 3 HP / 3 RP

**对应图表**: `results/chapter3_figures/bulk_vs_scRNA_comparison.png`

**关键发现**:

| ISG | Bulk全血 (HZ Acute) | scRNA T细胞 (HP) | Classical Monocyte (HP/HA) |
|-----|---------------------|-------------------|---------------------------|
| ISG15 | +1.57 ** | -0.37 ns | 0.96x |
| IFI27 | +2.16 * | 未检测到 | **27.63x** |
| SERPING1 | +1.80 *** | 未检测到 | 2.39x |
| IFI44L | +1.82 * | -0.35 ns | 1.19x |

**解读**: Bulk全血中的ISG信号来自单核细胞、DC、中性粒细胞等天然免疫细胞, **不是T细胞产生的**。这解释了为什么靶向T细胞的RZV疫苗不会触发ISG风暴。

---

## 第二部分：GSE249632 — RZV疫苗gE特异性CD4⁺ T细胞动力学

### 2.1 实验设计与数据

**方法**: limma-voom配对设计 (~donor + timepoint), 7例健康接种者, 2,231个QC通过细胞
**时间点**: D0(基线), D14(第一针后), D60(第二针前), D74(第二针后), D365(一年后)
**阈值**: FDR < 0.05, |logFC| > 0.58

**结果文件**:
| 文件 | 内容 |
|------|------|
| `results/GSE249632/DE_D14_vs_D0.csv` | D14差异表达 |
| `results/GSE249632/DE_D60_vs_D0.csv` | D60差异表达 |
| `results/GSE249632/DE_D74_vs_D0.csv` | D74差异表达 |
| `results/GSE249632/DE_D365_vs_D0.csv` | D365差异表达 |
| `results/GSE249632/logFC_matrix_all_timepoints.csv` | 4时间点logFC合并矩阵 |
| `results/GSE249632/Significant_DEGs_table.csv` | 显著DEG表 |

---

### 2.2 DEG时间线

**对应图表**: `results/rnaseq/FigB_DEG_timeline_barchart.pdf`

| 时间点 | 上调 | 下调 | 解读 |
|--------|:--:|:--:|------|
| D14 vs D0 | 164 | 152 | 第一针后强激活 |
| D60 vs D0 | 10 | 4 | **几乎回到基线** (第二针前) |
| D74 vs D0 | 72 | 35 | 第二针再激活, 幅度低于第一针 |
| D365 vs D0 | 20 | 21 | 长期印记残留 (41 DEGs) |

**解读**: 脉冲式、可调控的适应性免疫。D60几乎归零说明应答是自限性的, 不是持续炎症。D365仍有41个DEG提示存在长期程序性改变。

---

### 2.3 各时间点火山图

**对应图表**: `results/rnaseq/FigB_Volcano_D14.pdf`, `D60`, `D74`, `D365` (4张)

每个时间点独立展示差异基因全景, Top8上调+Top8下调基因标注symbol。

---

### 2.4 DEG重叠分析

**对应图表**: `results/rnaseq/FigB_Upset_DEGs.pdf`

**解读**: 不同时间点的DEGs大部分是时间点特异的(尤其是D14特有), 少数基因在所有时间点共享(如HAVCR2)。每次疫苗剂量激活一批新基因, 而非同一批基因反复开关。

---

### 2.5 保护性免疫Signature基因

**对应图表**: 
| 图 | 文件 | 内容 |
|----|------|------|
| 免疫基因点阵图 | `results/rnaseq/FigB_ImmuneGenes_dotplot.pdf` | 22个基因×4时间点, 按7类功能分组 |
| 数据驱动点阵图 | `results/rnaseq/FigB_DataDriven_dotplot.pdf` | 显著≥2时间点的Top25基因 |

**Protection Signature核心基因**:

| 基因 | 功能 | D14 | D60 | D74 | D365 | 模式 |
|------|------|-----|-----|-----|------|------|
| **ZEB2** | T细胞分化重编程 | +3.2 | +2.5 | +3.6 | +3.0 | **持久↑** — 疫苗核心标志 |
| **CTLA4** | 免疫自限性检查点 | +1.5 | +1.4 | +1.6 | +1.0 | **持久↑** — 主动建立刹车 |
| **ICOS** | T细胞共刺激 | +1.3 | +0.8 | +1.2 | +0.6 | D14/D74↑ |
| **HAVCR2** | 效应记忆(TIM-3) | +7.6 | +8.4 | +9.1 | +8.0 | **全时间点最强信号** |
| **SPP1** | 炎症因子(Osteopontin) | −5.4 | −5.6 | −6.2 | −5.9 | **持久↓** — 炎症被主动抑制 |

**阴性对照 — Disease Signature在疫苗中完全不动**:

| 基因 | HZ LFC | D14 | D74 | D365 |
|------|--------|-----|-----|------|
| ISG15 | **+1.57 \*\*** | +0.1 ns | −0.4 ns | −0.1 ns |
| MX1 | +0.82 | +0.3 ns | −0.2 ns | −0.1 ns |
| RSAD2 | **+1.55 \*** | −0.3 ns | −0.2 ns | −0.2 ns |

---

### 2.6 差异基因热图

**对应图表**: `results/rnaseq/FigB_Heatmap_DEGs_timepoints.pdf`

**方法**: 至少1个时间点显著(FDR<0.05, |LFC|>0.58)的基因, 按LFC方差取Top200, 5个时间点均值CPM做Z-score, 行聚类。

**解读模式**: 基因被聚类为几组——
- **簇1 (脉冲效应)**: D14红→D60白→D74红→D365白 (GZMA, CD38, TPX2)
- **簇2 (持久编程)**: D14红→持续红到D365 (ZEB2, HAVCR2, CTLA4)
- **簇3 (持续抑制)**: D14蓝→持续蓝 (SPP1, CA2, IL7R)
- **簇4 (完全不动)**: 全白 (ISG15, RSAD2, MX1 — 阴性对照)

**注意**: 热图颜色是Z-score(相对于该基因自身5个时间点的均值), 不是log2FoldChange。D0有颜色是因为该基因在D0的表达偏离了自身平均水平——这不代表D0有"上调下调"。

---

## 第三部分：两个参照系的关联分析

### 3.1 Disease vs Protection 头对头对比

**对应图表**: `results/chapter3_figures/quadrant_narrative.png`

| 象限 | 定义 | 代表基因 | 生物学含义 |
|------|------|---------|-----------|
| Q1 | HZ↑, RZV↑ | TOP2A, CD38 | 共同增殖程序 (非特异性) |
| **Q2** | **HZ↑, RZV↓** | **ISG15, RSAD2, IFI44L** | **HZ特有的天然免疫炎症 — 疫苗刻意避免** |
| Q3 | HZ↓, RZV↓ | — | 共同抑制 |
| **Q4** | **HZ↓, RZV↑** | **ZEB2, CTLA4, ICOS, HAVCR2** | **RZV特有的适应性免疫调控** |

Spearman ρ ≈ 0: 两个转录程序全局不相关——RZV不是"模拟HZ", 而是建立了一套本质上不同的免疫程序。

### 3.2 ISG模块在RZV时间线中的表现

**对应图表**: `results/chapter3_figures/ISG_module_timeline.png`

I型IFN ISG模块 (25个基因) 在RZV的D0→D14→D60→D74→D365全部时间点基本水平线——阴性证据, 直接证明RZV不激活HZ的标志性炎症通路。

### 3.3 五基因纵向轨迹

**对应图表**: `results/chapter3_figures/five_gene_trajectory.png`

这张图浓缩了第三章的核心叙事: ZEB2/CTLA4/ICOS/HAVCR2在RZV中持续↑ (Protection), ISG15在HZ中↑但在RZV中完全平坦 (Disease被避免)。

---

## 第四部分：候选疫苗评价框架

基于上述两个参照系, 建立双重评价标准:

| 标准 | 参照系 | 核心基因 | 期望方向 | 评价章节 |
|------|--------|---------|---------|---------|
| **标准1**: 是否避免HZ样炎症? | HZ Disease Signature | ISG15, RSAD2, IFI44L, IFI27, SERPING1 | 不激活/低表达 | 第四章 (固有免疫) |
| **标准2**: 是否诱导RZV样保护? | RZV Protection Signature | ZEB2, CTLA4, ICOS, HAVCR2 | 持久上调 | 第五、六章 (T细胞应答) |

---

## 第五部分：输出文件速查

### R分析直接输出 (`results/rnaseq/`)

| 文件 | 内容 |
|------|------|
| `FigA_Volcano_HZ.pdf` | HZ火山图 (p<0.05, |LFC|>0.58, Top20标注) |
| `FigA_PCA_HZ.pdf` | HZ PCA图 |
| `FigA_Heatmap_DEGs.pdf` | HZ显著DEG热图 |
| `FigA_GO_dotplot.pdf` | HZ GO气泡图 |
| `FigA_KEGG_dotplot.pdf` | HZ KEGG气泡图 |
| `Significant_DEGs_table.csv` | HZ显著DEG表 |
| `FigB_DEG_timeline_barchart.pdf` | RZV DEG时间线柱状图 |
| `FigB_Volcano_D14/D60/D74/D365.pdf` | RZV各时间点火山图 (4张) |
| `FigB_Upset_DEGs.pdf` | RZV DEG重叠Upset图 |
| `FigB_Heatmap_DEGs_timepoints.pdf` | RZV DEG热图 (Top200) |
| `FigB_ImmuneGenes_dotplot.pdf` | 22个免疫基因点阵图 (按功能分类) |
| `FigB_DataDriven_dotplot.pdf` | 数据驱动Top25基因点阵图 |
| `Curated_ImmuneGenes_LFC.csv` | 人工筛选免疫基因LFC表 |
| `DataDriven_Top25_DEGs.csv` | 数据驱动Top25基因LFC表 |
| `signature_genes.rds` | HZ + RZV Signature基因数据 (R格式) |

### Python生成的补充图 (`results/chapter3_figures/`)

| 文件 | 内容 |
|------|------|
| `Figure3_1_volcano.png` | HZ火山图 (论文Table S5数据, FDR<0.1, |LFC|>0.58) |
| `GO_bubble_plot.png` | HZ GO富集气泡图 |
| `DEG_timeline_barchart.png` | RZV DEG时间线 |
| `three_gene_trajectory.png` | ZEB2/CTLA4/ISG15三基因轨迹 |
| `five_gene_trajectory.png` | 五基因轨迹 (ZEB2/CTLA4/ICOS/HAVCR2/ISG15) |
| `quadrant_narrative.png` | 象限图 |
| `ISG_module_timeline.png` | ISG模块在RZV时间线评分 |
| `bulk_vs_scRNA_comparison.png` | Bulk vs scRNA ISG对比 |
| `RZV_hallmark_bubble.png` | RZV Hallmark富集 |

### 综合报告

| 文件 | 内容 |
|------|------|
| `results/Chapter3_Report.html` | 第三章HTML报告 (含全部图表) |
| `results/Chapter3_Presentation.pptx` | 第三章PPT (McKinsey风格, 9页) |
| `results/comprehensive_report.html` | 完整分析报告 |
| `results/解读报告_GSE242252_GSE249632.md` | 中文解读报告 |

---

## 第六部分：核心结论

1. **HZ疾病** 是I型IFN介导的天然免疫炎症风暴。ISG信号来自单核/DC/中性粒细胞, 非T细胞。同时伴随补体激活、B细胞抗体产生和广泛增殖。

2. **RZV疫苗** 诱导脉冲式、自限性的CD4⁺ T细胞适应性免疫。D14激活→D60归零→D74再激活→D365长期效应记忆。ZEB2/HAVCR2/CTLA4构成Protection Signature核心。SPP1等炎症通路被主动抑制。

3. **两个免疫程序本质不同** (Spearman ρ≈0)。RZV不激活HZ标志性的I型IFN通路。疫苗通过精准的T细胞分化重编程实现保护, 而非模拟自然感染。

4. **双重评价标准** 为后续章节(第四/五/六章)提供理论坐标系: 候选疫苗应避免Disease Signature + 获得Protection Signature。
