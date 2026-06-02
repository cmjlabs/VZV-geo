# 第三章 基于公共转录组数据构建HZ疫苗免疫评价的二维参照系

**数据版本: V2 | 日期: 2026-06-03**

---

## 3.1 研究设计与核心逻辑

### 3.1.1 核心问题

理想的HZ候选疫苗应诱导什么样的免疫应答？直觉上，"模拟自然感染"似乎合理——但HZ本身就是VZV再激活引起的疾病。我们需要更精确的标准。

### 3.1.2 分析策略

本章不试图比较HZ与RZV（两者细胞来源和测序平台不同，直接比较无意义）。取而代之的是，利用两个互补的数据集各自的信息优势，分别提取两种**独立的免疫参照系**：

| 数据集 | 样本 | 技术 | 信息优势 | 回答的问题 |
|--------|------|------|---------|-----------|
| **GSE242252** | 26例HZ急性期 vs 23例恢复期 | 全血Bulk RNA-seq | 捕捉系统性免疫全景，包含所有免疫细胞群体 | 哪些固有免疫和系统炎症特征需要**避免**？→ **Disease Reference** |
| **GSE249632** | 7例RZV接种者, 2,231个gE-CD4⁺ T细胞 | scRNA-seq (SMART-Seq v4) | 高分辨率追踪抗原特异性T细胞的内在分化程序 | 成功的T细胞应答应具备什么**质量**？→ **Protection Reference** |

> **方法学立场**: 两数据集分别锚定免疫评价的两个独立维度。HZ全血揭示"广度"——哪些通路在疾病中被全身性激活；RZV T细胞揭示"精度"——成功的疫苗如何在抗原特异性T细胞中建立持久编程。两者不直接比较，而是联合定义一个二维评价空间。

### 3.1.3 二维评价框架

```
                    适应性免疫质量 (T细胞编程)
                    需要获得 Protection Reference
                              ↑
                              │
          ┌───────────────────┼───────────────────┐
          │  低质量适应性免疫   │  高质量适应性免疫    │
          │  低固有免疫炎症     │  低固有免疫炎症      │
          │  (免疫不足)        │  ★ 理想疫苗         │
          ├───────────────────┼───────────────────┤
          │  低质量适应性免疫   │  高质量适应性免疫    │
          │  高固有免疫炎症     │  高固有免疫炎症      │
          │  (最差情况)        │  (类似HZ, 不安全)   │
          └───────────────────┼───────────────────┘
                              │
              固有免疫炎症 (系统性ISG)
              需要避免 Disease Reference
```

后续第四章沿纵轴（固有免疫）评价候选疫苗，第五/六章沿横轴（适应性免疫质量）评价。

---

## 3.2 HZ急性期系统性抗病毒转录特征（GSE242252）

### 3.2.1 分析设计与方法

| 项目 | 详情 |
|------|------|
| 数据来源 | GSE242252 (Vandoren et al. 2024, *J Infect Dis*) |
| 样本 | 26例HZ急性期全血 vs 23例恢复期全血 |
| 分析方法 | DESeq2, 非配对设计 (~timepoint) |
| 差异基因阈值 | p < 0.05 (nominal), |log2FC| > 0.58 |
| 注 | 采用名义p值而非FDR, 因样本量限制FDR过度保守; 结合|LFC|阈值双重过滤; 所有关键结论经原文Table S5独立验证一致 |

### 3.2.2 差异表达整体特征

HZ急性期 vs 恢复期共鉴定出 **352个上调、44个下调**基因 (p<0.05, |LFC|>0.58)。上调基因数量远超下调, 表明急性期以免疫激活为主。

**数据与图表**:
- 火山图: `results/rnaseq/FigA_Volcano_HZ.pdf`
- 显著基因表: `results/GSE242252/Significant_DEGs_table.csv`
- 完整结果: `results/GSE242252/DE_HZ_acute_vs_convalescent.csv`

**Top上调基因**:

| 基因 | log2FC | p-value | 功能 |
|------|--------|---------|------|
| IFI27 | +2.16 | 1.5e-04 | I型IFN诱导抗病毒ISG |
| PTTG1 | +2.01 | 3.0e-13 | 有丝分裂调控 |
| BATF2 | +1.95 | 9.5e-04 | IFN诱导转录因子 |
| IGHG4 | +1.94 | 2.6e-05 | 免疫球蛋白重链 (IgG4) |
| IFI44L | +1.82 | 8.5e-05 | I型IFN诱导抗病毒ISG |
| SERPING1 | +1.80 | 5.3e-07 | 补体C1抑制因子 |
| MZB1 | +1.68 | 3.9e-06 | 浆细胞分化标志 |
| ISG15 | +1.57 | 1.6e-05 | ISG化修饰/抗病毒 |
| RSAD2 | +1.55 | 1.2e-04 | 抗病毒效应蛋白 (Viperin) |
| TOP2A | +1.50 | 3.4e-13 | DNA拓扑异构酶/DNA复制 |

### 3.2.3 功能富集

**数据与图表**:
- GO BP: `results/rnaseq/FigA_GO_dotplot.pdf`, `results/GSE242252/GO_BP_enrichment.csv`
- KEGG: `results/rnaseq/FigA_KEGG_dotplot.pdf`, `results/GSE242252/KEGG_enrichment.csv`

GO和KEGG富集分析 (所有显著DEG合并) 突出以下生物学过程:
- Defense Response to Virus
- B Cell Receptor Signaling Pathway
- Antigen Receptor-Mediated Signaling
- Mitotic Spindle Organization

### 3.2.4 ISGs的细胞来源评估

**数据与图表**: `results/chapter3_figures/bulk_vs_scRNA_comparison.png`

HRA008316 (Zheng et al. 2024, PBMC scRNA-seq) 提供了细胞类型分辨率:

| 基因 | Bulk全血 (HZ Acute) | scRNA T细胞 (HP) | Classical Monocyte (HP/HA比) |
|------|---------------------|-------------------|------------------------------|
| IFI27 | +2.16 * | 未检测到 | **27.6x** |
| SERPING1 | +1.80 *** | 未检测到 | 2.4x |
| ISG15 | +1.57 ** | -0.37 ns | 0.96x |
| IFI44L | +1.82 * | -0.35 ns | 1.2x |

该分析表明: (i) 全血中观察到的ISG信号在PBMC scRNA-seq的T细胞中未出现, (ii) 单核细胞可能是ISG信号的主要贡献者之一。鉴于全血ISG来源的多样性 (中性粒细胞、pDC、DC、NK细胞均有贡献), 此处不做唯一来源推论。

### 3.2.5 小节

| 结论 | 内容 |
|------|------|
| **HZ急性期全血特征** | I型IFN相关系统性抗病毒程序激活, 伴随补体活化、B细胞抗体应答和免疫细胞增殖 |
| **Disease Reference** | ISG15, RSAD2, IFI44L, IFI27, SERPING1, MZB1 — 代表HZ急性期中被广泛激活的免疫程序 |
| **用途** | 作为**固有免疫激活"广度"的参照** — 候选疫苗不应在全身层面复现此广泛炎症特征 |

> **关键限定**: 该分析识别的是HZ疾病中的*系统性转录特征*, 不预设这些特征在疾病中是"有害的"或"病理性"的。它们在宿主抗病毒防御中具有保护功能。其作为参照系的意义在于: 疫苗应能在TLR4-MyD88通路启动免疫的同时, 避免触发RIG-I/MDA5→IRF→ISG级联反应。

---

## 3.3 RZV疫苗gE特异性CD4⁺ T细胞的转录编程（GSE249632）

### 3.3.1 分析设计与方法

| 项目 | 详情 |
|------|------|
| 数据来源 | GSE249632 (GEO) |
| 样本 | 7例健康接种者, 四聚体分选gE⁺CD4⁺ T细胞, 2,231个QC通过细胞 |
| 时间点 | D0 (基线), D14 (第一针后), D60 (第二针前), D74 (第二针后), D365 (一年后) |
| 分析方法 | Pseudobulk聚合 (donor × timepoint), limma-voom配对设计 (~donor + timepoint) |
| 阈值 | FDR < 0.05, |logFC| > 0.58 |

> **分析限定**: 本数据集仅包含gE特异性CD4⁺ T细胞, 因此所有结论严格限定于该细胞群体。关于RZV在其他细胞类型或全身水平是否激活固有免疫通路, 本章不做推断。

### 3.3.2 差异基因时间线: 脉冲式激活模式

**数据与图表**:
- DEG时间线柱状图: `results/rnaseq/FigB_DEG_timeline_barchart.pdf`
- 各时间点火山图: `results/rnaseq/FigB_Volcano_D14/D60/D74/D365.pdf`
- DEG重叠分析: `results/rnaseq/FigB_Upset_DEGs.pdf`

| 时间点 | 上调 | 下调 | 模式 |
|--------|:--:|:--:|------|
| D14 vs D0 | 164 | 152 | 第一针强转录应答 |
| D60 vs D0 | 10 | 4 | 应答消退, 接近基线 |
| D74 vs D0 | 72 | 35 | 第二针再激活 (幅度低于D14) |
| D365 vs D0 | 20 | 21 | 长期转录印记残留 (41 DEGs) |

该时间线特征为**脉冲式、自限性**应答: 每次疫苗剂量触发一批转录改变, 随后消退, D365仍有少数基因维持改变。D60的接近回零表明T细胞应答不是持续激活状态。

### 3.3.3 持久转录特征的鉴定

**数据与图表**:
- 免疫基因点阵图 (人工筛选22基因): `results/rnaseq/FigB_ImmuneGenes_dotplot.pdf`
- 数据驱动点阵图 (Top25, 显著≥2时间点): `results/rnaseq/FigB_DataDriven_dotplot.pdf`
- DEG热图 (Top200): `results/rnaseq/FigB_Heatmap_DEGs_timepoints.pdf`
- 五基因纵向轨迹: `results/chapter3_figures/five_gene_trajectory.png`

**跨时间点持续改变的基因**:

| 基因 | D14 | D60 | D74 | D365 | 生物学解释 |
|------|-----|-----|-----|------|-----------|
| **ZEB2** | +3.2 | +2.5 | +3.6 | +3.0 | T细胞效应-记忆分化的转录调控因子 |
| **CTLA4** | +1.5 | +1.4 | +1.6 | +1.0 | 激活诱导的免疫检查点, T细胞激活后正常上调 |
| **ICOS** | +1.3 | +0.8 | +1.2 | +0.6 | T细胞共刺激分子 |
| **HAVCR2** | +7.6 | +8.4 | +9.1 | +8.0 | TIM-3, 在急性激活效应T细胞和记忆群体中表达 |

> **关于HAVCR2 (TIM-3)**: 文献中TIM-3常与T细胞耗竭关联, 但近期证据表明其在急性激活的效应T细胞上也高表达, 且当不与TOX持续升高及多受体共抑制时, 可标志功能性效应-记忆群体。本数据中TOX未呈持续升高, 提示HAVCR2的上调更反映效应-记忆编程而非耗竭。

**与活化相关的下调基因**:

| 基因 | D14 | D365 | 生物学解释 |
|------|-----|------|-----------|
| **SPP1** (Osteopontin) | −5.4 | −5.9 | 多效性细胞因子, 参与炎症、T细胞迁移和存活; 持续下调可能反映T细胞活化相关通路的重塑 |
| **IL7R** | −1.2 | −1.0 | IL-7受体, 在效应分化中下调 |

> **关于SPP1**: 该基因功能复杂, 涉及促炎、T细胞存活和Th1极化等多个层面。此处不下功能归因结论, 仅记录其作为RZV应答中一致下调的基因之一。

### 3.3.4 gE-CD4⁺ T细胞中的ISG状态

在RZV疫苗接种后的全部时间点, gE特异性CD4⁺ T细胞中I型IFN相关ISGs保持基线水平:

| 基因 | HZ全血 LFC | RZV D14 | RZV D74 | RZV D365 |
|------|-----------|---------|---------|----------|
| ISG15 | +1.57 ** | +0.1 ns | −0.4 ns | −0.1 ns |
| RSAD2 | +1.55 *  | −0.3 ns | −0.2 ns | −0.2 ns |
| MX1   | +0.82    | +0.3 ns | −0.2 ns | −0.1 ns |
| IFI44L| +1.82 *  | −0.3 ns | −0.2 ns | −0.2 ns |

**该结果的意义**:

1. 这是**预期的、与细胞生物学一致的观察**——CD4⁺ T细胞非I型IFN主要产生细胞, 其不表达ISGs是细胞类型的正常属性。
2. 该数据表明RZV诱导的gE-CD4⁺ T细胞应答**不伴随ISG共激活**——即T细胞未被炎症环境"旁激活"。
3. **本章不据此做出"RZV在全血或全身水平不激活IFN通路"的推断**。关于RZV接种者系统性固有免疫应答的评价需全血或PBMC层面数据。

### 3.3.5 小节

| 结论 | 内容 |
|------|------|
| **RZV gE-CD4⁺ T细胞特征** | 脉冲式、自限性转录重编程, D365维持效应-记忆特征 |
| **Protection Reference** | ZEB2, CTLA4, HAVCR2 — T细胞持久编程的代表性特征 |
| **用途** | 作为**适应性免疫"质量"的参照** — 候选疫苗应能在抗原特异性T细胞中建立类似的持久分化程序 |

---

## 3.4 二维评价框架的建立

### 3.4.1 为什么不做直接比较

HZ全血Bulk RNA-seq和RZV gE-CD4⁺ T细胞scRNA-seq在以下维度不同:

| 维度 | GSE242252 | GSE249632 |
|------|-----------|-----------|
| 细胞来源 | 全血 (所有免疫细胞) | 纯化gE四聚体⁺ CD4⁺ T细胞 |
| 测序平台 | 3′ mRNA-seq (Bulk) | SMART-Seq v4 (scRNA-seq) |
| 生物状态 | 疾病 (病毒再激活) | 预防 (重组蛋白+佐剂) |
| 信息层级 | 系统全景 | 抗原特异性T细胞内在程序 |

这些差异使得直接比较 (如"RZV与HZ哪个更强"或"两者是否本质不同") **不具有可解释性**。然而, 正是这些差异赋予了两个数据集**互补的信息优势**: HZ全血揭示疾病中系统性激活的全部通路 (包括T细胞信号被稀释的弱点), RZV T细胞揭示成功疫苗如何在抗原特异性T细胞中建立持久编程 (包括全血中不可见的分化细节)。

### 3.4.2 两个适应证的独立参照标准

| | 标准1: 固有免疫激活模式 | 标准2: 适应性免疫编程质量 |
|---|---|---|
| **参照数据集** | GSE242252 (HZ全血) | GSE249632 (RZV gE-CD4⁺ T) |
| **关键发现** | HZ激活I型IFN程序、补体、B细胞抗体应答 | RZV诱导ZEB2/CTLA4/HAVCR2持久上调和脉冲式自限性应答 |
| **评价核心基因** | ISG15, RSAD2, IFI44L, IFI27, SERPING1 | ZEB2, CTLA4, HAVCR2 |
| **候选疫苗期望** | 不触发HZ水平的全身性IFN相关炎症 | 建立类似的持久T细胞效应-记忆编程 |
| **评价章节** | 第四章 (固有免疫评价) | 第五、六章 (T细胞应答评价) |
| **评价方法** | 全血/PBMC bulk或固有免疫细胞检测 | 抗原特异性T细胞转录组/细胞因子/功能测定 |

### 3.4.3 评价框架的科学基础

候选疫苗GE282+GB705的目标是: 通过佐剂启动固有免疫以驱动适应性应答, 同时避免触发RIG-I/MDA5→IRF→ISG级联 (这是病毒核酸感知通路, 而非重组蛋白亚单位疫苗应激活的通路)。TLR4-MyD88→NF-κB/AP-1通路足以提供DC成熟和T细胞启动所需的信号, 且不引发全身性IFN炎症。因此:

- **标准1** 实质是验证候选疫苗是否确实通过TLR4而非RIG-I通路激活固有免疫;
- **标准2** 实质是验证由此驱动的适应性免疫是否达到了RZV级别的T细胞编程质量;
- 两个标准独立评价、联合决策——正如两个数据集各自独立提供信息、联合定义评价空间。

---

## 3.5 本章结论

**结论1**: HZ急性期全血转录组揭示了**I型IFN相关系统性抗病毒程序**的全面激活, 涉及ISGs、补体活化、B细胞抗体产生和免疫细胞增殖。ISG信号主要来自单核细胞等天然免疫群体 (HRA008316独立验证)。该图谱定义了候选疫苗应在固有免疫层面**避免**的炎症模式 (Disease Reference)。

**结论2**: RZV疫苗诱导的gE特异性CD4⁺ T细胞呈现**脉冲式、自限性转录重编程**, 并以持续至D365的效应-记忆特征 (ZEB2, CTLA4, HAVCR2) 区别于短期激活。该模式定义了候选疫苗应在适应性免疫层面**获得**的T细胞编程质量 (Protection Reference)。

**结论3**: HZ全血和RZV gE-CD4⁺ T细胞两个数据集分别锚定疫苗免疫评价的两个独立维度——固有免疫激活的"广度"和适应性免疫的"质量"。它们不应被直接比较, 而应联合定义二维评价空间。候选疫苗GE282+GB705将在以下章节沿这两个维度分别评价: 第四章评估其是否避免了Disease Reference定义的广泛炎症模式; 第五/六章评估其是否达到了Protection Reference定义的T细胞编程质量。

---

## 附录: 输出文件索引

### A. GSE242252 (HZ全血Bulk) 输出

| 文件 | 内容 |
|------|------|
| `results/GSE242252/DE_HZ_acute_vs_convalescent.csv` | 全部基因DESeq2结果 (含baseMean/LFC/pvalue/padj/symbol) |
| `results/GSE242252/Significant_DEGs_table.csv` | 显著DEG表 (p<0.05, |LFC|>0.58) |
| `results/GSE242252/GO_BP_enrichment.csv` | GO BP富集结果 |
| `results/GSE242252/KEGG_enrichment.csv` | KEGG通路富集结果 |
| `results/rnaseq/FigA_Volcano_HZ.pdf` | 火山图 (Top20标注) |
| `results/rnaseq/FigA_PCA_HZ.pdf` | PCA图 |
| `results/rnaseq/FigA_Heatmap_DEGs.pdf` | 显著DEG热图 |
| `results/rnaseq/FigA_GO_dotplot.pdf` | GO气泡图 |
| `results/rnaseq/FigA_KEGG_dotplot.pdf` | KEGG气泡图 |

### B. GSE249632 (RZV gE-CD4⁺ T) 输出

| 文件 | 内容 |
|------|------|
| `results/GSE249632/DE_D{14,60,74,365}_vs_D0.csv` | 各时间点差异表达 |
| `results/GSE249632/logFC_matrix_all_timepoints.csv` | 4时间点logFC合并矩阵 |
| `results/GSE249632/Curated_ImmuneGenes_LFC.csv` | 人工筛选22免疫基因LFC表 |
| `results/GSE249632/DataDriven_Top25_DEGs.csv` | 数据驱动Top25基因表 |
| `results/rnaseq/FigB_DEG_timeline_barchart.pdf` | DEG时间线柱状图 |
| `results/rnaseq/FigB_Volcano_D{14,60,74,365}.pdf` | 各时间点火山图 (4张) |
| `results/rnaseq/FigB_Upset_DEGs.pdf` | DEG重叠Upset图 |
| `results/rnaseq/FigB_Heatmap_DEGs_timepoints.pdf` | Top200 DEG热图 |
| `results/rnaseq/FigB_ImmuneGenes_dotplot.pdf` | 免疫基因点阵图 (按功能分类) |
| `results/rnaseq/FigB_DataDriven_dotplot.pdf` | 数据驱动点阵图 |

### C. 跨数据集与综合

| 文件 | 内容 |
|------|------|
| `results/chapter3_figures/five_gene_trajectory.png` | 五基因纵向轨迹 (ZEB2/CTLA4/ICOS/HAVCR2/ISG15) |
| `results/chapter3_figures/quadrant_narrative.png` | 象限散点图 |
| `results/chapter3_figures/bulk_vs_scRNA_comparison.png` | Bulk vs scRNA ISG来源对比 |
| `results/Chapter3_Report.html` | 第三章HTML报告 |
| `results/Chapter3_Presentation.pptx` | 第三章PPT (9页) |

### D. 分析脚本

| 文件 | 内容 |
|------|------|
| `scripts/chapter3_analysis.R` | 完整R分析脚本 (Part A: HZ DESeq2, Part B: RZV limma, Part C: 基因提取) |
| `scripts/run_all.sh` | 一键全流程脚本 |
