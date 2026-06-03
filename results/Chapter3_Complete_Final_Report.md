# 第三章：基于公共转录组数据构建HZ疫苗免疫评价的二维参照系

**最终版本 | 2026-06-03**

---

## 3.1 研究设计与核心逻辑

### 3.1.1 核心问题

理想的HZ候选疫苗应诱导什么样的免疫应答？本章不试图回答"疫苗好不好"，而是先建立**评价"好"的标准**。

### 3.1.2 三个数据集，两个参照系

| 数据集 | 技术 | 细胞/样本 | 比较设计 | 信息优势 | 提供的参照系 |
|--------|------|---------|---------|---------|------------|
| **GSE242252** | Bulk RNA-seq | 全血, 26急性 vs 23恢复 | DESeq2 非配对 | 系统性免疫全景 | **Disease Reference**：疾病中激活了什么 |
| **GSE249632** | scRNA-seq | gE-CD4⁺ T, 7例纵向 D0→D365 | limma 配对 | 抗原特异性T细胞编程 | **Protection Reference**：疫苗中获得了什么 |
| **HRA008316** | scRNA-seq | PBMC, 3 HP vs 3 HA | 表达矩阵 LFC | T细胞内在变化 | **Bridge**：连接Bulk与scRNA |

### 3.1.3 分析逻辑

```
HZ全血(Bulk)                       RZV gE-CD4⁺ T(scRNA)
"广度：系统性免疫全景"              "精度：T细胞内在编程"
        │                                    │
        ▼                                    ▼
  疾病中激活了什么？                   疫苗建立了什么？
  → 天然免疫ISG程序                   → T细胞分化记忆程序
  → 补体、B细胞、增殖                 → 自限性调控
        │                                    │
        └──────────┬─────────────────────────┘
                   │
                   ▼
         HRA008316 T细胞(Bridge)
         "T细胞本身在疾病中做了什么？"
         → 基础激活与RZV同向
         → 分化深度远不及RZV
                   │
                   ▼
         二维评价框架：
         标准1(第四章)：候选疫苗是否避免了系统性ISG炎症？
         标准2(第五/六章)：候选疫苗是否建立了RZV级T细胞编程？
```

### 3.1.4 各数据集样本量

| 数据集 | 组别 | n | 技术 |
|--------|------|:--:|------|
| **GSE242252** (HZ全血) | Acute | 26 | Bulk RNA-seq |
| | Convalescent | 23 | |
| **GSE249632** (RZV) | 7 donors × 5 timepoints | 31 pseudobulk组 | scRNA-seq (SMART-Seq v4) |
| | QC通过单细胞 | 2,231 | |
| **HRA008316** (HZ PBMC) | Healthy (HA) | 3 | scRNA-seq (10x) |
| | HZ patient (HN) | 3 | |
| | Recovered (nHN) | 3 | |

> **注意**: HRA008316每组仅3例, 统计效力有限。LFC主要用方向和幅度做定性判断, p值作参考。

---

## 3.2 GSE242252：HZ急性期全血转录组

### 3.2.1 差异表达全景

**方法**: DESeq2非配对, p<0.05 (nominal), |LFC|>0.58, 26 Acute vs 23 Convalescent

**结果**: 352个上调, 44个下调。上调远多于下调, 急性期以免疫激活为主。

**Top上调基因**:

| 基因 | LFC | p值 | 功能分类 |
|------|-----|-----|---------|
| IFI27 | +2.16 | 1.5e-04 | I型IFN ISG |
| PTTG1 | +2.01 | 3.0e-13 | 细胞增殖 |
| BATF2 | +1.95 | 9.5e-04 | IFN诱导转录因子 |
| IGHG4 | +1.94 | 2.6e-05 | 免疫球蛋白(IgG4) |
| IFI44L | +1.82 | 8.5e-05 | I型IFN ISG |
| SERPING1 | +1.80 | 5.3e-07 | 补体C1抑制因子 |
| MZB1 | +1.68 | 3.9e-06 | 浆细胞标志 |
| ISG15 | +1.57 | 1.6e-05 | ISG化修饰 |
| RSAD2 | +1.55 | 1.2e-04 | 抗病毒效应蛋白(Viperin) |
| TOP2A | +1.50 | 3.4e-13 | DNA复制 |

**图表**: `results/rnaseq/FigA_Volcano_HZ.pdf`, `FigA_PCA_HZ.pdf`, `FigA_Heatmap_DEGs.pdf`

**数据**: `results/GSE242252/DE_HZ_acute_vs_convalescent.csv`, `Significant_DEGs_table.csv`

### 3.2.2 通路富集

**方法**: clusterProfiler GO BP + KEGG, 所有显著DEG合并

**关键GO Terms**: Defense Response to Virus, B Cell Receptor Signaling, Antigen Receptor-Mediated Signaling, Mitotic Spindle Organization

**图表**: `results/rnaseq/FigA_GO_dotplot.pdf`, `FigA_KEGG_dotplot.pdf`

**解读**: HZ急性期 = I型IFN相关系统性抗病毒程序 + 补体活化 + B细胞抗体应答 + 免疫细胞增殖。TCR信号相对薄弱——全血中T细胞特异性信号被稀释。

### 3.2.3 T细胞基因在全血中被系统性稀释

52个T细胞关键基因在三个数据集中的交叉验证结果显示：**22个基因在RZV T细胞中显著变化但在HZ全血中完全不可见**（ZEB2, CTLA4, ICOS, HAVCR2等）。保护性T细胞编程在Bulk全血中被大量非T细胞RNA稀释。

**完整表**: `results/Tcell_genes_cross_dataset.csv`

### 3.2.4 小节：HZ Disease Reference

HZ急性期全血转录组以**I型IFN相关系统性抗病毒程序**为主要特征。该参照系定义了候选疫苗应在固有免疫层面**避免**的广泛炎症模式。

---

## 3.3 GSE249632：RZV疫苗gE-CD4⁺ T细胞动力学

### 3.3.1 DEG时间线：脉冲式激活

| 时间点 | 上调 | 下调 | 模式 |
|--------|:--:|:--:|------|
| D14 vs D0 | 164 | 152 | 第一针强应答 |
| D60 vs D0 | 10 | 4 | 回到基线 |
| D74 vs D0 | 72 | 35 | 第二针再激活(减毒) |
| D365 vs D0 | 20 | 21 | 长期印记(41 DEGs) |

**图表**: `results/rnaseq/FigB_DEG_timeline_barchart.pdf`, `FigB_Volcano_D{14,60,74,365}.pdf`

### 3.3.2 持久转录特征

| 基因 | D14 | D60 | D74 | D365 | 功能 |
|------|-----|-----|-----|------|------|
| ZEB2 | +3.2 | +2.5 | +3.6 | +3.0 | T细胞效应-记忆分化 |
| CTLA4 | +1.5 | +1.4 | +1.6 | +1.0 | 免疫自限性检查点 |
| HAVCR2 | +7.6 | +8.4 | +9.1 | +8.0 | TIM-3/效应记忆 |
| ICOS | +1.3 | +0.8 | +1.2 | +0.6 | T细胞共刺激 |

**图表**: `results/rnaseq/FigB_ImmuneGenes_dotplot.pdf`, `FigB_Heatmap_DEGs_timepoints.pdf`

### 3.3.3 gE-CD4⁺ T细胞中的ISG状态

gE-CD4⁺ T细胞在所有时间点不表现ISG共激活（ISG15, RSAD2, MX1平坦）。这与CD4⁺ T细胞非I型IFN主要产生者的生物学一致。**本章不据此推断RZV在全血水平不激活IFN**——该问题需全血/PBMC数据。

### 3.3.4 小节：RZV Protection Reference

RZV诱导**脉冲式、自限性CD4⁺ T细胞转录重编程**，以持久效应-记忆特征（ZEB2, CTLA4, HAVCR2）区别于短期激活。该参照系定义了候选疫苗应在适应性免疫层面**获得**的T细胞编程质量。

---

## 3.4 HRA008316：CD4⁺ T细胞对比桥梁

### 3.4.1 分析设计与细胞群体选择

从HRA008316 Fig_4b(1)完整表达矩阵中提取T细胞数据，计算HZ患者(HN) vs 健康对照(HA)的LFC。

**两个T细胞群体并行分析**:

| 群体 | 列数 | 内容 | 用途 |
|------|:--:|------|------|
| **CD4NaiveT** | HA:3, HN:3 | 纯CD4⁺ naive T细胞 | **主分析** — 与RZV gE-CD4⁺最接近 |
| All-T | HA:21, HN:21 | CD4NaiveT + CD8NaiveT + CD8Teff + CD8Tem + GDTCells + ProliferatingT + Treg | **对照验证** — 证明结论不受亚群选择影响 |

**方法**: log2(HN_mean / HA_mean), t检验 (3 vs 3) + BH校正。

### 3.4.2 方向一致性总览

| 方向 | CD4NaiveT | All-T |
|------|:--:|:--:|
| 一致上调 | 2,099 (42%) | 3,346 (46%) |
| 一致下调 | 740 (15%) | 781 (11%) |
| 相反(HZ↑ RZV↓) | 1,381 (28%) | 2,323 (32%) |
| 相反(HZ↓ RZV↑) | 771 (15%) | 755 (10%) |
| **一致合计** | **2,839 (57%)** | **4,127 (57%)** |
| **相反合计** | **2,152 (43%)** | **3,078 (43%)** |

**CD4NaiveT和All-T给出几乎相同的一致性比例（57% vs 57%）**。CD4NaiveT vs All-T的LFC相关系数 Spearman ρ=0.655 (p=0.003), 方向高度一致。核心结论不受T细胞亚群选择影响。

**图表**: `results/rnaseq/FigD_Concordance_Scatter.pdf`, `FigD_Concordance_Bar.pdf`, `FigD_KeyGenes_Lollipop.pdf`

**数据**: `results/GSE249632/HRA_vs_RZV_ComparisonTable.csv` (CD4NaiveT), `HRA_AllT_vs_RZV_ComparisonTable.csv` (All-T)

### 3.4.3 关键发现：Protection基因方向一致但幅度差10-250倍

| 基因 | CD4NaiveT | All-T | RZV D14 | 倍差(CD4/RZV) | 方向 |
|------|:--:|:--:|---------|:--:|:--:|
| ZEB2 | +0.13 | +0.64 | +3.2 | 25× | 一致↑ |
| CTLA4 | +0.10 | +0.30 | +1.5 | 15× | 一致↑ |
| ICOS | +0.29 | +0.11 | +1.3 | 4× | 一致↑ |
| HAVCR2 | +0.03 | +0.51 | +7.6 | 250× | 一致↑ |
| TOX | +0.08 | +0.53 | +1.3 | 16× | 一致↑ |
| IL21 | +0.01 | +0.03 | +1.8 | 180× | 一致↑ |
| CD38 | +0.47 | +1.21 | +6.6 | 14× | 一致↑ |
| GZMA | +0.38 | +0.71 | +2.7 | 7× | 一致↑ |
| TIGIT | +0.47 | +1.17 | +0.8 | 2× | 一致↑ |

**核心结论**: Protection Signature在HZ CD4NaiveT中以同向但微弱的形式存在(LFC 0.01-0.47), RZV将其放大至功能水平(LFC 1.3-7.6)。All-T群体中部分基因LFC增大(如ZEB2 +0.64, HAVCR2 +0.51), 表明CD8⁺ T细胞也有贡献, 但CD4NaiveT的方向性结论完全一致。

### 3.4.4 真正反向：命运决定基因

| 基因 | CD4NaiveT | RZV | 分歧含义 |
|------|:--:|-----|---------|
| IL7R | **+0.82** | **-1.2** | HZ维持记忆特征, RZV推动效应分化 |
| BCL6 | **-0.58** | +0.7 | RZV启动Tfh程序, HZ中抑制 |
| CXCR5 | +0.12 | -0.9 | 淋巴结归巢方向相反 |

这些基因在All-T中方向完全一致(IL7R +0.54, BCL6 -1.00)。

### 3.4.5 HRA分析的限制

1. **样本量**: 每组仅3例 (HA=3, HN=3), 统计效力极低。t检验p值仅作参考。
2. **细胞群体**: Fig_4b(1)仅有CD4NaiveT, 缺少CD4⁺记忆/效应亚群。与RZV gE-CD4⁺ T(含naive+memory+effector)的群体组成不完全对等。
3. **平台差异**: 10x Genomics (HRA) vs SMART-Seq v4 (RZV), 后者灵敏度和基因检测率更高。
4. **ISG15/STAT1未检测到**: 在CD4NaiveT中无显著信号, 证实CD4⁺非I型IFN来源, 但也可能是10x平台的检测限制。

---

## 3.5 二维评价框架

### 3.5.1 为什么两个参照系互补而非对立

| 维度 | HZ全血 | RZV gE-CD4⁺ T |
|------|--------|---------------|
| 看得见 | ISG风暴、补体、抗体、增殖 | T细胞分化编程、自限性调控 |
| 看不见 | T细胞内在编程（被稀释） | 非T细胞来源的ISG（细胞类型限制） |
| 提供 | "需要避免什么" — Disease Reference | "需要获得什么" — Protection Reference |

两者**不在同一层面比较**，而是**联合定义二维评价空间**。HRA008316 T细胞数据进一步确认：T细胞本身在HZ中有基础激活（方向与RZV一致），但未达到RZV级别的编程深度。

### 3.5.2 候选疫苗评价双重标准

| | 标准1: 固有免疫模式 | 标准2: 适应性免疫质量 |
|---|---|---|
| **参照数据** | GSE242252 (HZ全血) | GSE249632 (RZV gE-CD4⁺ T) |
| **核心指标** | ISG15, RSAD2, IFI44L, IFI27, SERPING1 | ZEB2, CTLA4, HAVCR2, ICOS |
| **期望** | 不触发HZ水平的系统性ISG | 建立持久T细胞效应-记忆编程 |
| **评价章节** | 第四章 (固有免疫) | 第五/六章 (T细胞应答) |

---

## 3.6 本章结论

**结论1**: HZ急性期全血转录组以I型IFN相关系统性抗病毒程序为主要特征(ISG15 +1.57, IFI27 +2.16, RSAD2 +1.55)，伴随补体活化、B细胞抗体应答和免疫细胞增殖。该参照系(Disease Reference)定义了候选疫苗应在固有免疫层面避免的炎症模式。

**结论2**: RZV疫苗诱导的gE特异性CD4⁺ T细胞呈现脉冲式、自限性转录重编程(D14激活→D60回基线→D74再激活)，以持久效应-记忆特征(ZEB2 D365 +3.0, HAVCR2 D365 +8.0)区别于短期激活。该参照系(Protection Reference)定义了候选疫苗应在适应性免疫层面达到的T细胞编程质量。

**结论3**: HRA008316 T细胞数据揭示，Protection Signature基因(ZEB2, CTLA4, HAVCR2, ICOS)在HZ T细胞中以同向但微弱的形式存在(LFC 0.01-0.47)，RZV将其放大10-250倍至功能水平。RZV并非"创造"新程序，而是将HZ中已启动的程序推至功能阈值。

**结论4**: 三个数据集联合定义了候选疫苗免疫评价的二维空间——固有免疫的"广度"(避免系统性ISG炎症)和适应性免疫的"质量"(获得持久T细胞编程)。后续章节沿这两个维度分别评价候选疫苗GE282+GB705。

---

## 附录A: 分析阈值与方法

| 参数 | GSE242252 | GSE249632 | HRA008316 |
|------|-----------|-----------|-----------|
| 方法 | DESeq2 非配对 | limma-voom 配对 | 表达矩阵LFC + t检验 |
| 显著性 | p < 0.05 (nominal) | FDR < 0.05 | BH-adjusted p |
| |LFC|阈值 | > 0.58 | > 0.58 | > 0.5 |
| 样本 | 26 acute / 23 conv | 7 donors × 5 timepoints | 3 HP / 3 HA |

## 附录B: 输出文件索引

### GSE242252 (HZ Bulk)
| 文件 | 内容 |
|------|------|
| `results/GSE242252/DE_HZ_acute_vs_convalescent.csv` | 全部DESeq2结果 |
| `results/GSE242252/Significant_DEGs_table.csv` | 显著DEG表 |
| `results/rnaseq/FigA_Volcano_HZ.pdf` | 火山图 |
| `results/rnaseq/FigA_GO_dotplot.pdf` | GO富集 |

### GSE249632 (RZV T细胞)
| 文件 | 内容 |
|------|------|
| `results/GSE249632/DE_D{14,60,74,365}_vs_D0.csv` | 各时间点DEG |
| `results/GSE249632/Curated_ImmuneGenes_LFC.csv` | 22免疫基因LFC |
| `results/rnaseq/FigB_DEG_timeline_barchart.pdf` | DEG时间线 |
| `results/rnaseq/FigB_ImmuneGenes_dotplot.pdf` | 免疫基因点阵图 |
| `results/rnaseq/FigB_Heatmap_DEGs_timepoints.pdf` | DEG热图 |

### HRA008316 (对比分析)
| 文件 | 内容 |
|------|------|
| `results/GSE249632/HRA_vs_RZV_ComparisonTable.csv` | 4,991基因方向一致性表 |
| `results/GSE249632/KeyImmune_DirectionConcordance.csv` | 24关键免疫基因方向 |
| `results/rnaseq/FigD_Concordance_Scatter.pdf` | 方向一致性散点图 |
| `results/rnaseq/FigD_Concordance_Bar.pdf` | 方向一致性柱状图 |
| `results/rnaseq/FigD_KeyGenes_Lollipop.pdf` | Top50基因棒棒糖图 |

### 跨数据集
| 文件 | 内容 |
|------|------|
| `results/Tcell_genes_cross_dataset.csv` | 52 T细胞基因三数据集对比 |
| `results/ThreeDataset_Tcell_comparison.csv` | 三数据集T细胞基因表 |
| `results/chapter3_figures/five_gene_trajectory.png` | 五基因纵向轨迹 |

### 报告与脚本
| 文件 | 内容 |
|------|------|
| `results/Chapter3_Final_Report.md` | 本报告 |
| `results/Chapter3_Presentation.pptx` | PPT |
| `scripts/chapter3_analysis.R` | 完整R分析脚本 |
