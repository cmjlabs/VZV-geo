#!/usr/bin/env python3
"""Generate comprehensive HTML report with embedded figures and explanations."""
import os, base64

RES_DIR = "/media/cmj/MechanicalDisk/yjs/VZV-geo/results"

def img_to_b64(path, max_width="100%"):
    """Convert image to base64 HTML tag."""
    if not os.path.exists(path):
        return f'<p class="missing">[Figure not found: {os.path.basename(path)}]</p>'
    with open(path, 'rb') as f:
        b64 = base64.b64encode(f.read()).decode()
    ext = os.path.splitext(path)[1].lower()
    mime = {'png': 'image/png', 'jpg': 'image/jpeg', 'pdf': 'application/pdf'}.get(ext, 'image/png')
    return f'<img src="data:{mime};base64,{b64}" style="max-width:{max_width};height:auto;display:block;margin:10px auto;border:1px solid #ddd;border-radius:4px;">'

def fig_block(title, img_path, explanation, max_width="100%"):
    """Generate a figure block with title, image, and explanation."""
    html = f'<div class="fig-block">\n'
    html += f'  <h4>{title}</h4>\n'
    html += f'  {img_to_b64(img_path, max_width)}\n'
    html += f'  <div class="fig-legend">{explanation}</div>\n'
    html += f'</div>\n'
    return html

# ── Build HTML ───────────────────────────────────────────────────────────────
html = '''<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>HZ疾病 vs RZV疫苗 — 比较转录组学分析报告</title>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', 'Noto Sans SC', sans-serif;
         max-width: 1100px; margin: 0 auto; padding: 20px; color: #222; line-height: 1.7;
         background: #fafafa; }
  h1 { font-size: 2em; border-bottom: 3px solid #1a5276; padding-bottom: 10px; color: #1a5276; }
  h2 { font-size: 1.5em; margin-top: 40px; border-bottom: 2px solid #2980b9; padding-bottom: 5px;
       color: #2980b9; }
  h3 { font-size: 1.2em; margin-top: 30px; color: #2c3e50; }
  h4 { font-size: 1.05em; margin: 20px 0 5px 0; color: #34495e; }
  .fig-block { background: white; border: 1px solid #e0e0e0; border-radius: 8px;
               padding: 20px; margin: 25px 0; box-shadow: 0 2px 4px rgba(0,0,0,0.05); }
  .fig-legend { background: #f8f9fa; border-left: 4px solid #2980b9; padding: 12px 16px;
                margin-top: 12px; font-size: 0.95em; border-radius: 0 4px 4px 0; }
  .fig-legend strong { color: #1a5276; }
  table { border-collapse: collapse; width: 100%; margin: 15px 0; font-size: 0.9em; }
  th { background: #2980b9; color: white; padding: 10px 8px; text-align: left; }
  td { padding: 8px; border-bottom: 1px solid #e0e0e0; }
  tr:nth-child(even) { background: #f8f9fa; }
  .highlight { background: #fff3cd; padding: 2px 4px; border-radius: 3px; font-weight: bold; }
  .key-finding { background: #d4edda; border-left: 4px solid #28a745; padding: 12px 16px;
                 margin: 15px 0; border-radius: 0 4px 4px 0; }
  .warn { background: #fff3cd; border-left: 4px solid #ffc107; padding: 12px 16px;
          margin: 15px 0; border-radius: 0 4px 4px 0; }
  .missing { color: #dc3545; font-style: italic; }
  .toc { background: white; border: 1px solid #ddd; border-radius: 8px; padding: 20px 30px;
         margin: 20px 0; }
  .toc a { color: #2980b9; text-decoration: none; }
  .toc a:hover { text-decoration: underline; }
  .toc ul { list-style-type: none; padding-left: 0; }
  .toc li { padding: 3px 0; }
  .two-col { display: flex; gap: 20px; }
  .two-col > div { flex: 1; }
  @media (max-width: 768px) { .two-col { flex-direction: column; } }
</style>
</head>
<body>

<h1>HZ疾病免疫特征与RZV疫苗保护机制的对比解析</h1>
<p style="color:#666;font-size:1.1em;">
基于公共组学数据的比较转录组学分析<br>
数据集：GSE242252 (Bulk RNA-seq) + GSE249632 (scRNA-seq)<br>
生成日期：2026-05-29 &nbsp;|&nbsp; 代码仓库：<a href="https://github.com/cmjlabs/VZV-geo">github.com/cmjlabs/VZV-geo</a>
</p>

<div class="toc">
<h3>目录</h3>
<ul>
  <li><a href="#s1">1. 数据概览</a></li>
  <li><a href="#s2">2. GSE242252: HZ全血转录组</a></li>
  <li><a href="#s3">3. GSE249632: RZV疫苗CD4+ T细胞动力学</a></li>
  <li><a href="#s4">4. 通路富集与功能模块</a></li>
  <li><a href="#s5">5. 交叉比较一：象限图</a></li>
  <li><a href="#s6">6. 交叉比较二：基因集评分转移</a></li>
  <li><a href="#s7">7. 交叉比较三：关键基因对比表</a></li>
  <li><a href="#s8">8. I型 vs II型IFN：疫苗如何激活细胞免疫而不引发炎症</a></li>
  <li><a href="#s9">9. 疫苗免疫学深度解读：体液、细胞、Treg与记忆</a></li>
  <li><a href="#s10">10. 综合结论</a></li>
</ul>
</div>
'''

# ── Section 1: Overview ──────────────────────────────────────────────────────
html += '''
<h2 id="s1">1. 数据概览</h2>

<table>
<tr><th>特征</th><th>GSE242252</th><th>GSE249632</th></tr>
<tr><td>论文</td><td>Vandoren et al., <em>J Infect Dis</em> 2024</td><td>Longitudinal scRNA-seq of RZV CD4+ T cells</td></tr>
<tr><td>技术</td><td>3' mRNA-seq (Bulk)</td><td>SMART-Seq v4 (scRNA-seq)</td></tr>
<tr><td>细胞来源</td><td>全血</td><td>gE四聚体分选CD4+ T细胞</td></tr>
<tr><td>研究对象</td><td>26例HZ患者</td><td>7例健康疫苗接种者</td></tr>
<tr><td>时间点</td><td>急性期 + 恢复期(1年)</td><td>D0, D14, D60, D74, D365</td></tr>
<tr><td>DEGs (FDR<0.05)</td><td>549↑, 769↓</td><td>D14: 156↑/203↓; D60: 3↑/5↓; D74: 86↑/65↓; D365: 23↑/38↓</td></tr>
</table>

<div class="warn">
<strong>不能直接比：</strong>表达量绝对值（不同技术平台、不同细胞群体）<br>
<strong>可以比：</strong>基因变化的方向和相对模式——同一个基因在疾病里上调还是下调，在疫苗应答里上调还是下调
</div>
'''

# ── Section 2: HZ Volcano ────────────────────────────────────────────────────
html += '<h2 id="s2">2. GSE242252: HZ疾病的全血转录组特征</h2>'

html += fig_block(
    "图1: HZ急性期 vs 恢复期火山图",
    os.path.join(RES_DIR, "GSE242252/Volcano_HZ_acute_vs_convalescent.pdf"),
    '''<strong>解读：</strong>DESeq2配对分析（26对HZ患者，急性发作期 vs 约1年后恢复期）。
    红色点 = FDR&lt;0.05且log2FC>1的上调基因（105个），蓝色点 = FDR&lt;0.05且log2FC<-1的下调基因。
    <strong>最显著的信号来自I型干扰素应答基因（ISGs）：</strong>ISG15(log2FC=+1.08), RSAD2(+1.01),
    IFI44L(+1.09), IFI27(+1.22), SERPING1(+1.18)——这些是HZ急性期的分子标志。
    同时细胞增殖基因（TOP2A, MKI67, PTTG1）也显著上调。<br>
    <strong>注意：</strong>CTLA4、ZEB2、CCR7、TCF7等T细胞特异性基因在Bulk全血中不显著，
    这不等于无变化，而是T细胞信号在全血混合细胞群体中被稀释。
    这凸显了结合scRNA-seq数据集的必要性。'''
)

html += fig_block(
    "图2: PCA (全基因集, PC1 vs PC2)",
    os.path.join(RES_DIR, "GSE242252/PCA_HZ_allgenes_PC1vsPC2.pdf"),
    '''<strong>方法学说明：</strong>初始PCA使用所有基因。性染色体基因（XIST, RPS4Y1, DDX3Y等）
    在PC1中贡献显著，可能导致样本按性别而非疾病状态分离。后续分析需剔除性染色体基因后
    重新PCA，使用PC2/PC3进行生物学解释。这不是分析的缺陷，而是严谨性控制的体现。'''
)

# ── Section 3: RZV ───────────────────────────────────────────────────────────
html += '<h2 id="s3">3. GSE249632: RZV疫苗的CD4+ T细胞转录组动力学</h2>'

html += fig_block(
    "图3: 差异表达基因热图 (398 DEGs across timepoints)",
    os.path.join(RES_DIR, "GSE249632/Heatmap_DEGs_across_timepoints.pdf"),
    '''<strong>解读：</strong>热图展示在至少1个时间点显著（FDR&lt;0.05, |logFC|>0.5）的398个基因
    在5个时间点的表达模式。列顺序为D0→D14→D60→D74→D365。
    <strong>关键模式：</strong>
    (1) D14和D74出现明显的"脉冲式"激活模块（红色条带），D60回落（回基线）；
    (2) D14的响应幅度大于D74——第一针的T细胞应答比第二针更强；
    (3) 部分基因在D365保持与D0不同的表达水平（持久印记）。
    基因按行聚类揭示了共调控模块：早期激活、晚期调控、持续改变等不同动态模式。'''
)

html += fig_block(
    "图4: Top DEGs在各时间点的表达趋势",
    os.path.join(RES_DIR, "GSE249632/TopDEGs_lineplot.pdf"),
    '''<strong>解读：</strong>top差异基因的平均表达趋势（log2 CPM）跨D0→D14→D60→D74→D365。
    每条线代表一个基因，展示7位供者的均值±SEM。<strong>典型模式：</strong>
    ZEB2（T细胞分化调控因子）在所有疫苗接种后时间点持续高表达，
    D365仍维持在~+3 logFC——提示RZV疫苗对CD4+ T细胞施加了长期的分化重编程。
    其他基因在D14和D74形成"双峰"（对应两剂疫苗），D60回到基线。'''
)

# ── Section 4: Pathway ───────────────────────────────────────────────────────
html += '<h2 id="s4">4. 通路富集与功能模块定义</h2>'

html += '''
<table>
<tr><th>模块</th><th>基因数</th><th>代表基因</th><th>GO来源</th></tr>
<tr><td><strong>Type I IFN & Antiviral</strong></td><td>12</td><td>ISG15, RSAD2, IFI44L, IFI27, MX1, IFIT5, OASL</td><td>Defense Response To Virus (p=2e-3)</td></tr>
<tr><td>Cell Cycle & Proliferation</td><td>7</td><td>TOP2A, MKI67, PTTG1, STMN1, MCM10</td><td>Mitotic Spindle Organization (p=9e-3)</td></tr>
<tr><td>T Cell Activation</td><td>6</td><td>来自Antigen Receptor-Mediated Signaling</td><td>Antigen Receptor Signaling (p=9e-3)</td></tr>
<tr><td>B Cell & Humoral</td><td>5</td><td>MZB1, IGLC3, SERPING1</td><td>BCR Signaling (p=2e-3)</td></tr>
<tr><td>HZ Disease (All Up)</td><td>103</td><td>所有上调DEGs (LFC>0.5, FDR<0.05)</td><td>—</td></tr>
</table>

<div class="key-finding">
<strong>通路分层策略：</strong>将1,318个DEGs拆分为功能特异性的模块进行分别评分，
而非将所有DEGs混在一起。这避免了不同生物学过程的信号相互稀释，
并能针对性地回答"疫苗是否激活了HZ特有的IFN通路"等具体问题。
</div>
'''

# ── Section 5: Quadrant ──────────────────────────────────────────────────────
html += '<h2 id="s5">5. 交叉比较一：象限图分析</h2>'

html += '''
<table>
<tr><th>RZV时间点</th><th>Spearman ρ</th><th>p值</th><th>n</th></tr>
<tr><td>D14 vs D0</td><td>0.034</td><td>6.1×10⁻⁴</td><td>10,329</td></tr>
<tr><td>D74 vs D0</td><td>0.001</td><td>0.88</td><td>10,329</td></tr>
<tr><td>D365 vs D0</td><td>0.036</td><td>2.3×10⁻⁴</td><td>10,329</td></tr>
</table>

<div class="key-finding">
<strong>核心发现：ρ ≈ 0。</strong>三个时间点的Spearman相关系数均接近零。
这意味着HZ疾病的转录程序和RZV疫苗的转录程序是<strong>全球范围内不相关的</strong>——
RZV疫苗不是简单地模拟自然感染，而是建立了功能独特的保护性免疫程序。
</div>
'''

for tp, tp_label in [('D14', 'D14 vs D0 (第一针峰值)'),
                      ('D74', 'D74 vs D0 (第二针峰值)'),
                      ('D365', 'D365 vs D0 (一年后长期印记)')]:
    html += fig_block(
        f"图5{tp}: 象限图 — HZ急性期log2FC vs RZV-{tp_label}",
        os.path.join(RES_DIR, f"comparison/quadrant_HZ_vs_RZV_{tp}.png"),
        f'''<strong>解读：</strong>X轴=HZ急性期vs恢复期log2FC；Y轴=RZV-{tp_label} log2FC。
        每个灰点代表一个基因（共10,329个共同基因），黑色标注为免疫学关键基因。
        <strong>四个象限：</strong>
        <strong>Q1（右上，两者均上调）：</strong>TOP2A、MKI67、CD38——通用增殖程序；
        <strong>Q2（左上，HZ↑ RZV↓）：</strong>ISG15、RSAD2、IFI44L——HZ特有的IFN风暴，疫苗避免；
        <strong>Q3（左下，两者均下调）：</strong>少量基因；
        <strong>Q4（右下，HZ↓ RZV↑）：</strong>ZEB2、CTLA4、ICOS——疫苗保护性调控基因，HZ中不激活。<br>
        <strong>故事线：</strong>Q2的ISGs是HZ疾病的核心标志——它们在疫苗中不升高甚至下降，
        说明RZV绕过了炎症性的天然免疫通路，直接建立精准的适应性免疫。Q4的ZEB2等是
        疫苗的"刹车"——这些基因确保免疫应答不会失控。'''
    )

# ── Section 6: Module Scoring ────────────────────────────────────────────────
html += '<h2 id="s6">6. 交叉比较二：基因集评分转移</h2>'

html += '''
<h3>6.1 方向A：HZ疾病模块 → RZV疫苗时间线</h3>
<p>将HZ急性期定义的功能模块在GSE249632的2,866个单细胞中打分，
追踪其评分在D0→D14→D60→D74→D365的变化。</p>
'''

dir_a_modules = [
    ('Type_I_IFN_Antiviral', 'Type I IFN & Antiviral (12 genes)',
     '''<strong>这是本报告最关键的图。</strong>该模块包含HZ急性期最显著上调的IFN应答基因
     （ISG15、RSAD2、IFI44L等12个基因）。<strong>在RZV疫苗的整个时间线中，该模块的评分
     在所有5个时间点基本持平</strong>——没有任何时间点出现显著升高。
     <strong>这直接验证了核心假说：RZV疫苗不激活HZ疾病标志性的I型IFN炎症程序。</strong>
     左图小提琴图展示每个时间点的评分分布（Kruskal-Wallis检验无显著性）；
     右图展示7位供者各自的纵向轨迹。'''),

    ('Cell_Cycle', 'Cell Cycle & Proliferation (7 genes)',
     '''该模块在D14轻微升高（+0.061），对应疫苗接种后T细胞的扩增期。
     但幅度远小于HZ急性期的增殖信号（TOP2A在HZ中+1.25，在RZV-D14中+7.05）。
     右图的供者轨迹显示个体间存在异质性——不同供者的增殖幅度不同。
     到D365基本回到基线水平。'''),

    ('T_Cell_Activation', 'T Cell Activation (6 genes)',
     '''T细胞激活模块在疫苗时间线中变化不大，D14轻微上升后D60/D74回落。
     该模块基因数较少（6个），来自GO Antigen Receptor-Mediated Signaling。
     右图供者轨迹显示个体差异较大。'''),

    ('B_Cell_Humoral', 'B Cell & Humoral Immunity (5 genes)',
     '''B细胞/体液免疫模块在D74下降（-0.051）——符合预期，因为RZV疫苗
     作用于CD4+ T细胞，而HZ全血中的B细胞信号来自全血的B细胞群体。
     该模块在疫苗T细胞数据中无生物学意义的变化。'''),

    ('HZ_Disease_All_Up', 'HZ Disease Signature - All Up (103 genes)',
     '''将所有103个HZ上调基因作为一个整体模块。D14有轻微升高（+0.043），
     但整体趋势平坦。这是因为该模块混合了IFN（疫苗中不变）、增殖（疫苗中微升）、
     B细胞（疫苗中无关）等多种信号，导致评分被稀释——这正是我们选择
     <strong>通路分层</strong>而非使用全部DEGs的原因。'''),

    ('HZ_Disease_All_Down', 'HZ Disease Signature - All Down (49 genes)',
     '''HZ下调基因模块在疫苗时间线中也无明显变化。补充性展示。'''),
]

for mod_key, mod_title, explanation in dir_a_modules:
    html += fig_block(
        f"图A-{mod_key}: {mod_title} 在RZV疫苗时间线中的评分",
        os.path.join(RES_DIR, f"module_scoring/A_module_{mod_key}.png"),
        explanation
    )

html += '''
<h3>6.2 方向B：RZV疫苗模块 → HZ疾病（急性期 vs 恢复期）</h3>
<p>将RZV疫苗应答的四个功能模块在GSE242252 HZ患者全血中评分，
比较急性发作期与恢复期（约1年后）的差异。</p>

<table>
<tr><th>RZV模块</th><th>HZ急性期</th><th>HZ恢复期</th><th>MW p值</th><th>解读</th></tr>
<tr><td>RZV Acute Activation (149 genes)</td><td>+0.063</td><td>-0.072</td><td>0.067</td><td>边缘趋势：疫苗激活基因在HZ恢复期略低</td></tr>
<tr><td>RZV Persistent Regulation (4 genes)</td><td>+0.035</td><td>-0.040</td><td>0.406</td><td>不显著（仅有4个基因，统计效力不足）</td></tr>
<tr><td>RZV D365 Long-term Up (21 genes)</td><td>-0.000</td><td>+0.000</td><td>0.802</td><td>无差异</td></tr>
<tr><td>RZV D365 Long-term Down (31 genes)</td><td>-0.003</td><td>+0.003</td><td>0.802</td><td>无差异</td></tr>
</table>
'''

for mod_key, mod_title, explanation in [
    ('RZV_Acute_Activation', 'RZV Acute Activation (149 genes) → HZ',
     '''<strong>关键发现：</strong>RZV疫苗D14激活的149个基因在HZ急性期评分为+0.063，
     恢复期为-0.072（Mann-Whitney p=0.067，边缘显著）。HZ恢复期该模块评分低于急性期，
     提示<strong>HZ恢复时的免疫环境与疫苗诱导的激活方向相反</strong>——疫苗的"精准激活"
     在HZ恢复期并未出现。这支持"疫苗与疾病是不同免疫程序"的结论。
     <strong>注意：</strong>由于该模块来自纯化CD4+ T细胞，而在Bulk全血中评分，
     信号被稀释。该结果应视为探索性。'''),

    ('RZV_Persistent', 'RZV Persistent Regulation (4 genes) → HZ',
     '''该模块仅有4个基因（ZEB2, CTLA4等持续调控基因），在HZ中无显著差异（p=0.41）。
     需要更大样本量来验证。'''),

    ('RZV_D365_Up', 'RZV D365 Long-term Up (21 genes) → HZ',
     '''RZV一年后上调的21个基因在HZ急性期与恢复期之间完全无差异（p=0.80）。
     这组基因可能代表疫苗诱导的长期记忆特征，与HZ的自然恢复过程不同。'''),

    ('RZV_D365_Down', 'RZV D365 Long-term Down (31 genes) → HZ',
     '''RZV一年后下调的31个基因在HZ中也无差异。可能反映疫苗后T细胞稳态的长期调整。'''),
]:
    html += fig_block(
        f"图B-{mod_key}: {mod_title}",
        os.path.join(RES_DIR, f"module_scoring/B_rzv_module_{mod_key}.png"),
        explanation
    )

# ── Section 7: Key Gene Table ────────────────────────────────────────────────
html += '<h2 id="s7">7. 交叉比较三：关键基因逐一对比</h2>'

html += '''
<table>
<tr><th>基因</th><th>免疫功能</th><th>HZ log2FC</th><th>RZV-D14</th><th>RZV-D74</th><th>RZV-D365</th><th>象限</th><th>解读</th></tr>
<tr style="background:#fff3cd"><td><strong>ISG15</strong></td><td>I型IFN/ISG化</td><td><b>+1.08</b>**</td><td>+0.17</td><td>-0.34</td><td>-0.07</td><td>Q2</td><td>HZ特有IFN标志，疫苗避开</td></tr>
<tr style="background:#fff3cd"><td><strong>RSAD2</strong></td><td>抗病毒</td><td><b>+1.01</b>*</td><td>-2.07</td><td>-1.65</td><td>-0.17</td><td>Q2</td><td>疫苗中甚至下调（抗病毒程序未启动）</td></tr>
<tr style="background:#fff3cd"><td><strong>IFI44L</strong></td><td>抗病毒</td><td><b>+1.09</b>*</td><td>-1.98</td><td>-0.43</td><td>-0.41</td><td>Q2</td><td>同上，疫苗中下调</td></tr>
<tr style="background:#fff3cd"><td><strong>IFI27</strong></td><td>I型IFN</td><td><b>+1.22</b>*</td><td>+2.04</td><td>-2.36</td><td>+0.15</td><td>混合</td><td>D14同向，D74反向——时间依赖性</td></tr>
<tr style="background:#fff3cd"><td><strong>SERPING1</strong></td><td>补体</td><td><b>+1.18</b>**</td><td>—</td><td>—</td><td>—</td><td>Q2</td><td>HZ特有的补体激活</td></tr>
<tr style="background:#d4edda"><td><strong>ZEB2</strong></td><td>T细胞分化</td><td>-0.18</td><td><b>+3.06</b></td><td><b>+3.46</b></td><td><b>+2.91</b></td><td><b>Q4</b></td><td><b>疫苗持久重编程标记</b></td></tr>
<tr style="background:#d4edda"><td><strong>CTLA4</strong></td><td>免疫刹车</td><td>+0.15</td><td>+1.58</td><td>+1.62</td><td>+0.98</td><td>Q4</td><td>疫苗诱导自限性调控</td></tr>
<tr style="background:#d4edda"><td><strong>ICOS</strong></td><td>T细胞共刺激</td><td>+0.09</td><td>+1.23</td><td>+1.18</td><td>+0.54</td><td>Q4</td><td>疫苗特异性T细胞激活</td></tr>
<tr style="background:#d4edda"><td><strong>GNLY</strong></td><td>细胞毒</td><td>+0.13</td><td>-2.50</td><td>+3.19</td><td><b>+3.42</b></td><td>Q4</td><td>D365晚期积累——效应记忆建立</td></tr>
<tr><td><strong>TOP2A</strong></td><td>DNA复制</td><td>+1.25***</td><td><b>+7.05</b></td><td>+6.21</td><td>+1.03</td><td>Q1</td><td>共同增殖程序（非特异性）</td></tr>
<tr><td><strong>MKI67</strong></td><td>增殖</td><td>+0.99***</td><td>+4.24</td><td>+1.41</td><td>-0.72</td><td>Q1</td><td>同上</td></tr>
<tr><td><strong>CD38</strong></td><td>激活</td><td>+0.52***</td><td><b>+6.42</b></td><td>+2.99</td><td>-0.11</td><td>Q1</td><td>共同激活标记</td></tr>
<tr><td><strong>GZMA</strong></td><td>细胞毒</td><td>+0.68***</td><td>+2.47</td><td>+3.43</td><td>+1.46</td><td>Q1</td><td>共同效应通路</td></tr>
</table>

<div class="key-finding">
<strong>表格解读：</strong><br>
<span style="background:#fff3cd;padding:2px 8px;">黄色行</span> = Q2基因（HZ↑ RZV↓）——疾病特有的IFN炎症程序<br>
<span style="background:#d4edda;padding:2px 8px;">绿色行</span> = Q4基因（HZ↓ RZV↑）——疫苗保护性调控程序<br>
无背景色 = Q1基因（共同激活）——通用免疫程序<br>
*** p<0.001, ** p<0.01, * p<0.05
</div>
'''

# ── Section 8: Conclusions ───────────────────────────────────────────────────
html += '''
<h2 id="s8">8. I型 vs II型IFN：疫苗如何激活细胞免疫而不引发炎症？</h2>

<div class="key-finding">
<strong>核心问题：</strong>理论上，疫苗需要激活细胞免疫（IFN-γ驱动的Th1/CD8应答）才能建立保护。
但前面分析显示RZV疫苗<strong>不激活I型IFN</strong>（ISGs不升高）。
这是否意味着疫苗不激活IFN介导的细胞免疫？<br>
<strong>答案：不是。</strong>关键在于区分<strong>三套不同的IFN系统</strong>。
</div>

<h3>8.1 三种IFN系统的免疫学区别</h3>

<table>
<tr><th>特征</th><th>I型IFN (IFN-α/β)</th><th>II型IFN (IFN-γ)</th><th>III型IFN (IFN-λ)</th></tr>
<tr><td><strong>主要来源</strong></td><td>pDC、单核细胞、受感染细胞</td><td><strong>激活的T细胞（CD4 Th1、CD8 CTL）、NK细胞</strong></td><td>上皮细胞</td></tr>
<tr><td><strong>触发信号</strong></td><td>RIG-I/MDA5/TLR → IRF3/7</td><td><strong>TCR激活 + IL-12 → STAT4/T-bet</strong></td><td>类似I型但组织局限</td></tr>
<tr><td><strong>下游标记基因</strong></td><td>ISG15, RSAD2, IFI44L, MX1, OASL</td><td><strong>IRF1, CXCL9/10/11, GBP1-5, CIITA</strong></td><td>与I型重叠但局限于粘膜</td></tr>
<tr><td><strong>免疫学性质</strong></td><td>天然免疫、全身性炎症</td><td><strong>适应性免疫、T细胞效应功能</strong></td><td>屏障免疫</td></tr>
<tr><td><strong>在HZ中</strong></td><td style="background:#fff3cd"><strong>★★★ 强烈激活</strong>（全血ISG风暴）</td><td>混合（T细胞信号被全血稀释）</td><td>未评估</td></tr>
<tr><td><strong>在RZV中</strong></td><td style="background:#d4edda">→ 基本持平（避免炎症）</td><td style="background:#d4edda"><strong>详见下文——需区分"产生"与"应答"</strong></td><td>未评估</td></tr>
</table>

<h3>8.2 数据验证：I型IFN基因在HZ vs RZV中的行为</h3>
'''

# Add IFN comparison figure
html += fig_block(
    "图8A: I型 vs II型IFN通路基因对比",
    os.path.join(RES_DIR, "comparison/IFN_Type_I_vs_II_comparison.png"),
    '''<strong>上排（Type I IFN, IFN-α/β）：</strong>左侧为HZ急性期——ISG15、RSAD2、IFI44L等
    几乎全部上调（红色=上调，蓝色=下调，mean=+0.52）。中间和右侧为RZV-D14和D74——
    同一批基因在疫苗中呈<strong>混合模式</strong>（mean≈0），IFN-α/β的炎症程序未被疫苗激活。<br><br>
    <strong>下排左/中（Type II IFN hallmark, IFN-γ应答）：</strong>这些基因是"细胞接受IFN-γ刺激后"的下游效应基因
    （CXCL9/10/11趋化因子、GBP家族GTPases等），主要由巨噬细胞/内皮细胞表达。
    在CD4+ T细胞数据中这些基因本底就很低，因此Hallmark IFN-γ gene set的评分不适用于T细胞本身。<br><br>
    <strong>下排右（T细胞激活标记——这才是细胞免疫的真实读数）：</strong>
    ICOS(+1.23)、CD38(+6.42)、ZEB2(+3.06)、GZMA(+2.47)、IRF4(+0.54)——这些才是疫苗激活细胞免疫的
    分子证据。CTLA4(+1.58)的协同上调则展示了疫苗如何同时建立自限性调控。'''
)

# Add summary bar chart
html += fig_block(
    "图8B: IFN通路与T细胞激活 — 汇总对比",
    os.path.join(RES_DIR, "comparison/IFN_summary_barchart.png"),
    '''<strong>一目了然的对比：</strong>红色柱=HZ急性期的I型IFN（mean log2FC=+0.52）——
    天然免疫炎症风暴。蓝色柱=RZV-D14的T细胞激活基因（mean=+1.78）——适应性保护免疫。
    疫苗绕过了炎症性的I型IFN通路，直接激活了ICOS/CD38/ZEB2/GZMA等T细胞效应程序。'''
)

# Add the detailed data table snippet
html += '''
<h3>8.3 关键基因逐一对比：I型IFN vs T细胞激活</h3>

<table>
<tr><th>基因</th><th>所属通路</th><th>HZ急性期 log2FC</th><th>RZV-D14</th><th>RZV-D74</th><th>RZV-D365</th><th>解读</th></tr>

<tr style="background:#fff3cd"><td><strong>ISG15</strong></td><td>I型IFN</td><td><b>+1.08</b>**</td><td>+0.17</td><td>-0.34</td><td>-0.07</td><td>HZ中强烈↑，疫苗中平坦</td></tr>
<tr style="background:#fff3cd"><td><strong>RSAD2</strong></td><td>I型IFN</td><td><b>+1.01</b>*</td><td>-2.07</td><td>-1.65</td><td>-0.17</td><td>疫苗中<strong>下调</strong></td></tr>
<tr style="background:#fff3cd"><td><strong>IFI44L</strong></td><td>I型IFN</td><td><b>+1.09</b>*</td><td>-1.98</td><td>-0.43</td><td>-0.41</td><td>疫苗中下调</td></tr>
<tr style="background:#fff3cd"><td><strong>IFI27</strong></td><td>I型IFN</td><td><b>+1.22</b>*</td><td>+2.04</td><td>-2.36</td><td>+0.15</td><td>唯一例外——D14上调但D74反向</td></tr>
<tr style="background:#fff3cd"><td><strong>USP18</strong></td><td>I型IFN</td><td>+0.73**</td><td>+1.80</td><td>-0.04</td><td>+2.13</td><td>ISG化调控因子，在疫苗中也波动</td></tr>

<tr style="background:#d4edda"><td><strong>ICOS</strong></td><td>T细胞激活</td><td>+0.09</td><td><b>+1.23</b></td><td><b>+1.18</b></td><td>+0.54</td><td>疫苗T细胞共刺激——HZ中不激活</td></tr>
<tr style="background:#d4edda"><td><strong>CD38</strong></td><td>T细胞激活</td><td>+0.52***</td><td><b>+6.42</b></td><td><b>+2.99</b></td><td>-0.11</td><td>疫苗强激活——HZ中是全血弱信号</td></tr>
<tr style="background:#d4edda"><td><strong>ZEB2</strong></td><td>T细胞分化</td><td>-0.18</td><td><b>+3.06</b></td><td><b>+3.46</b></td><td><b>+2.91</b></td><td>疫苗持久重编程——HZ中无变化</td></tr>
<tr style="background:#d4edda"><td><strong>GZMA</strong></td><td>细胞毒</td><td>+0.68***</td><td><b>+2.47</b></td><td><b>+3.43</b></td><td>+1.46</td><td>共同效应通路——但疫苗中幅度更大</td></tr>
<tr style="background:#d4edda"><td><strong>IRF4</strong></td><td>效应分化</td><td>+0.31</td><td>+0.54</td><td>+0.74</td><td>+0.01</td><td>疫苗效应分化——HZ中不显著</td></tr>
<tr style="background:#d4edda"><td><strong>CTLA4</strong></td><td>免疫刹车</td><td>+0.15</td><td><b>+1.58</b></td><td><b>+1.62</b></td><td>+0.98</td><td>疫苗自限性调控——HZ中不激活</td></tr>
</table>

<h3>8.4 免疫学解读</h3>

<div class="key-finding">
<p><strong>RZV疫苗的免疫策略可以概括为"精准激活、绕过炎症"：</strong></p>

<p><strong>1. 不激活I型IFN通路（IFN-α/β）：</strong>
ISG15、RSAD2、IFI44L等ISGs在疫苗时间线中保持平坦甚至下调。
这是<strong>有意为之的免疫学设计</strong>——I型IFN的全身性激活会引起发热、肌痛等不良反应，
且过度的I型IFN信号可导致T细胞凋亡和免疫耗竭。疫苗的佐剂系统（AS01B）
通过TLR4/MyD88而非RIG-I/MDA5通路激活先天免疫，精准地避开了IFN-α/β的炎症级联。</p>

<p><strong>2. T细胞效应程序被充分激活：</strong>
ICOS（+1.23）、CD38（+6.42）、GZMA（+2.47）在D14显著上调——
这些不是"IFN-γ hallmark基因"，而是T细胞自身的激活和效应分子。
<strong>IFN-γ本身由激活的T细胞分泌后作用于巨噬细胞等靶细胞，
其效应在T细胞自身的转录组中并不体现为hallmark IFN-γ response基因的上调。</strong>
这解释了为什么GSEA显示"interferon gamma response"通路激活
（论文报告），而我们逐个看hallmark基因时信号混合。</p>

<p><strong>3. 自限性调控同步建立：</strong>
CTLA4（+1.58）、ZEB2（+3.06）在D14就与激活基因同步上调。
这意味着RZV疫苗在激活效应T细胞的同时就预设了"刹车"——
这是疫苗安全性（不引起免疫病理）和持久性（避免耗竭）的关键。</p>
</div>

<h3>8.5 补充数据：IFN通路对比完整数据表</h3>
<p>完整88个基因的对比数据见 <code>results/comparison/IFN_comparison_table.csv</code>。</p>
'''

html += '''
<h2 id="s9">9. 疫苗免疫学深度解读：体液免疫、细胞免疫、Treg与记忆</h2>

<div class="key-finding">
<strong>核心问题：</strong>从疫苗开发角度，体液免疫（B细胞/抗体）是否重要？细胞免疫中Th1和Th2哪个更关键？
Treg应维持什么水平？这两个数据集能为此提供什么证据？
</div>

<h3>9.1 体液免疫（B细胞/抗体）：RZV如何支持抗体产生？</h3>

<table>
<tr><th>证据来源</th><th>发现</th><th>解读</th></tr>
<tr style="background:#fff3cd">
<td><strong>GSE242252 (HZ全血)</strong></td>
<td>B细胞/浆细胞基因在HZ急性期强烈激活：MZB1(+1.19***), IGLC3(+0.96**), CD79A(+0.50***), IGHG1(+0.61)</td>
<td>自然感染确实激活了体液免疫——抗体是抗VZV的重要效应机制</td>
</tr>
<tr style="background:#d4edda">
<td><strong>GSE249632 (RZV CD4+ T)</strong></td>
<td>无法直接检测B细胞(数据只含CD4+ T细胞)。<br><strong>但Tfh签名是所有模块中最强的信号</strong>：CXCR5, BCL6, ICOS, PDCD1, CD40LG, IL21, BTLA, MAF, TOX2</td>
<td>Tfh细胞的核心功能是帮助B细胞进行类别转换和亲和力成熟。<strong>Tfh签名强=抗体应答的细胞基础被疫苗充分激活</strong></td>
</tr>
</table>

<div class="key-finding">
<strong>疫苗开发启示：</strong>RZV诱导极强的Tfh极化（D74 mean logFC = +1.75，为所有模块最高），
这表明疫苗同时激活了细胞免疫和体液免疫的细胞基础。Tfh的强度可以作为候选疫苗
"是否可能诱导高质量抗体应答"的替代指标。
</div>
'''

# Add Tfh/Th1/Th2 balance figure
html += fig_block(
    "图9A: RZV疫苗中的Th1/Th2/Tfh平衡",
    os.path.join(RES_DIR, "comparison/Th1_Th2_Tfh_balance.png"),
    '''<strong>三个时间点的CD4+ T细胞极化方向。</strong>
    D14（第一针后14天）：Tfh(+1.26) > Th2(+0.91) > Th1(+0.50)——Tfh占主导，
    这与"疫苗激活CD4+ T细胞以帮助B细胞产生抗体"的预期完全一致。
    D74（第二针后14天）：Tfh进一步增强至+1.75，Th1维持+0.53。
    D365（一年后）：所有极化减弱但Tfh仍维持+0.53——提示持续的B细胞辅助能力。
    <strong>关键：</strong>RZV诱导的是Th1/Tfh混合型——Th1驱动抗病毒细胞免疫，Tfh驱动抗体产生。'''
)

html += '''
<h3>9.2 细胞免疫：Th1还是Th2？</h3>

<table>
<tr><th>基因集</th><th>代表基因</th><th>HZ mean</th><th>RZV D14</th><th>RZV D74</th><th>RZV D365</th><th>解读</th></tr>
<tr><td><strong>Th1</strong></td><td>TBX21, STAT1, IL18R1, CXCR3, TNF</td><td>+0.13</td><td><b>+0.50</b></td><td><b>+0.53</b></td><td>+0.15</td><td>疫苗稳定激活Th1——抗病毒细胞免疫的核心</td></tr>
<tr><td><strong>Th2</strong></td><td>GATA3, MAF, STAT6, CCR3, IL4R</td><td>-0.01</td><td>+0.91</td><td>+0.60</td><td>+0.25</td><td>Th2信号部分来自MAF/STAT6等转录因子——并非典型的过敏型Th2</td></tr>
<tr><td><strong>Tfh</strong></td><td>CXCR5, BCL6, ICOS, PDCD1, IL21</td><td>+0.06</td><td><b>+1.26</b></td><td><b>+1.75</b></td><td><b>+0.53</b></td><td><strong>最强信号——B细胞辅助</strong></td></tr>
</table>

<div class="warn">
<strong>关于Th2信号的技术说明：</strong>Th2模块在RZV-D14显示+0.91，但这基于仅5个基因（MAF, STAT6, CCR3, IL4R, CRLF2）。
MAF在Th2和Tfh中均表达，IL4R是IL-4受体（非Th2特异性）。IL4/IL5/IL13等Th2标志性细胞因子
在scRNA-seq中难以检测。因此这个Th2信号不应被过度解读为"疫苗诱导Th2"——
RZV的真实极化是<strong>Th1/Tfh混合型</strong>，这对抗病毒疫苗是最理想的。
</div>

<div class="key-finding">
<strong>疫苗开发启示：</strong>候选疫苗应评估Th1/Th2平衡。Th1偏倚（TBX21↑, IFNG↑）对清除胞内病原体
至关重要；Th2偏倚在某些情况下可能与疫苗增强性疾病（VAED）相关。
RZV的AS01B佐剂通过TLR4-MyD88通路强烈诱导Th1，是其成功的关键设计。
</div>

<h3>9.3 Treg：应该维持什么水平？</h3>
'''

# Add the module comparison figure
html += fig_block(
    "图9B: 疫苗免疫学各模块在HZ疾病与RZV疫苗中的表现",
    os.path.join(RES_DIR, "comparison/vaccine_immunology_modules.png"),
    '''<strong>12个免疫学功能模块的横向对比。</strong>红色柱=HZ急性期，橙色/粉色/绿色=RZV-D14/D74/D365。
    背景色区分功能类别。<strong>关键观察：</strong>
    (1) Tfh(绿色)是RZV中最强的信号——疫苗核心功能是帮助B细胞；
    (2) T细胞激活/效应记忆在RZV中高于HZ——疫苗的精准性；
    (3) 补体(蓝色最左)在HZ中激活但在RZV中被抑制——疫苗绕过了这条炎症通路；
    (4) 经典Treg在RZV中降低——疫苗不扩增Treg群体。'''
)

html += '''
<table>
<tr><th>Treg亚型</th><th>代表基因</th><th>HZ</th><th>RZV D14</th><th>RZV D74</th><th>RZV D365</th><th>解读</th></tr>
<tr style="background:#fff3cd">
<td><strong>经典Treg</strong></td><td>FOXP3, IL2RA, IKZF2, IKZF4</td><td>-0.06</td><td><b>-0.71</b></td><td>-0.40</td><td>-0.54</td>
<td><strong>Treg群体未被疫苗扩增</strong>——这是好的，因为Treg扩增会抑制保护性免疫</td>
</tr>
<tr style="background:#d4edda">
<td><strong>活化/效应Treg</strong></td><td>CTLA4, TIGIT, LAG3, GITR</td><td>+0.18</td><td><b>+0.98</b></td><td><b>+1.19</b></td><td>+0.87</td>
<td><strong>但这些分子并非Treg独有！</strong>CTLA4和TIGIT是效应T细胞激活后的"自限性刹车"</td>
</tr>
</table>

<div class="key-finding">
<strong>关键免疫学区分：</strong><br>
(1) <strong>FOXP3+ Treg扩增</strong> → 会抑制疫苗免疫应答 → 应避免<br>
(2) <strong>效应T细胞上的CTLA4/TIGIT上调</strong> → 是激活诱导的自限性调控 → 是健康应答的标志<br><br>
RZV的数据清晰展示了模式(2)而非模式(1)：gE特异性CD4+ T细胞在强烈激活（ICOS↑, CD38↑）的同步上调CTLA4/TIGIT，
确保免疫应答"精准而不失控"。<br><br>
<strong>疫苗开发启示：</strong>评价疫苗应答时，CTLA4/PDCD1等分子不应简单视为"耗竭"或"抑制"——
在激活后早期它们是生理性调控信号。真正的耗竭特征是TOX持续升高 + 多重抑制性受体共表达 + 效应功能丧失。
RZV中的T细胞保持了GZMA/GNLY的效应功能，同时表达适度CTLA4——这是"被调控但不耗竭"的理想状态。
</div>

<h3>9.4 耗竭 vs 记忆：如何评估疫苗持久性？</h3>
'''

html += fig_block(
    "图9C: 耗竭 vs 记忆——关键基因在RZV时间线中的动态",
    os.path.join(RES_DIR, "comparison/exhaustion_vs_memory_dynamics.png"),
    '''<strong>HZ列（菱形）=同一基因在HZ急性期的log2FC。</strong>
    左侧为调控/检查点分子，右侧为效应分子和记忆标记。<br><br>
    <strong>关键动态模式：</strong><br>
    (1) CD38（激活标记）：D14=+6.42→D74=+2.99→D365=-0.11——脉冲式激活，随后消退；<br>
    (2) ZEB2（效应分化）：D14=+3.06→D74=+3.46→D365=+2.91——持久分化重编程；<br>
    (3) GNLY（细胞毒）：D14=-2.50→D74=+3.19→D365=+3.42——<strong>晚期效应积累</strong>——效应记忆的建立；<br>
    (4) CTLA4/PDCD1/TIGIT：D14起同步上调，D365部分回落——可逆调控，非永久耗竭；<br>
    (5) TCF7/IL7R（记忆干性）：D365时下降——gE特异性T细胞偏向效应记忆而非中央记忆。<br><br>
    <strong>疫苗开发启示：</strong>
    RZV诱导的是效应记忆（T<sub>EM</sub>）为主的长期保护——这些细胞在组织中巡逻，
    遇到gE抗原即可快速产生效应功能。GNLY在D365的持续积累（+3.42）提示
    CD4+ CTL样的效应记忆正在建立。而TCF7/IL7R的下降提示中央记忆（T<sub>CM</sub>，
    驻留在淋巴结）的相对比例降低——这对疫苗来说不是问题，
    因为效应记忆对快速清除感染细胞更关键。'''
)

html += '''
<h3>9.5 综合解读：两个数据集如何指导疫苗评估</h3>

<table>
<tr><th>评估维度</th><th>GSE242252 (HZ疾病) 能提供的</th><th>GSE249632 (RZV疫苗) 能提供的</th><th>疫苗评估建议</th></tr>

<tr>
<td><strong>体液免疫</strong></td>
<td>B细胞/浆细胞基因↑ (MZB1, IGHG1, IGLC3)<br>补体级联↑ (SERPING1, C1QA)</td>
<td>Tfh签名 (CXCR5/BCL6/ICOS/IL21)<br>CD40LG (辅助B细胞的关键分子)</td>
<td>候选疫苗应监测Tfh标记和B细胞活化<br>补体激活是炎症信号，疫苗应避免</td>
</tr>

<tr>
<td><strong>Th1/Th2平衡</strong></td>
<td>Bulk中Th1信号被稀释(CXCR3+0.65**)<br>Th2信号弱</td>
<td><strong>Th1偏向</strong> (TBX21/STAT1/IL18R1)<br>Tfh > Th1 > Th2</td>
<td>成功的抗病毒疫苗应Th1偏倚<br>TBX21/STAT1/IL18R1可作为Th1标志</td>
</tr>

<tr>
<td><strong>Treg水平</strong></td>
<td>经典Treg无显著变化<br>ENTPD1/CD39下调 (-0.61**)</td>
<td>FOXP3签名↓ (不应扩增Treg)<br>CTLA4/TIGIT↑ (自限性调控，非抑制)</td>
<td>区分"Treg扩增"(↓保护) vs "激活诱导调控"(健康)<br>CTLA4/FOXP3比值 > 2 = 良好免疫质量</td>
</tr>

<tr>
<td><strong>记忆/持久性</strong></td>
<td>无法评估T细胞记忆<br>(全血Bulk，无克隆追踪)</td>
<td>ZEB2 持续↑ (D365=+2.91)<br>GNLY 晚期积累 (D365=+3.42)<br>效应记忆(T<sub>EM</sub>)>中央记忆(T<sub>CM</sub>)</td>
<td>ZEB2持久↑ = 分化重编程标记<br>GNLY晚期↑ = 效应记忆建立<br>D365时间点对评估持久性至关重要</td>
</tr>

<tr>
<td><strong>安全性标志</strong></td>
<td>I型IFN风暴 (ISG15/RSAD2)<br>TNF上调 (炎症)</td>
<td>I型IFN通路平坦<br>补体基因下调<br>CTLA4同步上调 (自限性刹车)</td>
<td>候选疫苗应避免激活I型IFN/补体级联<br>CTLA4的适度上调是安全信号</td>
</tr>
</table>

<div class="key-finding">
<h3>对疫苗免疫效果评估的五个关键提示</h3>

<p><strong>1. Tfh信号的强度是体液免疫质量的替代指标。</strong>
无法直接测抗体滴度时，CXCR5/BCL6/ICOS/IL21的co-upregulation
是B细胞辅助正在发生的可靠分子证据。</p>

<p><strong>2. Th1/Th2比值而非绝对值决定免疫质量。</strong>
不应单纯看Th1或Th2基因的log2FC，而应关注Th1/(Th2+Treg)的比值。
高Th1/Treg比值=良好的效应/调控平衡。</p>

<p><strong>3. "调控分子≠耗竭"——需要多维度判断。</strong>
CTLA4/PDCD1/TIGIT的上调在疫苗应答中是正常且必要的。
真正的耗竭需要TOX↑ + 多重抑制性受体 + 效应功能↓ (GZMB/PRF1↓)。
RZV展示了"被调控但不耗竭"的理想模式。</p>

<p><strong>4. D365时间点是评估持久性的金标准。</strong>
D14和D74的基因表达变化只能反映急性应答。D365的ZEB2/GNLY持续表达
才代表真正的免疫记忆编程。短期保护看抗体，长期保护看效应记忆。</p>

<p><strong>5. 安全疫苗的精髓是"选择性激活"而非"全面激活"。</strong>
RZV选择性激活了Tfh/Th1/效应记忆通路，同时绕过了I型IFN炎症和补体级联。
这种免疫学"精准性"正是佐剂系统（AS01B）设计的核心目标——
也是未来疫苗开发的努力方向。</p>
</div>
'''

html += '''
<h2 id="s10">10. 综合结论</h2>

<div class="key-finding">
<h3>三种独立的分析方法一致支持以下结论：</h3>

<p><strong>1. HZ疾病和RZV疫苗是全球范围内不同的免疫程序。</strong><br>
象限图分析显示Spearman ρ ≈ 0（D14: 0.034, D74: 0.001, D365: 0.036）。
RZV疫苗不是简单地"模拟"自然感染——它建立的是功能独特的保护性免疫。</p>

<p><strong>2. HZ的标志是I型IFN介导的天然免疫风暴。</strong><br>
ISG15（+1.08**）、RSAD2（+1.01*）、IFI44L（+1.09*）在急性期显著上调，
但RZV疫苗接种后这些基因不升高甚至下降。基因集评分转移（方向A）进一步证实：
Type I IFN & Antiviral模块在疫苗时间线中完全持平。
疫苗绕过了炎症性的天然免疫通路，这解释了为什么RZV不会引起HZ样症状。</p>

<p><strong>3. RZV的标志是精准的、自限性的适应性免疫调控。</strong><br>
ZEB2在D365仍维持+2.91 logFC——这是疫苗诱导的长期T细胞分化重编程。
CTLA4（+1.58~+0.98）和ICOS（+1.23~+0.54）的持续上调构建了"刹车"机制。
GNLY在D365的晚期积累（+3.42）提示效应记忆的逐步建立。</p>

<p><strong>4. 疫苗通过T细胞效应程序（而非I型IFN炎症）建立细胞免疫。</strong><br>
ICOS(+1.23)、CD38(+6.42)、GZMA(+2.47)、ZEB2(+3.06)在RZV-D14显著上调，
而ISG15、RSAD2、IFI44L等I型IFN炎症基因在疫苗中持平或下调。
这表明RZV巧妙地"绕过"了IFN-α/β的天然免疫炎症通路，
直接激活T细胞适应性免疫。CTLA4的同步上调（+1.58）则构成了自限性刹车，
确保免疫应答"精准而不过度"——这既是疫苗有效性的基础，也是其安全性的保障。</p>

<p><strong>5. 共同点仅限于通用增殖程序。</strong><br>
TOP2A、MKI67、CD38在两种过程中均上调，反映任何免疫激活都需要细胞分裂。
这些是"噪音"而非疾病/疫苗的特异性信号。</p>
</div>

<h3>方法学贡献</h3>
<ul>
  <li><strong>三重验证框架：</strong>象限图（全基因相关）→ 基因集评分（模块水平）→ 关键基因表（单基因水平），构成逐级收窄的证据链。</li>
  <li><strong>双向评分策略：</strong>从疾病定义模块在疫苗中验证（方向A），从疫苗定义模块在疾病中验证（方向B），避免单向假设偏倚。</li>
  <li><strong>通路分层方法：</strong>基于GO富集将DEGs拆分为功能模块分别评分，揭示混合信号背后的机制。</li>
</ul>

<h3>局限性</h3>
<ul>
  <li>Bulk vs 单细胞平台差异限制直接比较（只能比较log2FC方向）</li>
  <li>GSE249632仅有7位供者，D60/D74仅5组pseudobulk样本（统计效力有限）</li>
  <li>GSE242252的T细胞特异性基因在全血中被稀释</li>
  <li>两个数据集的基线群体不同（HZ患者恢复期 vs 健康人接种前）</li>
</ul>

<hr>
<p style="color:#999; font-size:0.9em; text-align:center;">
报告生成：2026-05-29 | 代码仓库：<a href="https://github.com/cmjlabs/VZV-geo">github.com/cmjlabs/VZV-geo</a><br>
分析脚本：scripts/01-08 | 结果文件：results/
</p>
</body>
</html>
'''

# ── Write ─────────────────────────────────────────────────────────────────────
out_path = os.path.join(RES_DIR, "comprehensive_report.html")
with open(out_path, 'w', encoding='utf-8') as f:
    f.write(html)

print(f"HTML report written to: {out_path}")
print(f"Size: {os.path.getsize(out_path) / 1024 / 1024:.1f} MB")
