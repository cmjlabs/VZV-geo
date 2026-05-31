#!/usr/bin/env python3
"""Generate Chapter 3 HTML report — using paper's published Supplementary Table S5."""
import os, base64

RES = "/media/cmj/MechanicalDisk/yjs/VZV-geo/results"
FIG = os.path.join(RES, "chapter3_figures")

def img(path, w="95%"):
    if not os.path.exists(path): return f'<p class="miss">[{os.path.basename(path)}]</p>'
    with open(path,'rb') as f:
        b64 = base64.b64encode(f.read()).decode()
    mime = 'image/png'
    return f'<img src="data:{mime};base64,{b64}" style="max-width:{w};display:block;margin:10px auto;border:1px solid #ddd;border-radius:4px;">'

def fig(title, path, legend, w="95%"):
    return f'<div class="fig"><h4>{title}</h4>{img(path,w)}<div class="leg">{legend}</div></div>'

html = '''<!DOCTYPE html>
<html lang="zh-CN">
<head><meta charset="UTF-8"><meta name="viewport" content="width=device-width,initial-scale=1.0">
<title>第三章 HZ疾病特征与RZV保护性免疫特征分析</title>
<style>
body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI','Noto Sans SC',sans-serif;max-width:1100px;margin:0 auto;padding:20px;color:#222;line-height:1.8;background:#fafafa}
h1{font-size:1.8em;border-bottom:3px solid #1a5276;padding-bottom:10px;color:#1a5276}
h2{font-size:1.3em;margin-top:40px;border-bottom:2px solid #2980b9;padding-bottom:5px;color:#2980b9}
h3{font-size:1.1em;margin-top:30px;color:#2c3e50}
.fig{background:#fff;border:1px solid #e0e0e0;border-radius:8px;padding:20px;margin:25px 0;box-shadow:0 2px 4px rgba(0,0,0,.05)}
.leg{background:#f8f9fa;border-left:4px solid #2980b9;padding:12px 16px;margin-top:12px;font-size:.92em;border-radius:0 4px 4px 0}
.box{background:#fff;border:1px solid #ddd;border-radius:8px;padding:20px;margin:20px 0}
.key{background:#d4edda;border-left:4px solid #28a745;padding:12px 16px;margin:15px 0;border-radius:0 4px 4px 0}
.note{background:#fff3cd;border-left:4px solid #ffc107;padding:12px 16px;margin:15px 0;border-radius:0 4px 4px 0}
.ref{background:#d1ecf1;border-left:4px solid #17a2b8;padding:12px 16px;margin:15px 0;border-radius:0 4px 4px 0}
table{border-collapse:collapse;width:100%;margin:15px 0;font-size:.88em}
th{background:#2980b9;color:#fff;padding:10px 8px;text-align:left}
td{padding:8px;border-bottom:1px solid #e0e0e0}
tr:nth-child(even){background:#f8f9fa}
.miss{color:#dc3545;font-style:italic}
.q{font-weight:bold;color:#1a5276;margin-top:10px}
.a{color:#333;margin-left:20px;margin-bottom:15px}
</style></head><body>

<h1>第三章 基于公共组学数据的HZ疾病特征与RZV保护性免疫特征分析</h1>

<div class="ref">
<strong>核心科学问题：</strong><br>
<strong>Q1：</strong>HZ发生时宿主系统免疫发生什么变化？→ <em>Disease Signature</em><br>
<strong>Q2：</strong>RZV疫苗诱导了什么样的CD4⁺ T细胞状态？→ <em>Protection Signature</em><br>
<strong>Q3：</strong>这些特征能否作为后续评价候选疫苗的参照标准？→ <em>Evaluation Framework</em>
</div>

<div class="note">
<strong>本章不回答：</strong>HZ和RZV谁更强；两者是否本质不同（现有数据无法严格支持因果推断）。<br>
<strong>本章回答：</strong>两个参照系各自的免疫特征，为第四、五、六章提供评价坐标系。
</div>

<h2>3.1 研究设计与总体思路</h2>
<div class="box">
<p>水痘-带状疱疹病毒（VZV）再激活导致带状疱疹（HZ），而重组带状疱疹疫苗（RZV, Shingrix）能有效预防HZ。
一个理想候选疫苗应<strong>（1）诱导有效适应性免疫；（2）避免疾病相关炎症程序。</strong></p>
<p>为此，利用公共数据建立两个免疫参照系：</p>
<table>
<tr><th>参照系</th><th>数据集</th><th>技术</th><th>样本</th><th>数据来源</th></tr>
<tr><td><strong>HZ疾病参照系</strong></td><td>GSE242252</td><td>Bulk RNA-seq (全血)</td><td>26例HZ患者配对</td><td>Vandoren et al. 2024, <em>J Infect Dis</em></td></tr>
<tr><td><strong>RZV保护参照系</strong></td><td>GSE249632</td><td>scRNA-seq (gE-CD4⁺ T)</td><td>7例接种者, 2,231细胞</td><td>GEO公共数据</td></tr>
<tr><td><em>HZ细胞来源验证</em></td><td>HRA008316</td><td>scRNA-seq (PBMCs)</td><td>3健康/3HZ/3恢复</td><td>Zheng et al. 2024, <em>Commun Biol</em></td></tr>
</table>
</div>

<h2>3.2 HZ急性期系统性免疫特征分析</h2>
<p><strong>数据：</strong>GSE242252 (Vandoren et al. 2024)。26例HZ患者全血Bulk RNA-seq，急性发作期 vs 恢复期(约1年后)。差异表达数据来自<strong>Supplementary Table S5</strong>（文章正式发表的Top 200 DEGs）。</p>

<h3>3.2.1 差异表达基因分析</h3>
<div class="q">HZ急性期哪些基因发生改变？</div>

''' + fig("图3.1 HZ急性期 vs 恢复期火山图", os.path.join(FIG,"Figure3_1_volcano.png"),
'''<strong>数据来源：</strong>Vandoren et al. 2024 Supplementary Table S5 (Top 200 DEGs, Active HZ vs HZ after 1 year)。<br>
<strong>结果：</strong>FDR&lt;0.1且|LFC|>0.58的标准下共106个DEGs（101上调，5下调）。
7个最显著上调基因为ISGs，由天然免疫系统在病毒识别时产生：<br>
<span style="color:#B22222"><strong>IFI44L (LogFC=+2.42**), IFI27 (+2.32**), RSAD2 (+2.04**), SERPING1 (+1.98**), ISG15 (+1.94**), SIGLEC1 (+1.73**), IFI44 (+1.52**)</strong></span><br>
细胞增殖基因同样显著：<strong>TOP2A (+1.43***), PTTG1 (+1.96***)</strong>；
B细胞/浆细胞标志：<strong>MZB1 (+1.62***), IGLC2 (+1.43***)</strong>。<br>
<strong>回答：</strong>HZ急性期表现出显著的I型IFN抗病毒炎症特征，同时伴随细胞增殖和体液免疫激活。''') + '''


<h3>3.2.2 GO富集分析</h3>
<div class="q">这些基因变化反映了什么生物学过程？</div>

''' + fig("图3.2 HZ上调基因GO富集气泡图", os.path.join(FIG,"GO_bubble_plot.png"),
'''<strong>结果：</strong>GO Biological Process富集（Enrichr, GO_2023）：
Defense Response To Virus（抗病毒防御, p=2.0e-3, 8 genes）、
B Cell Receptor Signaling Pathway（BCR信号, p=2.0e-3, 5 genes）、
Antigen Receptor-Mediated Signaling（抗原受体信号, p=8.7e-3, 6 genes）、
Mitotic Spindle Organization（有丝分裂纺锤体, p=8.7e-3, 5 genes）。<br>
<strong>回答：</strong>HZ同时激活天然免疫（抗病毒）、体液免疫（BCR）和淋巴细胞扩增——多层面系统性炎症应答。''') + '''


<h3>3.2.3 HZ疾病标志提取</h3>
<div class="key">
<strong>HZ Disease Signature（来自 Vandoren et al. 2024 Table S5）：</strong>
<table>
<tr><th>类别</th><th>代表基因</th><th>功能</th><th>LogFC</th><th>padj</th></tr>
<tr><td>I型IFN</td><td>IFI44L</td><td>抗病毒ISG</td><td>+2.42</td><td>1.8e-04</td></tr>
<tr><td>I型IFN</td><td>IFI27</td><td>抗病毒ISG</td><td>+2.32</td><td>8.4e-04</td></tr>
<tr><td>I型IFN</td><td>RSAD2</td><td>抗病毒效应</td><td>+2.04</td><td>9.7e-04</td></tr>
<tr><td>I型IFN</td><td>ISG15</td><td>ISG化修饰</td><td>+1.94</td><td>1.3e-03</td></tr>
<tr><td>补体</td><td>SERPING1</td><td>补体激活</td><td>+1.98</td><td>8.0e-04</td></tr>
<tr><td>增殖</td><td>PTTG1</td><td>细胞周期</td><td>+1.96</td><td>9.4e-12</td></tr>
<tr><td>增殖</td><td>TOP2A</td><td>DNA复制</td><td>+1.43</td><td>4.3e-12</td></tr>
<tr><td>浆细胞</td><td>MZB1</td><td>抗体分泌</td><td>+1.62</td><td>7.4e-05</td></tr>
</table>
该Signature用于后续评价候选疫苗是否避免了HZ样炎症程序。
</div>

<h3>3.2.4 单细胞验证：疾病标志的细胞来源</h3>
<div class="q">这些ISGs来自哪些细胞类型？</div>

''' + fig("图3.3 Bulk全血 vs scRNA-seq PBMCs ISG表达对比", os.path.join(FIG,"bulk_vs_scRNA_comparison.png"),
'''<strong>数据：</strong>HRA008316 (Zheng et al. 2024), PBMC scRNA-seq, HZ患者 vs 健康对照。<br>
<strong>结果：</strong>(1) ISGs在Bulk全血中显著上调（文章Table S5），但在PBMC scRNA-seq的T细胞中不显著或反向（ISG15 LogFC=-0.37 ns）。
(2) Classical Monocyte中IFI27在HZ患者表达量为健康对照的<strong>27.6倍</strong>，SERPING1为2.4倍。
(3) ISGs作为细胞标志出现在Basophils/B cells/pDCs/Neutrophils/T cells等多个细胞类型。<br>
<strong>回答：</strong>Bulk全血中的ISG信号来自多种天然免疫细胞（单核细胞、pDC、中性粒细胞），<strong>而非T细胞</strong>。这解释了为什么靶向T细胞的RZV疫苗不触发ISG风暴。''') + '''


<h2>3.3 RZV诱导的gE特异性CD4⁺ T细胞特征</h2>
<p><strong>数据：</strong>GSE249632。SMART-Seq v4 scRNA-seq，7例健康疫苗接种者，2,231个QC通过的gE四聚体⁺CD4⁺ T细胞，纵向追踪D0→D14→D60→D74→D365。</p>

<h3>3.3.1 DEG动力学分析</h3>
<div class="q">RZV应答是持续炎症还是可控激活？</div>

''' + fig("图3.4 RZV疫苗CD4⁺ T细胞DEG时间线", os.path.join(FIG,"DEG_timeline_barchart.png"),
'''<strong>结果：</strong>
<table>
<tr><th>时间点</th><th>上调</th><th>下调</th><th>解读</th></tr>
<tr><td>D14 vs D0</td><td>165</td><td>178</td><td>第一针后强激活</td></tr>
<tr><td>D60 vs D0</td><td>9</td><td>4</td><td>几乎回到基线</td></tr>
<tr><td>D74 vs D0</td><td>68</td><td>40</td><td>第二针再激活（弱于第一针）</td></tr>
<tr><td>D365 vs D0</td><td>18</td><td>22</td><td>长期印记残留（40 DEGs）</td></tr>
</table>
<strong>回答：</strong>脉冲式、可调控的适应性免疫——非持续炎症。''') + '''


<h3>3.3.2 长期维持基因分析</h3>
<div class="q">疫苗留下了什么长期印记？</div>

''' + fig("图3.5 五基因纵向轨迹：HZ vs RZV", os.path.join(FIG,"five_gene_trajectory.png"),
'''<strong>本章核心图。</strong>灰色虚线=HZ急性期水平（来自Table S5），彩色折线=RZV疫苗接种后（D14/D60/D74/D365）log2FC变化。<br><br>
<strong>ZEB2（T细胞分化重编程）：</strong>HZ中无显著变化（-0.18 ns），RZV后持续高达一年 — <strong>疫苗特有标志</strong><br>
<strong>CTLA4（免疫自限性）：</strong>HZ中无变化（+0.14 ns），RZV后持续上调 — <strong>主动建立免疫刹车</strong><br>
<strong>ICOS（T细胞共刺激）：</strong>HZ中无变化（+0.09 ns），RZV后上调 — <strong>协同激活信号</strong><br>
<strong>HAVCR2/TIM-3（效应记忆）：</strong>D365仍维持 — <strong>长期效应记忆标志</strong><br>
<strong>ISG15（I型IFN）：</strong>HZ中显著上调（+1.94**），RZV中完全不动 — <strong>疫苗不触发炎症性IFN通路</strong><br><br>
<strong>回答：</strong>RZV的长期印记是精准的T细胞分化重编程和自限性调控，而非炎症。''') + '''


<h3>3.3.3 功能富集分析</h3>
<div class="q">疫苗诱导的基因对应什么功能？</div>

''' + fig("图3.6 RZV D14上调基因Hallmark富集", os.path.join(FIG,"RZV_hallmark_bubble.png"),
'''<strong>结果：</strong>MSigDB Hallmark富集（D14上调基因，FDR&lt;0.05）：
G2-M Checkpoint (p=3.2e-4, 8 genes), E2F Targets (p=9.4e-3, 6 genes),
Mitotic Spindle (p=2.9e-2, 5 genes), Myc Targets V1 (p=2.9e-2, 5 genes)。<br>
主要为细胞周期/增殖通路。T细胞分化记忆信号（ZEB2、CTLA4）是<strong>单个关键调控因子</strong>驱动，非通路级别重编程。<br>
<strong>回答：</strong>RZV主要重塑CD4⁺ T细胞增殖和记忆状态。保护性信号体现在关键基因的持久表达。''') + '''


<h3>3.3.4 RZV保护标志提取</h3>
<div class="key">
<strong>RZV Protection Signature：</strong>
<table>
<tr><th>类别</th><th>基因</th><th>功能</th><th>D14 LFC</th><th>D365 LFC</th></tr>
<tr><td>记忆分化</td><td>ZEB2</td><td>T细胞分化重编程</td><td>+3.06</td><td>+2.91</td></tr>
<tr><td>共刺激</td><td>ICOS</td><td>T细胞协同激活</td><td>+1.23</td><td>+0.54</td></tr>
<tr><td>免疫调节</td><td>CTLA4</td><td>自限性检查点</td><td>+1.58</td><td>+0.98</td></tr>
<tr><td>效应记忆</td><td>HAVCR2</td><td>TIM-3/长期效应</td><td>—</td><td>维持</td></tr>
</table>
该Signature用于后续评价候选疫苗是否诱导了类RZV保护性免疫。
</div>


<h2>3.4 HZ与RZV参照系的关联分析</h2>
<h3>3.4.1 两个参照系的生物学层级</h3>
<div class="box">
<table>
<tr><th>维度</th><th>HZ疾病参照系</th><th>RZV保护参照系</th></tr>
<tr><td>免疫类型</td><td>系统性天然免疫炎症</td><td>特异性CD4⁺适应性免疫</td></tr>
<tr><td>状态</td><td>疾病状态（被动应答）</td><td>保护状态（主动编程）</td></tr>
<tr><td>细胞来源</td><td>全血（多细胞）→ 单核/DC/中性粒</td><td>Tetramer⁺CD4⁺ T（抗原特异性）</td></tr>
<tr><td>核心信号</td><td>I型IFN、补体、增殖</td><td>T细胞分化、自限性调控</td></tr>
<tr><td>持续时间</td><td>急性（数周）</td><td>长期（≥1年）</td></tr>
<tr><td>数据来源</td><td>Vandoren et al. 2024 Table S5</td><td>GSE249632 (本分析)</td></tr>
</table>
<div class="note">两数据集不可直接进行机制比较（平台不同、细胞群体不同）。本章对比是<strong>特征层面的参照</strong>。</div>
</div>

<h3>3.4.2 候选疫苗评价框架</h3>
<div class="key">
<strong>基于两个参照系，建立候选疫苗GE282+GB705的双重评价标准：</strong>
<table>
<tr><th>标准</th><th>参照系</th><th>核心指标</th><th>期望方向</th><th>对应章节</th></tr>
<tr><td><strong>标准1：是否避免HZ样炎症？</strong></td><td>HZ Disease Signature</td><td>IFI44L, IFI27, RSAD2, ISG15, SERPING1</td><td>不激活/低表达</td><td>第四章（固有免疫）</td></tr>
<tr><td><strong>标准2：是否诱导RZV样CD4应答？</strong></td><td>RZV Protection Signature</td><td>ZEB2, CTLA4, ICOS, HAVCR2</td><td>持久上调</td><td>第五、六章</td></tr>
</table>
该框架将第三章从"公共数据库分析"转变为<strong>整篇论文的理论坐标系</strong>。
</div>


<h2>3.5 本章小结</h2>
<div class="box">
<h3>建立的两个参照系</h3>
<table>
<tr><th></th><th>HZ Disease Signature</th><th>RZV Protection Signature</th></tr>
<tr><td><strong>核心特征</strong></td><td>I型IFN风暴 / 补体激活 / 体液应答增强</td><td>CD4记忆重塑 / 长期分化维持 / 自限性调控</td></tr>
<tr><td><strong>代表基因</strong></td><td>IFI44L(+2.42), RSAD2(+2.04), ISG15(+1.94), SERPING1(+1.98)</td><td>ZEB2(+3.06→+2.91), ICOS, CTLA4, HAVCR2</td></tr>
<tr><td><strong>数据来源</strong></td><td>GSE242252 (Vandoren et al. 2024, Table S5)</td><td>GSE249632 (本分析, n=7供者)</td></tr>
<tr><td><strong>细胞来源验证</strong></td><td>HRA008316: ISGs来自单核/DC/中性粒, 非T细胞</td><td>gE四聚体⁺CD4⁺ T细胞 (抗原特异性)</td></tr>
</table>

<h3>核心结论</h3>
<ol>
<li>HZ急性期是<strong>I型IFN介导的天然免疫主导的系统性炎症</strong>（IFI44L, RSAD2, ISG15均LogFC>+1.9），ISG信号来自多种天然免疫细胞（单核细胞、pDC、中性粒细胞），而非T细胞。</li>
<li>RZV疫苗诱导的是<strong>精准的、自限性的CD4⁺ T细胞适应性免疫</strong>（ZEB2 D365=+2.91），其保护性不依赖于激活I型IFN通路（ISG15在疫苗中完全不动）。</li>
<li>两个参照系为后续实验章节提供了<strong>双重评价标准</strong>：<br>
<strong>标准1</strong> — 避免Disease Signature（第四章固有免疫评价）<br>
<strong>标准2</strong> — 获得Protection Signature（第五、六章T细胞应答评价）</li>
</ol>
</div>

<hr><p style="color:#999;text-align:center;font-size:.85em">
第三章报告 | HZ参照系数据：Vandoren et al. 2024 Supplementary Table S5 | 生成：2026-05-31 | <a href="https://github.com/cmjlabs/VZV-geo">github.com/cmjlabs/VZV-geo</a>
</p></body></html>'''

out = os.path.join(RES, "Chapter3_Report.html")
with open(out, 'w', encoding='utf-8') as f:
    f.write(html)
print(f"Written: {out} ({os.path.getsize(out)/1024:.0f} KB)")
