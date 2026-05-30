#!/usr/bin/env python3
"""Generate all cross-dataset comparison figures and tables."""
import os, sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HZ  = os.path.join(PROJ_ROOT, "results", "GSE242252")
RZV = os.path.join(PROJ_ROOT, "results", "GSE249632")
OUT = os.path.join(PROJ_ROOT, "results", "comparison")
os.makedirs(OUT, exist_ok=True)

# Load data
hz = pd.read_csv(os.path.join(HZ, "DE_HZ_annotated.csv"))
hz = hz[hz['symbol'].notna()]
rzv_raw = pd.read_csv(os.path.join(RZV, "logFC_matrix_all_timepoints.csv"))
c0 = rzv_raw.columns[0]
rzv_raw = rzv_raw.rename(columns={c0: 'gene_id'})
ens2sym = dict(zip(hz['ensembl_id_clean'].dropna(), hz['symbol'].dropna()))
rzv_raw['symbol'] = [ens2sym.get(str(x).split('.')[0], '') for x in rzv_raw['gene_id']]
rzv = rzv_raw[rzv_raw['symbol'] != ''].drop(columns=['gene_id']).groupby('symbol').mean()

hz_m = hz[['symbol','log2FoldChange','padj']].dropna().groupby('symbol').agg({'log2FoldChange':'mean','padj':'min'})
merged = hz_m.join(rzv, how='inner')
print(f"Common genes: {len(merged)}")

# Gene sets
hl_genes = ['ISG15','RSAD2','IFI44L','IFI44','IFI27','SERPING1','MX1','IFIT5','OASL',
            'GZMA','GNLY','IRF4','CTLA4','ZEB2','ICOS','TCF7','CCR7','CD38','TNFRSF9',
            'CXCR5','PDCD1','MKI67','TOP2A','MZB1','BATF2','PTTG1','CXCR3']
hl_genes = [g for g in hl_genes if g in merged.index]

# ── 1. Quadrant plots ────────────────────────────────────────────────────────
for tp, lb in [('D14','D14 vs D0'),('D74','D74 vs D0'),('D365','D365 vs D0')]:
    if tp not in merged.columns: continue
    df = merged[[tp,'log2FoldChange']].dropna()
    rho,p = spearmanr(df[tp], df['log2FoldChange'])
    fig,ax = plt.subplots(figsize=(10,9))
    ax.scatter(df['log2FoldChange'], df[tp], alpha=0.12, s=3, color='grey')
    q1=((df['log2FoldChange']>0)&(df[tp]>0)).sum()
    q2=((df['log2FoldChange']>0)&(df[tp]<0)).sum()
    q3=((df['log2FoldChange']<0)&(df[tp]<0)).sum()
    q4=((df['log2FoldChange']<0)&(df[tp]>0)).sum()
    xl=max(abs(df['log2FoldChange']).quantile(0.995),3)
    yl=max(abs(df[tp]).quantile(0.995),3)
    for lx,ly,lt in [(xl*.7,yl*.7,f'Q1\nBoth↑\n({q1})'),(-xl*.7,yl*.7,f'Q2\nHZ↑RZV↓\n({q2})'),
                      (-xl*.7,-yl*.7,f'Q3\nBoth↓\n({q3})'),(xl*.7,-yl*.7,f'Q4\nHZ↓RZV↑\n({q4})')]:
        ax.text(lx,ly,lt,ha='center',fontsize=10)
    for g in hl_genes:
        ax.annotate(g,(df.loc[g,'log2FoldChange'],df.loc[g,tp]),fontsize=7,fontweight='bold',color='black')
    ax.axhline(0,color='black',lw=0.5); ax.axvline(0,color='black',lw=0.5)
    ax.set_xlim(-xl,xl); ax.set_ylim(-yl,yl)
    ax.set_xlabel('HZ log2FC (acute vs convalescent)',fontsize=13)
    ax.set_ylabel(f'RZV log2FC ({lb})',fontsize=13)
    ax.set_title(f'HZ Disease vs RZV Vaccine: {lb}\nSpearman ρ={rho:.3f} (p={p:.2e}) | {len(df)} genes',fontsize=14)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"quadrant_HZ_vs_RZV_{tp}.png"), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Quadrant {tp}: ρ={rho:.3f}, p={p:.2e}, Q1={q1}, Q2={q2}, Q3={q3}, Q4={q4}")

# ── 2. Key gene comparison table ─────────────────────────────────────────────
kg = ['ISG15','RSAD2','IFI44L','IFI44','IFI27','SERPING1','SIGLEC1','MX1','IFIT5',
      'OASL','BATF2','GZMA','GZMB','PRF1','GNLY','IRF4','CTLA4','ZEB2','ICOS',
      'TCF7','CCR7','IL7R','CD38','TNFRSF9','CXCR5','CXCR3','PDCD1','MKI67','TOP2A','MZB1','PTTG1']
comp = merged.loc[merged.index.intersection(kg)].copy()
comp['HZ_sig'] = comp['padj'].apply(lambda x: '***' if x<0.001 else '**' if x<0.01 else '*' if x<0.05 else 'ns')
tps = [c for c in ['D14','D60','D74','D365'] if c in merged.columns]
comp[['log2FoldChange','padj','HZ_sig']+tps].round(3).to_csv(os.path.join(OUT, "key_gene_comparison.csv"))
print(f"Key gene table: {len(comp)} genes")

# ── 3. IFN comparison table ──────────────────────────────────────────────────
ifna = ['ISG15','RSAD2','IFI44L','IFI44','IFI27','IFI35','MX1','MX2',
        'OAS1','OAS2','OAS3','OASL','IFIT1','IFIT2','IFIT3','IFIT5',
        'IFITM1','IFITM2','IFITM3','IRF7','IRF9','STAT1','STAT2',
        'DDX58','DDX60','IFIH1','EIF2AK2','XAF1','SAMD9','SAMD9L',
        'USP18','ISG20','HERC5','HERC6','TRIM14','TRIM21','TRIM25',
        'TRIM26','UBE2L6','PLSCR1','RTP4','LGALS3BP']
ifng = ['IRF1','IRF8','STAT1','STAT3','STAT4','CXCL9','CXCL10','CXCL11',
        'GBP1','GBP2','GBP4','GBP5','CIITA','SOCS1','SOCS3',
        'TAP1','TAP2','PSMB8','PSMB9','PSMB10','ICAM1','B2M',
        'CCL2','CCL5','CD40','CD86','IL12RB1','IL15','IL18','CSF1','NOD1']
tcell = ['ICOS','CD38','IRF4','ZEB2','CTLA4','GZMA','GNLY','PRF1',
         'TCF7','CCR7','IL7R','TNFRSF9','CXCR5','PDCD1','BCL2L11']
rows = []
for genes,pathway in [(ifna,'Type I IFN'),(ifng,'Type II IFN hallmark'),(tcell,'T Cell Activation')]:
    for g in genes:
        row = {'gene':g,'pathway':pathway}
        hm = hz[hz['symbol']==g]
        row['HZ_log2FC'] = round(hm.iloc[0]['log2FoldChange'],3) if len(hm)>0 else np.nan
        row['HZ_padj'] = hm.iloc[0]['padj'] if len(hm)>0 else np.nan
        for tp in ['D14','D60','D74','D365']:
            row[f'RZV_{tp}'] = round(rzv.loc[g,tp],3) if g in rzv.index else np.nan
        rows.append(row)
pd.DataFrame(rows).to_csv(os.path.join(OUT, "IFN_comparison_table.csv"), index=False)
print(f"IFN table: {len(rows)} genes")

# ── 4. Th1/Th2/Tfh balance plot ──────────────────────────────────────────────
fig,axes = plt.subplots(1,3,figsize=(18,6))
for i,(tp,lb) in enumerate([('D14','D14 (1st dose)'),('D74','D74 (2nd dose)'),('D365','D365 (1 year)')]):
    ax=axes[i]
    th1_genes = ['TBX21','STAT1','STAT4','IL18R1','CXCR3','IL12RB2','TNF','IFNG','LTB','EOMES','RUNX3']
    th2_genes = ['GATA3','STAT6','CCR3','IL4R','MAF','CRLF2']
    tfh_genes = ['CXCR5','BCL6','PDCD1','ICOS','CD40LG','IL21','BTLA','SH2D1A','TOX2','ASCL2']
    vals = []
    labels = []
    for name,gs in [('Th1',th1_genes),('Th2',th2_genes),('Tfh',tfh_genes)]:
        found = [rzv.loc[g,tp] for g in gs if g in rzv.index]
        vals.append(np.mean(found) if found else 0)
        labels.append(name)
    colors = ['#ff7f0e','#1f77b4','#2ca02c']
    bars = ax.bar(labels, vals, color=colors, edgecolor='black', linewidth=0.5)
    ax.axhline(0,color='black',lw=0.5)
    for bar,val in zip(bars,vals):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.05 if val>=0 else bar.get_height()-0.15,
                f'{val:+.2f}', ha='center', fontweight='bold')
    ax.set_ylabel('Mean log2FC',fontsize=11)
    ax.set_title(f'RZV {lb}',fontsize=12,fontweight='bold')
fig.suptitle('RZV Vaccine: Th1/Th2/Tfh Balance in gE-Specific CD4+ T Cells',fontsize=14,fontweight='bold',y=1.02)
fig.tight_layout()
fig.savefig(os.path.join(OUT,"Th1_Th2_Tfh_balance.png"),dpi=150,bbox_inches='tight')
plt.close()
print("Th1/Th2/Tfh balance saved")

# ── 5. IFN summary bar chart ─────────────────────────────────────────────────
fig,ax = plt.subplots(figsize=(10,6))
t1_hz = [hz[hz['symbol']==g].iloc[0]['log2FoldChange'] for g in ifna if g in hz['symbol'].values]
t1_d14 = [rzv.loc[g,'D14'] for g in ifna if g in rzv.index]
t2_hz = [hz[hz['symbol']==g].iloc[0]['log2FoldChange'] for g in ifng if g in hz['symbol'].values]
t2_d14 = [rzv.loc[g,'D14'] for g in ifng if g in rzv.index]
tc_d14 = [rzv.loc[g,'D14'] for g in tcell if g in rzv.index]
cats = ['Type I IFN\n(HZ)','Type I IFN\n(RZV D14)','Type II IFN\n(HZ)','Type II IFN\n(RZV D14)','T Cell Act\n(RZV D14)']
means = [np.mean(t1_hz),np.mean(t1_d14),np.mean(t2_hz),np.mean(t2_d14),np.mean(tc_d14)]
colors = ['#E41A1C','#f4a582','#92c5de','#92c5de','#0571b0']
bars = ax.bar(range(len(cats)),means,color=colors,edgecolor='black',linewidth=0.5)
ax.axhline(0,color='black',lw=0.8)
ax.set_xticks(range(len(cats)))
ax.set_xticklabels(cats,fontsize=10)
ax.set_ylabel('Mean log2FC',fontsize=12)
ax.set_title('IFN Pathways and T Cell Activation: HZ vs RZV',fontsize=14,fontweight='bold')
for bar,val in zip(bars,means):
    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.02 if val>=0 else bar.get_height()-0.1,
            f'{val:+.2f}',ha='center',fontsize=10,fontweight='bold')
fig.tight_layout()
fig.savefig(os.path.join(OUT,"IFN_summary_barchart.png"),dpi=150,bbox_inches='tight')
plt.close()
print("IFN summary saved")
print("\nAll comparison figures regenerated.")
