#!/usr/bin/env python3
"""Generate Chapter 3 narrative figures for the thesis."""
import os, pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import spearmanr

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(PROJ, "results")
HZ_DIR  = os.path.join(RES, "GSE242252")
RZV_DIR = os.path.join(RES, "GSE249632")
OUT_DIR = os.path.join(RES, "chapter3_figures")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Load data ─────────────────────────────────────────────────────────────────
hz = pd.read_csv(os.path.join(HZ_DIR, "DE_HZ_annotated.csv"))
hz = hz[hz['symbol'].notna()]

rzv_de = {}
for tp in ['D14','D60','D74','D365']:
    df = pd.read_csv(os.path.join(RZV_DIR, f"DE_{tp}_vs_D0.csv"))
    df['ensembl_clean'] = df['gene_id'].str.split('.').str[0]
    rzv_de[tp] = df

rzv_lfc = pd.read_csv(os.path.join(RZV_DIR, "logFC_matrix_all_timepoints.csv"))
c0 = rzv_lfc.columns[0]
rzv_lfc = rzv_lfc.rename(columns={c0: 'gene_id'})
ens2sym = dict(zip(hz['ensembl_id_clean'].dropna(), hz['symbol'].dropna()))
rzv_lfc['symbol'] = [ens2sym.get(str(x).split('.')[0], '') for x in rzv_lfc['gene_id']]
rzv = rzv_lfc[rzv_lfc['symbol'] != ''].drop(columns=['gene_id']).groupby('symbol').mean()

sns_set = False
try:
    import seaborn as sns
    sns.set_style("whitegrid")
    sns_set = True
except ImportError:
    pass

# ══════════════════════════════════════════════════════════════════════════════
# Result 2: GO enrichment bubble plot
# ══════════════════════════════════════════════════════════════════════════════
print("Generating Result 2: GO bubble plot...")
go_up = pd.read_csv(os.path.join(RES, "pathway_enrichment", "GO_BP_upregulated.csv"))
go_up['n_genes'] = go_up['Genes'].apply(lambda x: len(x.split(';')))
go_up = go_up.sort_values('Adjusted P-value')
top_go = go_up[go_up['n_genes'] >= 3].head(8)

fig, ax = plt.subplots(figsize=(10, 5))
colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(top_go)))
for i, (_, row) in enumerate(top_go.iterrows()):
    ax.scatter(-np.log10(row['Adjusted P-value']), row['Term'][:70],
               s=row['n_genes']*30, color=colors[i], edgecolor='black', linewidth=0.5,
               zorder=5)
ax.axvline(-np.log10(0.05), color='red', linestyle='--', alpha=0.5, label='p=0.05')
ax.set_xlabel('-log10(Adjusted P-value)', fontsize=12)
ax.set_title('GO Biological Process: HZ Acute Up-regulated Genes', fontsize=13, fontweight='bold')
# Legend for bubble size
for n in [3, 5, 8]:
    ax.scatter([], [], s=n*30, color='#808080', edgecolor='black', linewidth=0.5, label=f'{n} genes')
ax.legend(title='Gene count', loc='lower right', fontsize=8, title_fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "GO_bubble_plot.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  GO bubble plot saved.")

# ══════════════════════════════════════════════════════════════════════════════
# Result 3: DEG timeline bar chart
# ══════════════════════════════════════════════════════════════════════════════
print("Generating Result 3: DEG timeline bar chart...")
tps = ['D14','D60','D74','D365']
up_counts = []
down_counts = []
for tp in tps:
    df = rzv_de[tp]
    up = (df['adj.P.Val'] < 0.05) & (df['logFC'] > 0)
    down = (df['adj.P.Val'] < 0.05) & (df['logFC'] < 0)
    up_counts.append(up.sum())
    down_counts.append(down.sum())

fig, ax = plt.subplots(figsize=(8, 5))
x = np.arange(len(tps))
width = 0.35
bars_up = ax.bar(x - width/2, up_counts, width, label='Up-regulated', color='#E41A1C', edgecolor='black', linewidth=0.5)
bars_down = ax.bar(x + width/2, down_counts, width, label='Down-regulated', color='#377EB8', edgecolor='black', linewidth=0.5)
for bar, val in zip(bars_up, up_counts):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 3, str(val), ha='center', fontsize=11, fontweight='bold')
for bar, val in zip(bars_down, down_counts):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 3, str(val), ha='center', fontsize=11, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(['D14 vs D0\n(1st dose)', 'D60 vs D0\n(pre-2nd)', 'D74 vs D0\n(2nd dose)', 'D365 vs D0\n(1 year)'], fontsize=10)
ax.set_ylabel('Number of DEGs (FDR<0.05)', fontsize=12)
ax.set_title('RZV Vaccine: CD4+ T Cell Transcriptional Response', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.set_ylim(0, max(up_counts + down_counts) * 1.2)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "DEG_timeline_barchart.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  DEG timeline saved.")

# ══════════════════════════════════════════════════════════════════════════════
# Result 4: 3-gene trajectory — ZEB2, CTLA4, ISG15
# ══════════════════════════════════════════════════════════════════════════════
print("Generating Result 4: 3-gene trajectory...")
# Get HZ values
hz_vals = {}
for g in ['ZEB2','CTLA4','ISG15']:
    match = hz[hz['symbol'] == g]
    if len(match) > 0:
        hz_vals[g] = {
            'lfc': match.iloc[0]['log2FoldChange'],
            'padj': match.iloc[0]['padj'],
            'sig': '***' if match.iloc[0]['padj'] < 0.001 else ('**' if match.iloc[0]['padj'] < 0.01 else ('*' if match.iloc[0]['padj'] < 0.05 else 'ns'))
        }

# Get RZV trajectory values
rzv_traj = {}
for g in ['ZEB2','CTLA4','ISG15']:
    if g in rzv.index:
        rzv_traj[g] = {tp: rzv.loc[g, tp] for tp in ['D14','D60','D74','D365'] if tp in rzv.columns}

fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
gene_info = [
    ('ZEB2', '#2ca02c', 'T cell differentiation\nreprogramming'),
    ('CTLA4', '#ff7f0e', 'Immune self-limitation\n(checkpoint)'),
    ('ISG15', '#d62728', 'Type I IFN response\n(innate inflammation)'),
]
for ax, (gene, color, label) in zip(axes, gene_info):
    # RZV trajectory
    if gene in rzv_traj:
        tp_order = ['D14','D60','D74','D365']
        rzv_ys = [rzv_traj[gene].get(tp, np.nan) for tp in tp_order]
        ax.plot(range(4), rzv_ys, 'o-', color=color, linewidth=2.5, markersize=10, label='RZV vaccine')
        # Add value labels
        for i, (tp, y) in enumerate(zip(tp_order, rzv_ys)):
            if not np.isnan(y):
                ax.annotate(f'{y:+.1f}', (i, y), textcoords="offset points", xytext=(0, 12),
                           ha='center', fontsize=9, fontweight='bold', color=color)
    # HZ value as horizontal line
    if gene in hz_vals:
        hz_y = hz_vals[gene]['lfc']
        sig = hz_vals[gene]['sig']
        ax.axhline(y=hz_y, color='#808080', linestyle='--', linewidth=1.5, alpha=0.7,
                   label=f'HZ acute: {hz_y:+.2f} {sig}')
        ax.fill_between([-0.3, 3.3], hz_y-0.1, hz_y+0.1, alpha=0.08, color='#808080')
    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.set_xticks(range(4))
    ax.set_xticklabels(['D14', 'D60', 'D74', 'D365'], fontsize=11)
    ax.set_ylabel('log2 Fold Change vs D0', fontsize=11)
    ax.set_title(f'{gene}\n({label})', fontsize=12, fontweight='bold', color=color)
    ax.legend(fontsize=8, loc='best')

fig.suptitle('Key Genes: HZ Disease vs RZV Vaccine Response', fontsize=15, fontweight='bold', y=1.02)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "three_gene_trajectory.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  3-gene trajectory saved.")

# ══════════════════════════════════════════════════════════════════════════════
# Result 1: Updated volcano with only 8 narrative genes
# ══════════════════════════════════════════════════════════════════════════════
print("Generating Result 1: Narrative volcano plot...")
res_df = hz[hz['padj'].notna()].copy()
res_df['sig'] = 'NS'
res_df.loc[(res_df['padj'] < 0.05) & (res_df['log2FoldChange'] > 1), 'sig'] = 'Up'
res_df.loc[(res_df['padj'] < 0.05) & (res_df['log2FoldChange'] < -1), 'sig'] = 'Down'
n_up = (res_df['padj'] < 0.05) & (res_df['log2FoldChange'] > 0)
n_down = (res_df['padj'] < 0.05) & (res_df['log2FoldChange'] < 0)
narrative_genes = ['ISG15','IFI44L','RSAD2','IFIT5','TOP2A','PTTG1','MZB1','SERPING1']

fig, ax = plt.subplots(figsize=(10, 8))
for s, c, alpha, size in [('NS', '#d9d9d9', 0.3, 0.5), ('Down', '#377EB8', 0.5, 0.8), ('Up', '#E41A1C', 0.5, 0.8)]:
    subset = res_df[res_df['sig'] == s]
    ax.scatter(subset['log2FoldChange'], -np.log10(subset['padj']), c=c, alpha=alpha, s=size, label=s)

x_max = max(abs(res_df['log2FoldChange']).max(), 3) * 1.05
ax.set_xlim(-x_max, x_max)
ax.axvline(0, color='#808080', linewidth=0.5)
ax.axvline(-1, color='#808080', linestyle='--', alpha=0.3)
ax.axvline(1, color='#808080', linestyle='--', alpha=0.3)
ax.axhline(-np.log10(0.05), color='#808080', linestyle='--', alpha=0.3)

for g in narrative_genes:
    match = res_df[res_df['symbol'] == g]
    if len(match) > 0:
        x, y = match.iloc[0]['log2FoldChange'], -np.log10(match.iloc[0]['padj'])
        p = match.iloc[0]['padj']
        sig = '***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else ''
        ax.annotate(f'{g} {sig}', (x, y), fontsize=10, fontweight='bold', color='black',
                   arrowprops=dict(arrowstyle='->', color='black', lw=0.8),
                   textcoords="offset points", xytext=(8, 8))

ax.set_xlabel('log2 Fold Change (acute vs convalescent)', fontsize=13)
ax.set_ylabel('-log10(adjusted p-value)', fontsize=13)
ax.set_title(f'GSE242252: HZ Acute vs Convalescent\n{n_up.sum()} up, {n_down.sum()} down (FDR<0.05) | 23 paired patients',
             fontsize=14, fontweight='bold')
ax.legend(fontsize=9, loc='upper right')
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "volcano_narrative.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  Narrative volcano saved.")

# ══════════════════════════════════════════════════════════════════════════════
# Result 5: Quadrant plot (D14) with narrative annotations
# ══════════════════════════════════════════════════════════════════════════════
print("Generating Result 5: Narrative quadrant plot...")
hz_m = hz[['symbol','log2FoldChange','padj']].dropna().groupby('symbol').agg({'log2FoldChange':'mean','padj':'min'})
merged = hz_m.join(rzv, how='inner')
tp = 'D14'
df = merged[[tp,'log2FoldChange']].dropna()
rho, p = spearmanr(df[tp], df['log2FoldChange'])

q2_genes = ['ISG15','RSAD2','IFI44L','IFIT5','SERPING1']  # HZ↑, RZV↓
q4_genes = ['ZEB2','CTLA4','ICOS']  # HZ↓/ns, RZV↑
q1_genes = ['TOP2A','CD38']  # Both↑

fig, ax = plt.subplots(figsize=(10, 9))
ax.scatter(df['log2FoldChange'], df[tp], alpha=0.1, s=2, color='#808080')
q1 = ((df['log2FoldChange']>0)&(df[tp]>0)).sum()
q2 = ((df['log2FoldChange']>0)&(df[tp]<0)).sum()
q3 = ((df['log2FoldChange']<0)&(df[tp]<0)).sum()
q4 = ((df['log2FoldChange']<0)&(df[tp]>0)).sum()
xl = max(abs(df['log2FoldChange']).quantile(0.995), 3)
yl = max(abs(df[tp]).quantile(0.995), 3)

# Quadrant annotations with biological meaning
for lx, ly, txt, color in [
    (xl*0.7, yl*0.7, 'Q1: Shared proliferation\nTOP2A, CD38', '#7f7f7f'),
    (-xl*0.7, yl*0.7, 'Q2: HZ-specific\nInnate inflammation\nISG15, RSAD2,\nIFI44L, IFIT5,\nSERPING1', '#d62728'),
    (-xl*0.7, -yl*0.7, 'Q3: Both↓', '#7f7f7f'),
    (xl*0.7, -yl*0.7, 'Q4: RZV-specific\nAdaptive regulation\nZEB2, CTLA4, ICOS', '#2ca02c'),
]:
    ax.text(lx, ly, txt, ha='center', fontsize=9, color=color, fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# Color-coded gene labels
for g in q2_genes:
    if g in df.index:
        ax.annotate(g, (df.loc[g,'log2FoldChange'], df.loc[g,tp]),
                   fontsize=8, fontweight='bold', color='#d62728')
for g in q4_genes:
    if g in df.index:
        ax.annotate(g, (df.loc[g,'log2FoldChange'], df.loc[g,tp]),
                   fontsize=8, fontweight='bold', color='#2ca02c')
for g in q1_genes:
    if g in df.index:
        ax.annotate(g, (df.loc[g,'log2FoldChange'], df.loc[g,tp]),
                   fontsize=8, fontweight='bold', color='#7f7f7f')

ax.axhline(0, color='black', lw=0.5); ax.axvline(0, color='black', lw=0.5)
ax.set_xlim(-xl, xl); ax.set_ylim(-yl, yl)
ax.set_xlabel('HZ log2FC (acute vs convalescent)', fontsize=13)
ax.set_ylabel('RZV log2FC (D14 vs D0)', fontsize=13)
ax.set_title(f'HZ Disease vs RZV Vaccine: Global Transcriptome Comparison\nSpearman ρ = {rho:.3f} | {len(df)} genes | Q1={q1}, Q2={q2}, Q3={q3}, Q4={q4}',
             fontsize=13, fontweight='bold')
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "quadrant_narrative.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  Narrative quadrant saved.")

print(f"\nAll Chapter 3 figures saved to {OUT_DIR}")
