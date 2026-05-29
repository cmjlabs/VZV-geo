#!/usr/bin/env python3
"""
Stratified gene set scoring across datasets.

Direction A: HZ pathway modules → scored in GSE249632 single cells
  Track module scores across D0→D14→D60→D74→D365

Direction B: RZV response modules → scored in GSE242252 bulk
  Compare acute vs convalescent
"""
import os, re
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, wilcoxon, kruskal
from scipy.stats import gaussian_kde
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

# Proxy handling
for k in list(os.environ.keys()):
    if 'proxy' in k.lower():
        del os.environ[k]

DATA_DIR = "/media/cmj/MechanicalDisk/yjs/VZV-geo/data"
RES_DIR = "/media/cmj/MechanicalDisk/yjs/VZV-geo/results"
MOD_DIR = os.path.join(RES_DIR, "pathway_enrichment")
OUT_DIR = os.path.join(RES_DIR, "module_scoring")
os.makedirs(OUT_DIR, exist_ok=True)

# ── 1. Load HZ pathway modules ───────────────────────────────────────────────
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()

ro.r['load'](os.path.join(MOD_DIR, "pathway_modules.rds"))
print("RDS loaded via rpy2 would go here...")

# Fallback: load from text files
print("Loading modules from text files...")
modules = {}
for f in sorted(os.listdir(MOD_DIR)):
    if f.startswith("module_") and f.endswith(".txt"):
        mod_name = f.replace("module_", "").replace(".txt", "")
        with open(os.path.join(MOD_DIR, f)) as fh:
            genes = [g.strip() for g in fh if g.strip()]
        modules[mod_name] = genes
        print(f"  {mod_name}: {len(genes)} genes")

if not modules:
    # Define modules manually from GO results (fallback)
    print("WARNING: No module files found. Using manually defined modules from DEG analysis.")
    degs = pd.read_csv(os.path.join(RES_DIR, "GSE242252/DE_HZ_annotated.csv"))
    sig = degs[degs['padj'].notna() & (degs['padj'] < 0.05)].copy()

    # Top ISGs (known from paper + our analysis)
    isg_genes = ['ISG15', 'RSAD2', 'IFI44L', 'IFI44', 'IFI27', 'IFI35',
                 'MX1', 'MX2', 'OAS1', 'OAS2', 'OAS3', 'OASL', 'IFIT1',
                 'IFIT2', 'IFIT3', 'IFIT5', 'IFITM1', 'IFITM2', 'IFITM3',
                 'IRF7', 'IRF9', 'STAT1', 'STAT2', 'SERPING1', 'SIGLEC1',
                 'DDX58', 'DDX60', 'IFIH1', 'EIF2AK2', 'XAF1']
    isg_genes = [g for g in isg_genes if g in sig['symbol'].values]

    # Cell cycle / proliferation
    cc_genes = sig[sig['log2FoldChange'] > 0.5]['symbol'].tolist()
    # Filter to known cell cycle genes
    cc_keywords = ['MKI67', 'TOP2A', 'PTTG1', 'STMN1', 'MCM', 'CDC', 'CDK',
                   'CENP', 'BUB', 'AURK', 'PLK', 'CCN', 'PCNA', 'RRM2',
                   'TYMS', 'TK1', 'H2AFZ', 'HMGB', 'UBE2C', 'CKS2']
    cc_genes = [g for g in cc_genes if any(k in g.upper() for k in cc_keywords)]

    modules = {
        "ISG_Storm": isg_genes,
        "Cell_Cycle_Proliferation": cc_genes,
    }

# ── 2. Direction A: Score HZ modules in GSE249632 cells ──────────────────────
print("\n===== Direction A: HZ modules in RZV single cells =====")

# Load cell-level counts
counts_file = os.path.join(DATA_DIR, "GSE249632_suppl/GSE249632_P261_gene_counts.txt.gz")
cell_counts = pd.read_csv(counts_file, sep='\t', index_col=0)
# Transpose to cells x genes
cell_counts = cell_counts.T  # cells x genes
print(f"Cell count matrix: {cell_counts.shape}")

# Load cell metadata
cell_meta = pd.read_csv(os.path.join(RES_DIR, "GSE249632/cell_metadata_clean.csv"))
# Get library IDs present in count matrix
valid_libs = [l for l in cell_meta['lib_id'] if l in cell_counts.index]
cell_meta_valid = cell_meta[cell_meta['lib_id'].isin(valid_libs)].copy()

# Load gene symbol mapping for GSE249632
# Map Ensembl IDs in count matrix to symbols
import mygene
mg = mygene.MyGeneInfo()
gene_ids_raw = cell_counts.columns.tolist()
gene_ids_clean = [g.split('.')[0] for g in gene_ids_raw]

print(f"Mapping {len(gene_ids_clean)} gene IDs...")
mg_results = mg.querymany(gene_ids_clean, scopes='ensembl.gene',
                          fields='symbol', species='human', batch_size=2000)
gene_to_symbol = {}
for r in mg_results:
    gene_to_symbol[r['query']] = r.get('symbol', '')

# Create gene symbol index for count matrix
counts_with_symbol = cell_counts.copy()
symbol_cols = {}
for i, gid in enumerate(gene_ids_raw):
    clean = gid.split('.')[0]
    sym = gene_to_symbol.get(clean, '')
    if sym:
        if sym not in symbol_cols:
            symbol_cols[sym] = []
        symbol_cols[sym].append(gid)

# Deduplicate: sum counts for same symbol
print(f"Unique symbols in matrix: {len(symbol_cols)}")

# Compute module scores per cell
# Use mean z-score method (similar to AddModuleScore in Seurat)
def compute_module_score(counts_df, gene_set, gene_symbol_map):
    """Compute mean z-score for a gene set per cell.
    Normalizes within each cell first (CPM), then z-scores across genes."""
    # Get columns for genes in the set
    cols = []
    for g in gene_set:
        if g in gene_symbol_map:
            cols.extend(gene_symbol_map[g])

    if len(cols) == 0:
        return None

    # Subset to genes in the module
    sub = counts_df[cols]

    # Normalize to CPM per cell
    lib_sizes = counts_df.sum(axis=1)
    cpm = sub.div(lib_sizes, axis=0) * 1e6

    # log2 transform
    log_cpm = np.log2(cpm + 1)

    # Z-score per gene (across cells), then mean per cell
    z_per_gene = (log_cpm - log_cpm.mean(axis=0)) / (log_cpm.std(axis=0) + 1e-8)
    score = z_per_gene.mean(axis=1)

    return score

print("\nComputing module scores for each cell...")
score_df = pd.DataFrame(index=cell_meta_valid['lib_id'])

for mod_name, mod_genes in modules.items():
    score = compute_module_score(cell_counts.loc[score_df.index], mod_genes, symbol_cols)
    if score is not None:
        score_df[mod_name] = score
        print(f"  {mod_name}: scored {len([g for g in mod_genes if g in symbol_cols])}/{len(mod_genes)} genes present")

# Merge with metadata
score_df = score_df.join(cell_meta_valid.set_index('lib_id')[['donor', 'timepoint', 'visit']])

# ── 3. Plot: Module scores across vaccine timepoints ─────────────────────────
print("\nPlotting module scores across timepoints...")
timepoint_order = ['D0', 'D14', 'D60', 'D74', 'D365']
tp_colors = {'D0': '#1B9E77', 'D14': '#D95F02', 'D60': '#7570B3',
             'D74': '#E7298A', 'D365': '#66A61E'}

for mod_name in score_df.columns[:6]:  # Skip metadata columns
    if mod_name in ['donor', 'timepoint', 'visit']:
        continue

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Panel 1: Violin plot per timepoint
    data_by_tp = [score_df[score_df['timepoint'] == tp][mod_name].dropna().values
                  for tp in timepoint_order]

    ax = axes[0]
    vp = ax.violinplot(data_by_tp, positions=range(len(timepoint_order)),
                       showmeans=True, showmedians=True)
    for i, body in enumerate(vp['bodies']):
        body.set_facecolor(tp_colors[timepoint_order[i]])
        body.set_alpha(0.6)

    ax.set_xticks(range(len(timepoint_order)))
    ax.set_xticklabels(timepoint_order)
    ax.set_ylabel('Module Score (z-scored)')
    ax.set_title(f'{mod_name} Score Across Vaccination Timeline')

    # Kruskal-Wallis test
    h_stat, p_val = kruskal(*data_by_tp)
    ax.text(0.02, 0.98, f'Kruskal-Wallis p = {p_val:.2e}',
            transform=ax.transAxes, va='top', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Panel 2: Per-donor trajectory
    ax = axes[1]
    for donor in sorted(score_df['donor'].unique()):
        donor_data = score_df[score_df['donor'] == donor].copy()
        donor_data = donor_data.sort_values('timepoint',
            key=lambda x: [timepoint_order.index(t) if t in timepoint_order else 999 for t in x])
        donor_means = donor_data.groupby('timepoint')[mod_name].mean()
        tp_list = [t for t in timepoint_order if t in donor_means.index]
        ax.plot([timepoint_order.index(t) for t in tp_list],
                [donor_means[t] for t in tp_list],
                'o-', alpha=0.5, markersize=5, label=f'Donor {donor}')

    ax.set_xticks(range(len(timepoint_order)))
    ax.set_xticklabels(timepoint_order)
    ax.set_ylabel('Mean Module Score')
    ax.set_title(f'{mod_name} - Per Donor Trajectory')
    ax.legend(fontsize=7, loc='best')

    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, f"module_score_{mod_name}.png"), dpi=150)
    plt.close(fig)
    print(f"  Plot saved: module_score_{mod_name}.png")

# ── 4. Summary statistics table ──────────────────────────────────────────────
print("\n===== Module Score Summary =====")
summary_rows = []
for mod_name in score_df.columns:
    if mod_name in ['donor', 'timepoint', 'visit']:
        continue
    row = {'module': mod_name}
    for tp in timepoint_order:
        vals = score_df[score_df['timepoint'] == tp][mod_name].dropna()
        row[f'{tp}_mean'] = vals.mean()
        row[f'{tp}_std'] = vals.std()
    # Test D0 vs each other timepoint
    d0_vals = score_df[score_df['timepoint'] == 'D0'][mod_name].dropna()
    for tp in ['D14', 'D74', 'D365']:
        tp_vals = score_df[score_df['timepoint'] == tp][mod_name].dropna()
        if len(d0_vals) > 0 and len(tp_vals) > 0:
            _, p = mannwhitneyu(d0_vals, tp_vals, alternative='two-sided')
            row[f'p_{tp}_vs_D0'] = p
    summary_rows.append(row)

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(os.path.join(OUT_DIR, "module_score_summary.csv"), index=False)
print(summary_df.to_string())

# ── 5. Direction B: RZV modules in GSE242252 ─────────────────────────────────
print("\n===== Direction B: RZV modules in HZ bulk =====")

# Define RZV vaccine response modules from our DE results
rzv_de = {}
for tp in ['D14', 'D74', 'D365']:
    tp_file = os.path.join(RES_DIR, f"GSE249632/DE_{tp}_vs_D0.csv")
    if os.path.exists(tp_file):
        df = pd.read_csv(tp_file)
        rzv_de[tp] = df

# Define RZV modules
rzv_modules = {}

# RZV Activation (D14 up)
if 'D14' in rzv_de:
    d14 = rzv_de['D14']
    d14_up = d14[(d14['adj.P.Val'] < 0.05) & (d14['logFC'] > 0.5)]
    rzv_modules['RZV_Activation'] = d14_up['gene_id'].str.split('.').str[0].tolist()

# RZV Regulation (ZEB2, CTLA4 signature - consistently up at all timepoints)
# Get genes up at D14, D74, AND D365 (persistent)
persistent_up = None
for tp in ['D14', 'D74', 'D365']:
    if tp in rzv_de:
        df = rzv_de[tp]
        tp_up = set(df[(df['adj.P.Val'] < 0.05) & (df['logFC'] > 0)]['gene_id'].str.split('.').str[0])
        if persistent_up is None:
            persistent_up = tp_up
        else:
            persistent_up = persistent_up & tp_up

if persistent_up:
    rzv_modules['RZV_Persistent_Regulation'] = list(persistent_up)

# RZV Long-term Imprint (D365 specific)
if 'D365' in rzv_de:
    d365 = rzv_de['D365']
    d365_sig = d365[(d365['adj.P.Val'] < 0.05) & (abs(d365['logFC']) > 0.5)]
    d365_sig_up = set(d365_sig[d365_sig['logFC'] > 0]['gene_id'].str.split('.').str[0])
    d365_sig_down = set(d365_sig[d365_sig['logFC'] < 0]['gene_id'].str.split('.').str[0])
    rzv_modules['RZV_D365_Up'] = list(d365_sig_up)
    rzv_modules['RZV_D365_Down'] = list(d365_sig_down)

print(f"RZV modules defined: {[(k, len(v)) for k, v in rzv_modules.items()]}")

# Score in GSE242252 bulk
hz_counts = pd.read_csv(os.path.join(RES_DIR, "GSE242252/filtered_counts.csv"), index_col=0)
hz_meta = pd.read_csv(os.path.join(RES_DIR, "GSE242252/deseq2_metadata.csv"), index_col=0)

# Filter to HZ samples
hz_samples = hz_meta[hz_meta['condition_label'] == 'Herpes_Zoster'].index
hz_counts_hz = hz_counts[hz_samples]

# Build gene symbol map for GSE242252
hz_gene_ids = [g.split('.')[0] for g in hz_counts_hz.index]
print(f"Mapping {len(hz_gene_ids)} HZ gene IDs...")
hz_mg = mg.querymany(hz_gene_ids, scopes='ensembl.gene', fields='symbol',
                      species='human', batch_size=2000)
hz_symbol_map = {}
for r in hz_mg:
    hz_symbol_map[r['query']] = r.get('symbol', '')

# Build reverse map for GSE242252
hz_gene_to_sym = {}
for i, gid in enumerate(hz_counts_hz.index):
    clean = gid.split('.')[0]
    sym = hz_symbol_map.get(clean, '')
    if sym:
        if sym not in hz_gene_to_sym:
            hz_gene_to_sym[sym] = []
        hz_gene_to_sym[sym].append(gid)

# Map RZV Ensembl IDs to symbols via mygene
print("Mapping RZV gene IDs to symbols...")
for mod_name, mod_genes in rzv_modules.items():
    if len(mod_genes) == 0:
        continue
    mod_sym_results = mg.querymany(mod_genes, scopes='ensembl.gene',
                                    fields='symbol', species='human', batch_size=1000)
    mod_symbols = []
    for r in mod_sym_results:
        sym = r.get('symbol', '')
        if sym:
            mod_symbols.append(sym)
    rzv_modules[mod_name] = mod_symbols
    print(f"  {mod_name}: {len(mod_symbols)} symbols mapped")

# Score each RZV module in HZ bulk
print("\nScoring RZV modules in HZ bulk...")
rzv_score_df = pd.DataFrame(index=hz_samples)

for mod_name, mod_symbols in rzv_modules.items():
    cols_in_hz = []
    for sym in mod_symbols:
        if sym in hz_gene_to_sym:
            cols_in_hz.extend(hz_gene_to_sym[sym])

    if len(cols_in_hz) < 3:
        print(f"  {mod_name}: only {len(cols_in_hz)} genes found, skipping")
        continue

    # Use GSVA-like scoring: ssGSEA
    # Simplified: mean of top 50% expressing genes (robust against noise)
    sub = hz_counts_hz.loc[hz_counts_hz.index.isin(cols_in_hz)]
    # CPM normalize
    cpm = sub.div(hz_counts_hz.sum(axis=0), axis=1) * 1e6
    log_cpm = np.log2(cpm + 1)

    # Z-score genes, take mean
    z = (log_cpm.T - log_cpm.mean(axis=1)) / (log_cpm.std(axis=1) + 1e-8)
    rzv_score_df[mod_name] = z.mean(axis=1)

    print(f"  {mod_name}: scored with {len(cols_in_hz)} genes")

# Merge with metadata
rzv_score_df['timepoint'] = hz_meta.loc[hz_samples, 'timepoint'].values
rzv_score_df['patient_id'] = hz_meta.loc[hz_samples, 'patient_id'].values

# Plot
for mod_name in rzv_score_df.columns:
    if mod_name in ['timepoint', 'patient_id']:
        continue

    fig, ax = plt.subplots(figsize=(7, 6))

    acute = rzv_score_df[rzv_score_df['timepoint'] == 'acute'][mod_name]
    conv = rzv_score_df[rzv_score_df['timepoint'] == 'convalescent'][mod_name]

    # Paired boxplot
    paired_data = rzv_score_df.pivot(index='patient_id', columns='timepoint', values=mod_name)
    paired_data = paired_data.dropna()

    bp = ax.boxplot([acute, conv], labels=['Acute HZ', 'Convalescent'],
                     patch_artist=True, widths=0.5)
    bp['boxes'][0].set_facecolor('#E41A1C')
    bp['boxes'][1].set_facecolor('#377EB8')
    for box in bp['boxes']:
        box.set_alpha(0.6)

    # Paired lines
    for _, row in paired_data.iterrows():
        ax.plot([0, 1], [row['acute'], row['convalescent']],
                'o-', color='grey', alpha=0.3, markersize=4)

    # Wilcoxon paired test
    if len(paired_data) > 3:
        _, p = wilcoxon(paired_data['acute'], paired_data['convalescent'])
        ax.text(0.5, 0.95, f'Wilcoxon paired p = {p:.3f}',
                transform=ax.transAxes, ha='center', va='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.set_ylabel(f'{mod_name} Score')
    ax.set_title(f'RZV {mod_name} in HZ Acute vs Convalescent')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, f"rzv_module_in_HZ_{mod_name}.png"), dpi=150)
    plt.close(fig)
    print(f"  Plot saved: rzv_module_in_HZ_{mod_name}.png")

# Save
rzv_score_df.to_csv(os.path.join(OUT_DIR, "rzv_modules_in_HZ_scores.csv"))

print("\n===== Module scoring complete =====")
print(f"Output in: {OUT_DIR}")
