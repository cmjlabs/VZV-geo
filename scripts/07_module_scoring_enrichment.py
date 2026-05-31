#!/usr/bin/env python3
"""
Complete pathway-stratified cross-dataset gene set scoring.
1. GO/KEGG enrichment on HZ DEGs → define functional modules
2. Direction A: Score HZ modules in GSE249632 cells across vaccine timeline
3. Direction B: Score RZV modules in GSE242252 bulk (acute vs convalescent)
"""
import re
import os, sys, json
for k in list(os.environ.keys()):
    if 'proxy' in k.lower():
        del os.environ[k]

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, wilcoxon, kruskal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("whitegrid")

# ── Paths ────────────────────────────────────────────────────────────────────
DATA_DIR = "/media/cmj/MechanicalDisk/yjs/VZV-geo/data"
RES_DIR = "/media/cmj/MechanicalDisk/yjs/VZV-geo/results"
MOD_DIR = os.path.join(RES_DIR, "pathway_enrichment")
OUT_DIR = os.path.join(RES_DIR, "module_scoring")
os.makedirs(MOD_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

# ── 1. Load HZ DEGs ─────────────────────────────────────────────────────────
print("=" * 70)
print("STEP 1: Load HZ DEGs and run pathway enrichment")
print("=" * 70)

hz_degs = pd.read_csv(os.path.join(RES_DIR, "GSE242252/DE_HZ_annotated.csv"))
hz_degs['ensembl_clean'] = hz_degs['gene_id'].str.split('.').str[0]

sig = hz_degs[hz_degs['padj'].notna() & (hz_degs['padj'] < 0.05)].copy()
# Use ALL FDR<0.05 genes for enrichment (no LFC cutoff — more genes = better enrichment power)
sig_up = sig[sig['log2FoldChange'] > 0]
sig_down = sig[sig['log2FoldChange'] < 0]

print(f"Total DEGs: {len(sig)} (up: {len(sig_up)}, down: {len(sig_down)})")

# ── 2. GO enrichment with gseapy ────────────────────────────────────────────
print("\nRunning GO BP enrichment...")
import gseapy as gp

# Use gene symbols for enrichment
bg_genes = hz_degs['symbol'].dropna().unique().tolist()
up_genes = sig_up['symbol'].dropna().unique().tolist()
down_genes = sig_down['symbol'].dropna().unique().tolist()

# Enrichr GO BP
try:
    enr_up = gp.enrichr(
        gene_list=up_genes,
        gene_sets=['GO_Biological_Process_2023'],
        organism='human',
        outdir=None,
        no_plot=True
    )
    enr_up_res = enr_up.results
    enr_up_res.to_csv(os.path.join(MOD_DIR, "GO_BP_upregulated.csv"), index=False)
    print(f"GO terms (up): {len(enr_up_res)}")
    print("\nTop 15 GO BP (up):")
    for _, row in enr_up_res.head(15).iterrows():
        print(f"  {row['Term'][:70]:70s} p={row['Adjusted P-value']:.2e} genes={len(row['Genes'].split(';'))}")
except Exception as e:
    print(f"GO enrichment error: {e}")
    enr_up_res = pd.DataFrame()

try:
    enr_down = gp.enrichr(
        gene_list=down_genes,
        gene_sets=['GO_Biological_Process_2023'],
        organism='human',
        outdir=None,
        no_plot=True
    )
    enr_down_res = enr_down.results
    enr_down_res.to_csv(os.path.join(MOD_DIR, "GO_BP_downregulated.csv"), index=False)
    print(f"\nGO terms (down): {len(enr_down_res)}")
    print("\nTop 10 GO BP (down):")
    for _, row in enr_down_res.head(10).iterrows():
        print(f"  {row['Term'][:70]:70s} p={row['Adjusted P-value']:.2e}")
except Exception as e:
    print(f"GO enrichment error: {e}")
    enr_down_res = pd.DataFrame()

# ── 3. Define functional gene modules ────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 2: Define functional gene modules from GO results")
print("=" * 70)

def merge_overlapping_modules(modules, overlap_threshold=0.5):
    """Merge modules whose gene sets overlap by > overlap_threshold (Jaccard)."""
    keys = list(modules.keys())
    merged = {}
    skip = set()
    for i, k1 in enumerate(keys):
        if k1 in skip:
            continue
        g1 = set(modules[k1]['genes'])
        mt = modules[k1]['terms'].copy()
        for j, k2 in enumerate(keys):
            if j <= i or k2 in skip:
                continue
            g2 = set(modules[k2]['genes'])
            jac = len(g1 & g2) / len(g1 | g2) if len(g1 | g2) > 0 else 0
            if jac > overlap_threshold:
                g1 |= g2
                mt.extend(modules[k2]['terms'])
                skip.add(k2)
        merged[k1] = {'name': modules[k1]['name'], 'genes': sorted(g1),
                      'n_terms': len(mt), 'terms': mt}
    return merged

def top_go_modules(enr_df, n_top=10, min_genes=3, direction='up',
                   gene_universe=None, max_background_dilution=0.3):
    """Extract gene modules from top N significant GO terms (data-driven).
    Skips terms that are diluted in the gene universe (too few of their genes in our DEG list).
    """
    if enr_df is None or len(enr_df) == 0:
        return {}
    df = enr_df.sort_values('Adjusted P-value').copy()

    # Print all candidates for transparency
    print(f"\n  Top {n_top} {direction}-regulated GO terms (by adj p-value):")
    for i, (_, row) in enumerate(df.head(n_top).iterrows()):
        n_genes = len(row['Genes'].split(';'))
        print(f"    {i+1:2d}. [{row['Adjusted P-value']:.1e}] {row['Term'][:60]:60s} ({n_genes} genes)")

    modules = {}
    taken = set()
    skipped_generic = 0
    for i, row in df.iterrows():
        if len(modules) >= n_top:
            break
        term = row['Term']
        genes = [g.strip() for g in row['Genes'].split(';')]

        # Skip if too few genes
        if len(genes) < min_genes:
            continue

        # Skip overly generic terms (common GO traps)
        generic_patterns = [
            'positive regulation of transcription',
            'negative regulation of transcription',
            'regulation of transcription',
            'positive regulation of gene expression',
            'regulation of gene expression',
            'positive regulation of RNA',
            'negative regulation of RNA',
            'regulation of RNA',
            'positive regulation of cellular biosynthetic',
            'positive regulation of nitrogen compound',
            'positive regulation of macromolecule',
            'cellular response to oxygen-containing',
            'positive regulation of nucleobase',
        ]
        term_lower = term.lower()
        if any(p in term_lower for p in generic_patterns):
            skipped_generic += 1
            continue

        # Skip if most genes already covered
        new_genes = set(genes) - taken
        if len(new_genes) < max(3, len(genes) * 0.4):
            continue

        key = re.sub(r'[^a-zA-Z0-9_]', '_', term)[:60]
        key = re.sub(r'_+', '_', key).strip('_')
        if key in modules:
            key = f"{key}_{i}"
        modules[key] = {'name': term, 'genes': sorted(genes),
                        'n_terms': 1, 'terms': [term]}
        taken.update(genes)

    if skipped_generic:
        print(f"  Skipped {skipped_generic} overly generic GO terms ({direction})")

    # If too few modules retained, supplement with top terms by p-value (no generic filter)
    if len(modules) < 3:
        print(f"  Only {len(modules)} modules after filtering — supplementing with top terms")
        for i, row in df.iterrows():
            if len(modules) >= max(n_top, 5):
                break
            term = row['Term']
            genes = [g.strip() for g in row['Genes'].split(';')]
            if len(genes) < min_genes:
                continue
            key = re.sub(r'[^a-zA-Z0-9_]', '_', term)[:60]
            key = re.sub(r'_+', '_', key).strip('_')
            if key not in modules:
                modules[key] = {'name': term, 'genes': sorted(genes),
                                'n_terms': 1, 'terms': [term]}

    modules = merge_overlapping_modules(modules)
    print(f"  {direction}-regulated: {len(modules)} modules retained")
    return modules

# Data-driven module definition: top N GO terms by significance
TOP_N_UP, TOP_N_DOWN = 12, 8
MIN_GENES = 3  # initial per-term minimum, will filter < 5 after merging

modules_up = top_go_modules(enr_up_res, TOP_N_UP, MIN_GENES, 'up') if len(enr_up_res) > 0 else {}
modules_down = top_go_modules(enr_down_res, TOP_N_DOWN, MIN_GENES, 'down') if len(enr_down_res) > 0 else {}

modules_up['HZ_Disease_All_Up'] = {
    'name': 'HZ Disease Signature (All Up-regulated DEGs)',
    'genes': sorted(up_genes),
    'n_terms': 'N/A',
    'terms': [f'All {len(up_genes)} up-regulated DEGs (FDR<0.05, LFC>0.5)']
}
modules_down['HZ_Disease_All_Down'] = {
    'name': 'HZ Disease Signature (All Down-regulated DEGs)',
    'genes': sorted(down_genes),
    'n_terms': 'N/A',
    'terms': [f'All {len(down_genes)} down-regulated DEGs (FDR<0.05, LFC<-0.5)']
}

all_modules = {**modules_up, **modules_down}

# Filter out modules with too few genes (statistically unstable mean z-score)
MIN_MODULE_GENES = 5
all_modules = {
    k: v for k, v in all_modules.items()
    if len(v['genes']) >= MIN_MODULE_GENES or k.startswith('HZ_Disease')
}
small_modules = {k for k in {**modules_up, **modules_down} if k not in all_modules}
if small_modules:
    print(f"\nDropped modules with < {MIN_MODULE_GENES} genes: {small_modules}")
    # Remove stale txt files from previous runs
    for mod_key in small_modules:
        stale = os.path.join(MOD_DIR, f"module_{mod_key}.txt")
        if os.path.exists(stale):
            os.remove(stale)

# Save module gene lists
for mod_key, mod_info in all_modules.items():
    fname = f"module_{mod_key}.txt"
    with open(os.path.join(MOD_DIR, fname), 'w') as f:
        f.write('\n'.join(mod_info['genes']))

# Save modules as JSON
with open(os.path.join(MOD_DIR, "modules.json"), 'w') as f:
    json.dump({k: {'name': v['name'], 'n_genes': len(v['genes']),
                    'genes': v['genes'], 'terms': v['terms']}
               for k, v in all_modules.items()}, f, indent=2)

print("\nDefined modules:")
for mod_key, mod_info in all_modules.items():
    print(f"  {mod_info['name']:40s} : {len(mod_info['genes']):4d} genes  [{mod_key}]")

# ── 4. Direction A: Score HZ modules in RZV cells ────────────────────────────
print("\n" + "=" * 70)
print("STEP 3: Direction A - HZ modules scored in RZV single cells")
print("=" * 70)

# Load cell counts
print("Loading cell count matrix...")
counts_file = os.path.join(DATA_DIR, "GSE249632_suppl/GSE249632_P261_gene_counts.txt.gz")
cell_counts = pd.read_csv(counts_file, sep='\t', index_col=0).T
print(f"  Matrix: {cell_counts.shape[1]} genes x {cell_counts.shape[0]} cells")

# Build gene symbol → column map for the count matrix
# Use cached Ensembl→symbol mapping from HZ annotation (avoids mygene API dependency)
gene_ids = cell_counts.columns.tolist()
gene_clean = [g.split('.')[0] for g in gene_ids]

# Build mapping from cached annotation (already loaded as ens2sym in Step 1)
sym_map = ens2sym  # from Step 1 HZ annotation loading

# Fallback: also try old HZ annotated file for additional mappings
annot_file = os.path.join(RES_DIR, "GSE242252", "DE_HZ_annotated.csv")
if os.path.exists(annot_file):
    annot = pd.read_csv(annot_file)
    for _, row in annot.iterrows():
        eid = str(row.get('ensembl_id_clean', ''))
        sym = str(row.get('symbol', ''))
        if eid and sym and sym != 'nan' and eid not in sym_map:
            sym_map[eid] = sym

# Build symbol → column names map
sym_to_cols = {}
for i, gid in enumerate(gene_ids):
    clean = gid.split('.')[0]
    sym = sym_map.get(clean, '')
    if sym:
        if sym not in sym_to_cols:
            sym_to_cols[sym] = []
        sym_to_cols[sym].append(gene_ids[i])

print(f"  Mapped {len(sym_to_cols)} unique gene symbols (from cached annotation)")

# Cell metadata
cell_meta = pd.read_csv(os.path.join(RES_DIR, "GSE249632/cell_metadata_clean.csv"))
cell_meta_valid = cell_meta[cell_meta['lib_id'].isin(cell_counts.index)].copy()

# Compute module scores (mean z-score method = AddModuleScore equivalent)
def compute_module_scores(counts, modules, sym_to_cols):
    """Compute module score for each cell for each module."""
    lib_sizes = counts.sum(axis=1)
    cpm = counts.div(lib_sizes, axis=0) * 1e6
    log_cpm = np.log2(cpm + 1)

    scores = pd.DataFrame(index=counts.index)
    for mod_key, mod_info in modules.items():
        cols = []
        for g in mod_info['genes']:
            if g in sym_to_cols:
                cols.extend(sym_to_cols[g])
        cols = list(set(cols))  # deduplicate
        if len(cols) < 3:
            continue
        sub = log_cpm[cols]
        z = (sub - sub.mean(axis=0)) / (sub.std(axis=0) + 1e-8)
        scores[mod_key] = z.mean(axis=1)
    return scores

# Use all cells for scoring
common_libs = list(set(cell_counts.index) & set(cell_meta_valid['lib_id']))
counts_sub = cell_counts.loc[common_libs]
meta_sub = cell_meta_valid[cell_meta_valid['lib_id'].isin(common_libs)].set_index('lib_id')

print(f"Scoring {len(all_modules)} modules in {len(common_libs)} cells...")
cell_scores = compute_module_scores(counts_sub, all_modules, sym_to_cols)
cell_scores = cell_scores.join(meta_sub[['donor', 'timepoint', 'visit']])

print(f"  Modules scored: {len(cell_scores.columns) - 3}")

# ── Plot Direction A ─────────────────────────────────────────────────────────
timepoint_order = ['D0', 'D14', 'D60', 'D74', 'D365']
tp_colors = {'D0': '#1B9E77', 'D14': '#D95F02', 'D60': '#7570B3',
             'D74': '#E7298A', 'D365': '#66A61E'}

print("\nGenerating Direction A plots...")
for mod_key in cell_scores.columns:
    if mod_key in ['donor', 'timepoint', 'visit']:
        continue

    mod_name = all_modules.get(mod_key, {}).get('name', mod_key)

    fig, axes = plt.subplots(1, 2, figsize=(18, 6))

    # Panel 1: Violin
    data_by_tp = [cell_scores[cell_scores['timepoint'] == tp][mod_key].dropna().values
                  for tp in timepoint_order]

    ax = axes[0]
    parts = ax.violinplot(data_by_tp, positions=range(len(timepoint_order)),
                          showmeans=True, showmedians=True, widths=0.7)
    for i, body in enumerate(parts['bodies']):
        body.set_facecolor(tp_colors.get(timepoint_order[i], 'grey'))
        body.set_alpha(0.6)
    ax.set_xticks(range(len(timepoint_order)))
    ax.set_xticklabels(timepoint_order, fontsize=11)
    ax.set_ylabel('Module Score', fontsize=12)
    ax.set_title(f'{mod_name}', fontsize=12, fontweight='bold')

    try:
        h_stat, p_val = kruskal(*data_by_tp)
        ax.text(0.02, 0.98, f'Kruskal-Wallis p = {p_val:.2e}',
                transform=ax.transAxes, va='top', fontsize=9,
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7))
    except:
        pass

    # Panel 2: Per-donor mean trajectory
    ax = axes[1]
    for donor in sorted(cell_scores['donor'].dropna().unique()):
        dd = cell_scores[cell_scores['donor'] == donor]
        dd_tp = dd.groupby('timepoint')[mod_key].mean()
        pts = [(timepoint_order.index(t), dd_tp[t]) for t in dd_tp.index if t in timepoint_order]
        if len(pts) >= 2:
            xs, ys = zip(*sorted(pts))
            ax.plot(xs, ys, 'o-', alpha=0.5, markersize=6, linewidth=1.5,
                    label=f'Donor {int(donor)}')

    ax.set_xticks(range(len(timepoint_order)))
    ax.set_xticklabels(timepoint_order, fontsize=11)
    ax.set_ylabel('Mean Module Score', fontsize=12)
    ax.set_title(f'{mod_name} - Donor Trajectories', fontsize=12)
    if len(cell_scores['donor'].unique()) <= 10:
        ax.legend(fontsize=7, loc='best', ncol=2)

    fig.tight_layout()
    safe_name = mod_key.replace('/', '_').replace(' ', '_')
    fig.savefig(os.path.join(OUT_DIR, f"A_module_{safe_name}.png"), dpi=150, bbox_inches='tight')
    plt.close(fig)

# ── Summary table for Direction A ────────────────────────────────────────────
print("\nDirection A summary statistics:")
summary_a = []
for mod_key in cell_scores.columns:
    if mod_key in ['donor', 'timepoint', 'visit']:
        continue
    row = {'Module': mod_key}
    for tp in timepoint_order:
        vals = cell_scores[cell_scores['timepoint'] == tp][mod_key].dropna()
        row[f'{tp}'] = f"{vals.mean():.4f} ± {vals.std():.4f}"

    # D0 vs D14, D0 vs D74, D0 vs D365
    d0 = cell_scores[cell_scores['timepoint'] == 'D0'][mod_key].dropna()
    for tp in ['D14', 'D60', 'D74', 'D365']:
        tp_vals = cell_scores[cell_scores['timepoint'] == tp][mod_key].dropna()
        if len(d0) > 0 and len(tp_vals) > 0:
            u, p = mannwhitneyu(d0, tp_vals, alternative='two-sided')
            row[f'p_{tp}_vs_D0'] = f"{p:.2e}"
    summary_a.append(row)

summary_a_df = pd.DataFrame(summary_a)
summary_a_df.to_csv(os.path.join(OUT_DIR, "A_module_scores_summary.csv"), index=False)
print(summary_a_df[['Module'] + [c for c in summary_a_df.columns if 'D' in c[:2]]].to_string())

# ── 5. Direction B: RZV modules scored in HZ bulk ────────────────────────────
print("\n" + "=" * 70)
print("STEP 4: Direction B - RZV modules scored in HZ bulk")
print("=" * 70)

# Load RZV DE results to define vaccine response modules
rzv_de_all = {}
for tp in ['D14', 'D60', 'D74', 'D365']:
    tp_file = os.path.join(RES_DIR, f"GSE249632/DE_{tp}_vs_D0.csv")
    if os.path.exists(tp_file):
        rzv_de_all[tp] = pd.read_csv(tp_file)

# Map RZV genes to symbols
rzv_gene_ids = set()
for df in rzv_de_all.values():
    rzv_gene_ids.update(df['gene_id'].str.split('.').str[0].tolist())

print(f"Mapping {len(rzv_gene_ids)} RZV gene IDs to symbols (cached)...")
# Use the symbol mapping already built (ens2sym from Step 1 + annotation file)
rzv_sym_map = sym_map  # reuse the comprehensive cached mapping

# Define RZV modules
rzv_modules = {}

# RZV Activation (D14 up-regulated)
if 'D14' in rzv_de_all:
    d14 = rzv_de_all['D14']
    d14_up_ids = d14[(d14['adj.P.Val'] < 0.05) & (d14['logFC'] > 0.5)]['gene_id'].str.split('.').str[0]
    genes = [rzv_sym_map.get(g, '') for g in d14_up_ids]
    genes = [g for g in genes if g]
    if genes:
        rzv_modules['RZV_Acute_Activation'] = {
            'name': 'RZV Acute Activation (D14 up)',
            'genes': sorted(set(genes)),
            'terms': [f'{len(genes)} genes up-regulated at D14 vs D0']
        }

# RZV Persistent Regulation (genes up at D14 AND D74 AND D365)
persistent_ids = None
for tp in ['D14', 'D74', 'D365']:
    if tp in rzv_de_all:
        df = rzv_de_all[tp]
        tp_up = set(df[(df['adj.P.Val'] < 0.05) & (df['logFC'] > 0)]['gene_id'].str.split('.').str[0])
        if persistent_ids is None:
            persistent_ids = tp_up
        else:
            persistent_ids &= tp_up

if persistent_ids:
    genes = [rzv_sym_map.get(g, '') for g in persistent_ids]
    genes = [g for g in genes if g]
    rzv_modules['RZV_Persistent'] = {
        'name': 'RZV Persistent Regulation (D14+D74+D365 up)',
        'genes': sorted(set(genes)),
        'terms': [f'{len(genes)} genes up at all timepoints']
    }

# RZV Long-term Imprint (D365 specific)
if 'D365' in rzv_de_all:
    d365 = rzv_de_all['D365']
    d365_sig = d365[(d365['adj.P.Val'] < 0.05) & (abs(d365['logFC']) > 0.5)]
    d365_up = set(d365_sig[d365_sig['logFC'] > 0]['gene_id'].str.split('.').str[0])
    d365_down = set(d365_sig[d365_sig['logFC'] < 0]['gene_id'].str.split('.').str[0])

    up_genes = [rzv_sym_map.get(g, '') for g in d365_up]
    up_genes = [g for g in up_genes if g]
    down_genes = [rzv_sym_map.get(g, '') for g in d365_down]
    down_genes = [g for g in down_genes if g]

    if up_genes:
        rzv_modules['RZV_D365_Up'] = {
            'name': 'RZV Long-term Up (D365)',
            'genes': sorted(set(up_genes)),
            'terms': [f'{len(up_genes)} genes up at D365']
        }
    if down_genes:
        rzv_modules['RZV_D365_Down'] = {
            'name': 'RZV Long-term Down (D365)',
            'genes': sorted(set(down_genes)),
            'terms': [f'{len(down_genes)} genes down at D365']
        }

print("RZV modules:")
for k, v in rzv_modules.items():
    print(f"  {v['name']:45s} : {len(v['genes']):4d} genes")

# Load HZ bulk data
hz_counts = pd.read_csv(os.path.join(RES_DIR, "GSE242252/filtered_counts.csv"), index_col=0)
hz_meta_raw = pd.read_csv(os.path.join(RES_DIR, "GSE242252/deseq2_metadata.csv"), index_col=0)
hz_samples = hz_meta_raw[hz_meta_raw['condition_label'] == 'Herpes_Zoster'].index
hz_counts_hz = hz_counts[hz_samples]
hz_meta_hz = hz_meta_raw.loc[hz_samples]

# Build gene symbol map for HZ counts
print("Building HZ gene symbol map...")
hz_gene_sym = {}
for i, gid in enumerate(hz_counts_hz.index):
    clean = gid.split('.')[0]
    # We already have this mapping from the earlier symbol mapping step
    hz_gene_sym[gid] = hz_degs.loc[hz_degs['gene_id'] == gid, 'symbol'].values
    hz_gene_sym[gid] = hz_gene_sym[gid][0] if len(hz_gene_sym[gid]) > 0 else ''

# Alternative: use the annotated file
hz_annot = pd.read_csv(os.path.join(RES_DIR, "GSE242252/DE_HZ_annotated.csv"))
hz_sym_lookup = dict(zip(hz_annot['gene_id'], hz_annot['symbol']))

# Score RZV modules in HZ bulk
print("\nScoring RZV modules in HZ bulk...")
# CPM normalize
hz_cpm = hz_counts_hz.div(hz_counts_hz.sum(axis=0), axis=1) * 1e6
hz_log_cpm = np.log2(hz_cpm + 1)

rzv_scores = pd.DataFrame(index=hz_samples)
for mod_key, mod_info in rzv_modules.items():
    mod_symbols = mod_info['genes']
    cols = []
    for sym in mod_symbols:
        matching = hz_annot[hz_annot['symbol'] == sym]['gene_id']
        for gid in matching:
            if gid in hz_log_cpm.index:
                cols.append(gid)
    cols = list(set(cols))

    if len(cols) < 3:
        print(f"  {mod_key}: only {len(cols)} genes mapped, skipping")
        continue

    sub = hz_log_cpm.loc[cols]
    z = (sub.T - sub.mean(axis=1)) / (sub.std(axis=1) + 1e-8)
    rzv_scores[mod_key] = z.mean(axis=1)
    print(f"  {mod_key}: {len(cols)} genes scored")

rzv_scores['timepoint'] = hz_meta_hz['timepoint'].values
rzv_scores['patient_id'] = hz_meta_hz['patient_id'].values

# ── Plot Direction B ─────────────────────────────────────────────────────────
print("\nGenerating Direction B plots...")
for mod_key in rzv_scores.columns:
    if mod_key in ['timepoint', 'patient_id']:
        continue

    mod_name = rzv_modules.get(mod_key, {}).get('name', mod_key)

    fig, ax = plt.subplots(figsize=(7, 6))

    acute = rzv_scores[rzv_scores['timepoint'] == 'acute'][mod_key].dropna()
    conv = rzv_scores[rzv_scores['timepoint'] == 'convalescent'][mod_key].dropna()

    bp = ax.boxplot([acute, conv], labels=['Acute HZ', 'Convalescent (1 yr)'],
                     patch_artist=True, widths=0.5, showfliers=True)
    bp['boxes'][0].set_facecolor('#E41A1C')
    bp['boxes'][1].set_facecolor('#377EB8')
    for box in bp['boxes']:
        box.set_alpha(0.7)

    # Add paired lines
    paired = rzv_scores.pivot(index='patient_id', columns='timepoint', values=mod_key)
    paired = paired.dropna()
    for _, row in paired.iterrows():
        ax.plot([0, 1], [row['acute'], row['convalescent']],
                'o-', color='grey', alpha=0.3, markersize=4)

    if len(paired) > 3:
        _, p = wilcoxon(paired['acute'], paired['convalescent'])
        sig_str = f'p = {p:.4f}' + (' *' if p < 0.05 else '')
        ax.text(0.5, 0.96, sig_str, transform=ax.transAxes, ha='center', va='top',
                fontsize=11, fontweight='bold',
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7))

    ax.set_ylabel(f'{mod_name} Score', fontsize=11)
    ax.set_title(f'RZV Module: {mod_name}\nScored in HZ Bulk (Acute vs Convalescent)',
                 fontsize=10, fontweight='bold')
    fig.tight_layout()
    safe_name = mod_key.replace('/', '_').replace(' ', '_')
    fig.savefig(os.path.join(OUT_DIR, f"B_rzv_module_{safe_name}.png"), dpi=150, bbox_inches='tight')
    plt.close(fig)

# Save
rzv_scores.to_csv(os.path.join(OUT_DIR, "B_rzv_modules_in_HZ.csv"), index=False)

print("\nDirection B summary:")
for mod_key in rzv_scores.columns:
    if mod_key in ['timepoint', 'patient_id']:
        continue
    acute = rzv_scores[rzv_scores['timepoint'] == 'acute'][mod_key]
    conv = rzv_scores[rzv_scores['timepoint'] == 'convalescent'][mod_key]
    _, p = mannwhitneyu(acute, conv, alternative='two-sided') if len(acute) > 0 and len(conv) > 0 else (None, 1)
    print(f"  {mod_key:30s}  Acute: {acute.mean():+.4f}  Conv: {conv.mean():+.4f}  p={p:.4f}")

print("\n" + "=" * 70)
print("Module scoring complete!")
print(f"Output in: {OUT_DIR}")
