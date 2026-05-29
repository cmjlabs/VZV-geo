#!/usr/bin/env python3
"""
Process GSE249632 scRNA-seq data:
- Parse metadata (donor, visit, tetramer) from GEO, match to count matrix library IDs
- QC filter (paper criteria: >=500k reads, >=70% alignment, median cov = 1)
- Build pseudobulk profiles per donor x visit
"""
import os, re, gzip
import pandas as pd
import numpy as np
import GEOparse

DATA_DIR = "/media/cmj/MechanicalDisk/yjs/VZV-geo/data"
OUT_DIR = "/media/cmj/MechanicalDisk/yjs/VZV-geo/results/GSE249632"
os.makedirs(OUT_DIR, exist_ok=True)

# ── 1. Read count matrix ─────────────────────────────────────────────────────
print("Reading count matrix...")
counts_file = os.path.join(DATA_DIR, "GSE249632_suppl", "GSE249632_P261_gene_counts.txt.gz")
df = pd.read_csv(counts_file, sep='\t', index_col=0)
print(f"Matrix: {df.shape[0]} genes x {df.shape[1]} cells")

# ── 2. Parse GEO metadata and map to library IDs ─────────────────────────────
print("Parsing GEO metadata...")
gse = GEOparse.get_GEO(geo="GSE249632", destdir=DATA_DIR, silent=True)

cell_meta = []
for gsm_id, gsm in gse.gsms.items():
    title = gsm.metadata.get('title', [''])[0]

    # Extract library ID from title (e.g., "2_DR0701_Baseline_lib40388")
    lib_match = re.search(r'lib(\d+)', title)
    lib_id = f"lib{lib_match.group(1)}" if lib_match else None

    # Parse characteristics
    chars = {}
    for c in gsm.metadata.get('characteristics_ch1', []):
        if ':' in c:
            key, val = c.split(':', 1)
            chars[key.strip().lower()] = val.strip()

    # Get QC info from description
    desc = gsm.metadata.get('description', [''])[0]

    cell_meta.append({
        'gsm': gsm_id,
        'lib_id': lib_id,
        'title': title,
        'donor': chars.get('donorid', ''),
        'visit': chars.get('visit', ''),
        'tetramer': chars.get('tetramer', ''),
        'cell_type': chars.get('cell type', ''),
        'qc_pass': chars.get('qc pass', ''),
    })

meta_df = pd.DataFrame(cell_meta)

# Map visit names to standard format
visit_map = {
    'Baseline': 'D0',
    'Day 14': 'D14',
    'Day 60': 'D60',
    'Day 74': 'D74',
    'Day 365': 'D365'
}
meta_df['timepoint'] = meta_df['visit'].map(visit_map)

print(f"\nMetadata entries: {len(meta_df)}")
print(f"Unique lib_ids: {meta_df['lib_id'].nunique()}")
print(f"lib_ids matching count matrix columns: {meta_df['lib_id'].isin(df.columns).sum()}")

# Check which lib_ids don't match
missing_libs = set(meta_df['lib_id'].dropna()) - set(df.columns)
if missing_libs:
    print(f"lib_ids NOT in count matrix: {len(missing_libs)}")
    if len(missing_libs) < 10:
        print(f"  {missing_libs}")

# ── 3. Summary statistics ───────────────────────────────────────────────────
print("\n=== Sample Distribution ===")
print("By donor:")
print(meta_df['donor'].value_counts().sort_index().to_string())
print("\nBy timepoint:")
print(meta_df['timepoint'].value_counts().to_string())
print("\nBy QC pass:")
print(meta_df['qc_pass'].value_counts().to_string())
print("\nBy tetramer:")
print(meta_df['tetramer'].value_counts().to_string())

# ── 4. Filter to QC-pass cells with valid library ID ─────────────────────────
meta_valid = meta_df[meta_df['lib_id'].notna() & meta_df['lib_id'].isin(df.columns)].copy()
print(f"\nCells with valid lib_id in count matrix: {len(meta_valid)}")

# QC filter
meta_qc = meta_valid[meta_valid['qc_pass'] == 'TRUE'].copy()
print(f"Cells passing QC: {len(meta_qc)}")

# ── 5. Build pseudobulk per donor x timepoint ────────────────────────────────
# Use ONLY QC-pass cells (matches paper: >=500k reads, >=70% alignment, median cov=1)
meta_use = meta_qc.copy()

# Subset count matrix to available cells
available_libs = [l for l in meta_use['lib_id'] if l in df.columns]
counts_sub = df[available_libs]
meta_use = meta_use[meta_use['lib_id'].isin(available_libs)].copy()

# Create pseudobulk groups
meta_use['pb_group'] = meta_use['donor'] + '_' + meta_use['timepoint']
pb_groups = meta_use['pb_group'].unique()
print(f"\nPseudobulk groups (donor_timepoint): {len(pb_groups)}")

# Sum counts per group
pb_matrix = np.zeros((counts_sub.shape[0], len(pb_groups)), dtype=np.int32)
pb_meta_rows = []
for i, grp in enumerate(pb_groups):
    mask = meta_use['pb_group'] == grp
    libs = meta_use.loc[mask, 'lib_id'].tolist()
    pb_matrix[:, i] = counts_sub[libs].sum(axis=1).values.flatten()
    pb_meta_rows.append({
        'group': grp,
        'donor': meta_use.loc[mask, 'donor'].iloc[0],
        'timepoint': meta_use.loc[mask, 'timepoint'].iloc[0],
        'n_cells': mask.sum()
    })

pb_df = pd.DataFrame(pb_matrix, index=counts_sub.index, columns=pb_groups)
pb_meta_df = pd.DataFrame(pb_meta_rows)

print(f"Pseudobulk matrix: {pb_df.shape[0]} genes x {pb_df.shape[1]} samples")
print(f"Cells per pseudobulk sample: median={pb_meta_df['n_cells'].median():.0f}, range={pb_meta_df['n_cells'].min()}-{pb_meta_df['n_cells'].max()}")

print("\nPseudobulk groups per timepoint:")
print(pb_meta_df.groupby('timepoint').size().to_string())

# ── 6. Save ──────────────────────────────────────────────────────────────────
counts_sub.to_csv(os.path.join(OUT_DIR, "cell_counts_qcpass.csv"))
pb_df.to_csv(os.path.join(OUT_DIR, "pseudobulk_counts.csv"))
pb_meta_df.to_csv(os.path.join(OUT_DIR, "pseudobulk_metadata.csv"), index=False)
meta_use.to_csv(os.path.join(OUT_DIR, "cell_metadata_clean.csv"), index=False)

print(f"\nAll files saved to {OUT_DIR}")
