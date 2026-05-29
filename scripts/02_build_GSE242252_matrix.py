#!/usr/bin/env python3
"""
Build count matrix from 80 individual GSE242252 sample files.
Parse sample metadata: CO=Control, HZ=Herpes Zoster, ERC=External Ref Control
-1Y suffix = 1 year later (convalescent for HZ, follow-up for controls)
"""
import os, re, gzip
import pandas as pd
import numpy as np

DATA_DIR = "/media/cmj/MechanicalDisk/yjs/VZV-geo/data"
SUPP_DIR = os.path.join(DATA_DIR, "GSE242252_suppl")
OUT_DIR = "/media/cmj/MechanicalDisk/yjs/VZV-geo/results/GSE242252"
os.makedirs(OUT_DIR, exist_ok=True)

# ── 1. Scan all sample files ─────────────────────────────────────────────────
sample_dirs = sorted(os.listdir(SUPP_DIR))
sample_files = []
for d in sample_dirs:
    full = os.path.join(SUPP_DIR, d)
    if os.path.isdir(full):
        files = os.listdir(full)
        if files:
            sample_files.append((d, os.path.join(full, files[0])))

print(f"Found {len(sample_files)} sample files.")

# ── 2. Parse metadata from GSM ID / filename ─────────────────────────────────
def parse_metadata(dir_name):
    """Parse condition, patient ID, and timepoint from directory name."""
    # Format: Supp_GSMnnnnnnn_CONDITION_PATIENT[_1Y]_Sxx
    name = dir_name.replace("Supp_", "")
    parts = name.split("_")
    gsm_id = parts[0]
    condition = parts[1]  # CO, HZ, or ERC

    # Check for -1Y suffix in patient part
    patient_str = parts[2]
    timepoint = "acute"  # default
    if patient_str.endswith("-1Y"):
        patient_id = patient_str.replace("-1Y", "")
        timepoint = "convalescent"
    else:
        patient_id = patient_str

    return {
        "gsm": gsm_id,
        "dir_name": dir_name,
        "condition": condition,
        "condition_label": {"CO": "Control", "HZ": "Herpes_Zoster",
                           "ERC": "External_Control"}.get(condition, condition),
        "patient_id": patient_id,
        "timepoint": timepoint,
        "is_paired": condition in ("CO", "HZ")
    }

meta_list = [parse_metadata(d) for d, f in sample_files]
meta_df = pd.DataFrame(meta_list)

print("\nSample distribution:")
print(meta_df.groupby(["condition_label", "timepoint"]).size().to_string())

# ── 3. Read and merge all gene count files ───────────────────────────────────
print("\nReading gene count files...")
all_genes = None
count_dict = {}

for (dir_name, filepath), meta in zip(sample_files, meta_list):
    counts = {}
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                gene_id = parts[0]
                # Remove version number from Ensembl ID (e.g., ENSG00000173611.18 -> ENSG00000173611)
                gene_id_clean = gene_id.split('.')[0]
                counts[gene_id_clean] = int(parts[1])

    count_dict[meta["gsm"]] = counts
    if all_genes is None:
        all_genes = sorted(counts.keys())

# Build matrix
print(f"\nBuilding count matrix: {len(all_genes)} genes x {len(count_dict)} samples")
mat = np.zeros((len(all_genes), len(count_dict)), dtype=np.int32)
gene_index = {g: i for i, g in enumerate(all_genes)}
gsm_ids = [m["gsm"] for m in meta_list]

for j, gsm in enumerate(gsm_ids):
    counts = count_dict[gsm]
    for gene, idx in gene_index.items():
        mat[idx, j] = counts.get(gene, 0)

count_matrix = pd.DataFrame(mat, index=all_genes, columns=gsm_ids)
count_matrix.index.name = "gene_id"

# ── 4. Filter low-count genes (at least 10 counts in >= 3 samples) ──────────
keep = (count_matrix >= 10).sum(axis=1) >= 3
print(f"Genes passing min-count filter: {keep.sum()} / {len(keep)}")
count_matrix_filt = count_matrix.loc[keep]

# ── 5. Save ──────────────────────────────────────────────────────────────────
count_matrix.to_csv(os.path.join(OUT_DIR, "raw_counts.csv"))
count_matrix_filt.to_csv(os.path.join(OUT_DIR, "filtered_counts.csv"))
meta_df.to_csv(os.path.join(OUT_DIR, "sample_metadata.csv"), index=False)

# Build a DESeq2-ready metadata
deseq_meta = meta_df.copy()
deseq_meta["group"] = deseq_meta["condition_label"] + "_" + deseq_meta["timepoint"]
deseq_meta = deseq_meta.set_index("gsm")
deseq_meta.to_csv(os.path.join(OUT_DIR, "deseq2_metadata.csv"))

print("\nCount matrix shapes:")
print(f"  Raw: {count_matrix.shape}")
print(f"  Filtered: {count_matrix_filt.shape}")
print(f"\nFiles saved in {OUT_DIR}")
