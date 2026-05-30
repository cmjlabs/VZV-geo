#!/usr/bin/env python3
"""Annotate DESeq2 results with gene symbols using cached Ensembl→symbol mapping."""
import os, sys
import pandas as pd

# Project root
PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(PROJ_ROOT, "results", "GSE242252")

# Load new DESeq2 output
new = pd.read_csv(os.path.join(RES, "DE_HZ_acute_vs_convalescent.csv"))
ensembl_ids = new['gene_id'].str.split('.').str[0].tolist()

# Build symbol map from old annotation or use mygene fallback
old_file = os.path.join(RES, "DE_HZ_annotated.csv")
ens2sym = {}
if os.path.exists(old_file):
    old = pd.read_csv(old_file)
    for _, row in old.iterrows():
        eid = str(row.get('ensembl_id_clean', ''))
        sym = str(row.get('symbol', ''))
        if eid and sym and sym != 'nan':
            ens2sym[eid] = sym
    print(f"Loaded {len(ens2sym)} symbol mappings from existing annotation")

# Try mygene for unmapped genes
unmapped = [e for e in set(ensembl_ids) if e not in ens2sym]
if unmapped:
    try:
        import mygene
        mg = mygene.MyGeneInfo()
        results = mg.querymany(unmapped, scopes='ensembl.gene', fields='symbol',
                               species='human', batch_size=1000)
        for r in results:
            sym = r.get('symbol', '')
            if sym:
                ens2sym[r['query']] = sym
        print(f"  MyGene mapped {sum(1 for r in results if r.get('symbol'))} additional genes")
    except Exception as e:
        print(f"  MyGene unavailable: {e}")

# Annotate
new['ensembl_id_clean'] = ensembl_ids
new['symbol'] = [ens2sym.get(eid, '') for eid in ensembl_ids]

n_sig = sum((new['padj'].notna()) & (new['padj'] < 0.05))
n_up = sum((new['padj'].notna()) & (new['padj'] < 0.05) & (new['log2FoldChange'] > 0))
n_down = sum((new['padj'].notna()) & (new['padj'] < 0.05) & (new['log2FoldChange'] < 0))
print(f"Genes: {len(new)} total, {n_sig} sig (FDR<0.05): {n_up} up, {n_down} down")
print(f"With symbols: {new['symbol'].notna().sum()} / {len(new)}")

new.to_csv(os.path.join(RES, "DE_HZ_annotated.csv"), index=False)

# Comparison file
lfc = new[['gene_id', 'ensembl_id_clean', 'symbol', 'log2FoldChange', 'padj']].copy()
lfc.columns = ['gene_id', 'ensembl_id', 'symbol', 'log2FC_HZ', 'padj_HZ']
lfc.to_csv(os.path.join(RES, "HZ_log2FC_for_comparison.csv"), index=False)
print("Saved: DE_HZ_annotated.csv, HZ_log2FC_for_comparison.csv")
