#!/usr/bin/env python3
"""Download GSE242252 and GSE249632 from GEO."""
import os
import GEOparse

DATA_DIR = "/media/cmj/MechanicalDisk/yjs/VZV-geo/data"

os.makedirs(DATA_DIR, exist_ok=True)

# --- GSE242252: Bulk RNA-seq, HZ whole blood ---
print("Downloading GSE242252...")
gse_bulk = GEOparse.get_GEO(geo="GSE242252", destdir=DATA_DIR, silent=True)
print(f"  Title: {gse_bulk.metadata['title'][0]}")
print(f"  Samples: {len(gse_bulk.gsms)}")
# Save expression matrix if available
try:
    gse_bulk.download_supplementary_files(directory=os.path.join(DATA_DIR, "GSE242252_suppl"))
    print("  Supplementary files downloaded.")
except Exception as e:
    print(f"  Supplementary download note: {e}")

# --- GSE249632: scRNA-seq, RZV CD4+ T cells ---
print("\nDownloading GSE249632...")
gse_sc = GEOparse.get_GEO(geo="GSE249632", destdir=DATA_DIR, silent=True)
print(f"  Title: {gse_sc.metadata['title'][0]}")
print(f"  Samples: {len(gse_sc.gsms)}")
try:
    gse_sc.download_supplementary_files(directory=os.path.join(DATA_DIR, "GSE249632_suppl"))
    print("  Supplementary files downloaded.")
except Exception as e:
    print(f"  Supplementary download note: {e}")

print("\nDone. Check data/ directory.")
