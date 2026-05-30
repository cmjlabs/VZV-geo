#!/bin/bash
# VZV-geo full pipeline — run from project root
# Usage: bash scripts/run_all.sh
set -e
cd "$(dirname "$0")/.."
echo "============================================"
echo "VZV-geo Full Pipeline — $(date)"
echo "============================================"

# Step 1: Build GSE242252 count matrix
echo ""
echo "=== Step 1/8: Build GSE242252 count matrix ==="
python3 scripts/02_build_GSE242252_matrix.py

# Step 2: GSE242252 DESeq2
echo ""
echo "=== Step 2/8: GSE242252 DESeq2 ==="
Rscript scripts/03_GSE242252_DESeq2.R

# Step 3: Gene symbol annotation
echo ""
echo "=== Step 3/8: Gene symbol annotation ==="
python3 scripts/09_annotate_genes.py

# Step 4: GSE249632 pseudobulk (QC-only)
echo ""
echo "=== Step 4/8: GSE249632 pseudobulk ==="
python3 scripts/04_GSE249632_process.py

# Step 5: GSE249632 DEA
echo ""
echo "=== Step 5/8: GSE249632 DEA ==="
Rscript scripts/05_GSE249632_pseudobulk_DEA.R

# Step 6: Pathway enrichment + module scoring + cross-dataset comparison
echo ""
echo "=== Step 6/8: Module scoring + enrichment ==="
unset ALL_PROXY all_proxy HTTP_PROXY http_proxy HTTPS_PROXY https_proxy
python3 scripts/07_module_scoring_enrichment.py

# Step 7: Comparison figures (quadrant, IFN, vaccine immunology)
echo ""
echo "=== Step 7/8: Comparison figures ==="
unset ALL_PROXY all_proxy HTTP_PROXY http_proxy HTTPS_PROXY https_proxy
python3 scripts/10_comparison_figures.py

# Step 8: HTML report
echo ""
echo "=== Step 8/8: HTML report ==="
python3 scripts/08_generate_html_report.py

echo ""
echo "============================================"
echo "Pipeline complete — $(date)"
echo "Output: results/comprehensive_report.html"
echo "============================================"
