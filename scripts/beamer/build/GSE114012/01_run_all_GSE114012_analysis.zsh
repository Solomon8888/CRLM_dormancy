#!/usr/bin/env zsh
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
cd "$PROJECT_ROOT"

RESULT_ROOT="results/ngs/GSE114012"
GSEA_QS2_CACHE_ROOT="temporary/ngs/GSE114012/GSEA_qs2_cache"

rm -rf \
  "${RESULT_ROOT}/tables" \
  "${RESULT_ROOT}/plots" \
  "${RESULT_ROOT}/intersect" \
  "${RESULT_ROOT}/TF" \
  "${RESULT_ROOT}/TF_summary" \
  "${GSEA_QS2_CACHE_ROOT}"
mkdir -p "$RESULT_ROOT"

Rscript scripts/GSE114012/00_sample_clustering_heatmap.R
Rscript scripts/GSE114012/01_limma_differential_expression.R
Rscript scripts/GSE114012/02_intersect_significant_genes.R
Rscript scripts/GSE114012/03_volcano_plot.R
Rscript scripts/GSE114012/04_multiple_volcano_plot.R
Rscript scripts/GSE114012/05_top_deg_gene_heatmap.R
Rscript scripts/GSE114012/06_gsea_analysis.R
Rscript scripts/GSE114012/07_gsea_plot.R
Rscript scripts/GSE114012/08_tf_enrichment_analysis.R
Rscript scripts/GSE114012/09_integrate_tf_enrichment_results.R
