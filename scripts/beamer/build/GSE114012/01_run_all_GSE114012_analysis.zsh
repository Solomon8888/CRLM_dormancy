#!/usr/bin/env zsh
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
cd "$PROJECT_ROOT"

RESULT_ROOT="results/ngs/GSE114012"
GSEA_QS2_CACHE_ROOT="temporary/ngs/GSE114012/GSEA_qs2_cache"
RSCRIPT_BIN="${RSCRIPT_BIN:-$(command -v Rscript)}"

if [[ -z "$RSCRIPT_BIN" || ! -x "$RSCRIPT_BIN" ]]; then
  printf "Cannot find executable Rscript. Set RSCRIPT_BIN=/path/to/Rscript and rerun.\n" >&2
  exit 1
fi

R_SCRIPTS=(
  "scripts/GSE114012/00_sample_clustering_heatmap.R"
  "scripts/GSE114012/01_limma_differential_expression.R"
  "scripts/GSE114012/02_intersect_significant_genes.R"
  "scripts/GSE114012/03_volcano_plot.R"
  "scripts/GSE114012/04_multiple_volcano_plot.R"
  "scripts/GSE114012/05_top_deg_gene_heatmap.R"
  "scripts/GSE114012/06_gsea_analysis.R"
  "scripts/GSE114012/07_gsea_plot.R"
  "scripts/GSE114012/08_tf_enrichment_analysis.R"
  "scripts/GSE114012/09_integrate_tf_enrichment_results.R"
)

run_r_script() {
  local script="$1"
  local index="$2"
  local total="$3"

  if [[ ! -f "$script" ]]; then
    printf "Missing script: %s\n" "$script" >&2
    return 1
  fi

  printf "\n============================================================\n"
  printf "[%02d/%02d] Running: %s\n" "$index" "$total" "$script"
  printf "============================================================\n"
  "$RSCRIPT_BIN" "$script"
}

rm -rf \
  "${RESULT_ROOT}/tables" \
  "${RESULT_ROOT}/plots" \
  "${RESULT_ROOT}/intersect" \
  "${RESULT_ROOT}/TF" \
  "${RESULT_ROOT}/TF_summary" \
  "${GSEA_QS2_CACHE_ROOT}" \
  omnipathr-log \
  omipathr-log \
  manuscripts
mkdir -p "$RESULT_ROOT"

printf "\nRunning GSE114012 full analysis pipeline...\n"
printf "Project root: %s\n" "$PROJECT_ROOT"
printf "Rscript binary: %s\n" "$RSCRIPT_BIN"

total_scripts="${#R_SCRIPTS[@]}"

for (( i = 1; i <= total_scripts; i++ )); do
  run_r_script "${R_SCRIPTS[$i]}" "$i" "$total_scripts"
done

printf "\nGSE114012 full analysis pipeline finished.\n"
