#!/usr/bin/env zsh
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
cd "$PROJECT_ROOT"

LOG_ROOT="temporary/GSE114012/full_analysis_logs/$(date +%Y%m%d_%H%M%S)"
mkdir -p "$LOG_ROOT"
RESULT_ROOT="results/ngs/GSE114012"
GSEA_QS2_CACHE_ROOT="temporary/ngs/GSE114012/GSEA_qs2_cache"

print_step() {
  printf '\n[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$1"
}

run_r_script() {
  local script_path="$1"
  local script_name
  script_name="$(basename "$script_path" .R)"

  print_step "Running ${script_path}"
  /usr/bin/time -p Rscript "$script_path" > "${LOG_ROOT}/${script_name}.log" 2>&1
  print_step "Finished ${script_path}"
}

run_parallel_group() {
  local group_name="$1"
  shift

  print_step "Starting group: ${group_name}"
  local -a pids=()
  local script_path
  for script_path in "$@"; do
    (
      run_r_script "$script_path"
    ) &
    pids+=("$!")
  done

  local pid
  for pid in "${pids[@]}"; do
    wait "$pid"
  done
  print_step "Finished group: ${group_name}"
}

START_TS="$(date +%s)"

print_step "GSE114012 full analysis started"
print_step "Logs: ${LOG_ROOT}"

print_step "Cleaning previous GSE114012 analysis outputs"
rm -rf \
  "${RESULT_ROOT}/tables" \
  "${RESULT_ROOT}/plots" \
  "${RESULT_ROOT}/intersect" \
  "${RESULT_ROOT}/TF" \
  "${RESULT_ROOT}/TF_summary" \
  "${GSEA_QS2_CACHE_ROOT}"
mkdir -p "$RESULT_ROOT"

# 00与01互不依赖，可并行启动；01是后续分析的核心依赖。
run_parallel_group "quality_control_and_deg" \
  "scripts/GSE114012/00_sample_clustering_heatmap.R" \
  "scripts/GSE114012/01_limma_differential_expression.R"

# 02-06均依赖01；其中06内部会自行并行运行GSEA。
run_parallel_group "post_deg_tables_and_plots" \
  "scripts/GSE114012/02_intersect_significant_genes.R" \
  "scripts/GSE114012/03_volcano_plot.R" \
  "scripts/GSE114012/04_multiple_volcano_plot.R" \
  "scripts/GSE114012/05_top_deg_gene_heatmap.R" \
  "scripts/GSE114012/06_gsea_analysis.R"

# 07依赖06；08依赖01/02。二者均较重，放在同一批次里让各自脚本内部调度核心。
run_parallel_group "gsea_plots_and_tf_enrichment" \
  "scripts/GSE114012/07_gsea_plot.R" \
  "scripts/GSE114012/08_tf_enrichment_analysis.R"

# 09依赖08的全部TF输出。
run_r_script "scripts/GSE114012/09_integrate_tf_enrichment_results.R"

END_TS="$(date +%s)"
RUNTIME=$((END_TS - START_TS))

print_step "GSE114012 full analysis finished in ${RUNTIME}s"
print_step "Logs: ${LOG_ROOT}"
