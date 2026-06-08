#!/usr/bin/env zsh
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
cd "$PROJECT_ROOT"

# quickanalysis脚本自身会按各自逻辑清理前次运行结果与中间文件。
# 这里不额外重定向输出，保留每个R脚本的单行进度条、完成路径和耗时统计。

export QUICKANALYSIS_PARALLEL_BACKEND="${QUICKANALYSIS_PARALLEL_BACKEND:-${PARALLEL_RUNTIME_BACKEND:-auto}}"
if [[ -n "${PARALLEL_RUNTIME_WORKERS:-}" && -z "${QUICKANALYSIS_PARALLEL_WORKERS:-}" ]]; then
  export QUICKANALYSIS_PARALLEL_WORKERS="$PARALLEL_RUNTIME_WORKERS"
fi

R_SCRIPTS=(
  "scripts/quickanalysis/01_tcgaplot_quick_analysis.R"
  "scripts/quickanalysis/02_slcptac_quick_analysis.R"
  "scripts/quickanalysis/03_local_tcga_01a_median_de_gsea.R"
  "scripts/quickanalysis/04_local_tcga_01a_correlation_gsea.R"
  "scripts/quickanalysis/05_local_gtex_median_de_gsea.R"
  "scripts/quickanalysis/06_local_gtex_correlation_gsea.R"
)

for script in "${R_SCRIPTS[@]}"; do
  printf "\n==== Running %s ====\n" "$script"
  Rscript "$script"
done

printf "\nAll quickanalysis scripts finished.\n"
