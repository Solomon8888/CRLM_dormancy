#!/usr/bin/env zsh
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
cd "$PROJECT_ROOT"

source scripts/functions/zsh_runtime_functions.zsh

# GSE114012数据集Beamer TeX构建入口。
# 仅刷新scripts/beamer/beamer_report.tex与sections/GSE114012，不执行LaTeX编译。
LOG_FILE="temporary/beamer/logs/GSE114012_build_tex.log"
print_step_header "Build GSE114012 Beamer TeX"
run_with_spinner \
  "Generating Beamer sources for GSE114012" \
  "$LOG_FILE" \
  python3 scripts/functions/beamer_dataset_report_builder.py
