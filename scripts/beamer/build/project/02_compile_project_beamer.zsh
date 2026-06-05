#!/usr/bin/env zsh
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
cd "$PROJECT_ROOT"

source scripts/functions/zsh_runtime_functions.zsh

# 项目级Beamer编译入口。
# 中间文件保存在temporary/beamer；最终PDF保存在results/reports/beamer。
LOG_FILE="temporary/beamer/logs/project_compile.log"
print_step_header "Compile Project Beamer PDF"
run_with_spinner \
  "Compiling existing Beamer TeX to PDF" \
  "$LOG_FILE" \
  python3 scripts/functions/beamer_dataset_report_builder.py --compile-only
