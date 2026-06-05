#!/usr/bin/env zsh
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
cd "$PROJECT_ROOT"

# 项目级Beamer编译入口。
# 中间文件保存在temporary/beamer；最终PDF保存在results/reports/beamer。
python3 scripts/functions/beamer_dataset_report_builder.py --compile
