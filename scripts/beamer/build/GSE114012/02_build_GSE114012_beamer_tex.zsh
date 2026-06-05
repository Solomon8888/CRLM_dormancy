#!/usr/bin/env zsh
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
cd "$PROJECT_ROOT"

# GSE114012数据集Beamer TeX构建入口。
# 仅刷新scripts/beamer/beamer_report.tex与sections/GSE114012，不执行LaTeX编译。
python3 scripts/functions/beamer_dataset_report_builder.py
