#!/usr/bin/env zsh
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
cd "$PROJECT_ROOT"

source scripts/functions/zsh_runtime_functions.zsh

# 项目级Beamer构建入口。
# 后续新增数据集时，在这里继续串联各数据集的02号构建脚本即可。
LOG_FILE="temporary/beamer/logs/project_build_tex.log"
print_step_header "Build Project Beamer TeX"
run_with_spinner \
  "Building project Beamer source tree" \
  "$LOG_FILE" \
  zsh scripts/beamer/build/GSE114012/02_build_GSE114012_beamer_tex.zsh
