#!/usr/bin/env zsh
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
cd "$PROJECT_ROOT"

source scripts/functions/zsh_runtime_functions.zsh

# 仅清理可再生成内容：
#   1. results 下全部分析结果、图片和最终报告
#   2. temporary 下全部中间文件、日志、LaTeX构建文件和临时表格
#   3. Beamer构建脚本自动生成的tex源文件
#   4. R包或历史模板流程误落在项目根目录的临时/冗余目录
# 不删除 scripts、data、原始输入文件、R函数脚本或构建脚本本身。
clean_generated_outputs() {
  rm -rf \
    results \
    temporary \
    manuscripts \
    omnipathr-log \
    omipathr-log \
    scripts/beamer/beamer_report.tex \
    scripts/beamer/sections \
    scripts/beamer/generated_tables

  find . -path "./.git" -prune -o -type d -name "__pycache__" -exec rm -rf {} +
  find . -path "./.git" -prune -o -name ".DS_Store" -delete

  mkdir -p results temporary
}

LOG_FILE="$(mktemp -t crlm_clean_generated_outputs.XXXXXX.log)"
print_step_header "Clean Generated Outputs"
if run_with_spinner "Removing regenerable outputs and temporary files" "$LOG_FILE" clean_generated_outputs; then
  rm -f "$LOG_FILE"
fi
