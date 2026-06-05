#!/usr/bin/env zsh
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
cd "$PROJECT_ROOT"

# 仅清理可再生成内容：
#   1. results 下全部分析结果、图片和最终报告
#   2. temporary 下全部中间文件、日志、LaTeX构建文件和临时表格
#   3. Beamer构建脚本自动生成的tex源文件
# 不删除 scripts、data、原始输入文件、R函数脚本或构建脚本本身。
rm -rf \
  results \
  temporary \
  scripts/beamer/beamer_report.tex \
  scripts/beamer/sections \
  scripts/beamer/generated_tables

find . -path "./.git" -prune -o -type d -name "__pycache__" -exec rm -rf {} +
find . -path "./.git" -prune -o -name ".DS_Store" -delete

mkdir -p results temporary
