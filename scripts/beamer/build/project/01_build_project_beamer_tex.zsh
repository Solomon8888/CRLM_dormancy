#!/usr/bin/env zsh
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
cd "$PROJECT_ROOT"

# 项目级Beamer构建入口。
# 后续新增数据集时，在这里继续串联各数据集的02号构建脚本即可。
zsh scripts/beamer/build/GSE114012/02_build_GSE114012_beamer_tex.zsh
