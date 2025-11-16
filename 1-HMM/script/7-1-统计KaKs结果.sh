#!/bin/bash

# 汇总 Ka/Ks 结果

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="${SCRIPT_DIR}/7-1-统计KaKs结果.py"

# ==================== 用户可修改区域 ====================
KAKS_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/kaks"
OUTPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/kaks_stats"
# ==================== 用户可修改区域 ====================

cmd=(python3 "$PYTHON_SCRIPT" --kaks-dir "$KAKS_DIR" --output-dir "$OUTPUT_DIR")

echo "运行命令: ${cmd[*]}"
"${cmd[@]}"
