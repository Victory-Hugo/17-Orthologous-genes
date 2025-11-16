#!/bin/bash

# 绘制 Ka/Ks 分布图

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="${SCRIPT_DIR}/7-2-绘制KaKs分布.py"

# ==================== 用户可修改区域 ====================
SUMMARY_FILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/kaks_stats/KaKs_summary.tsv"
OUTPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/kaks_plots"
METHOD=""   # 留空表示全部方法
# ==================== 用户可修改区域 ====================

cmd=(python3 "$PYTHON_SCRIPT" --summary "$SUMMARY_FILE" --output-dir "$OUTPUT_DIR")
[[ -n "$METHOD" ]] && cmd+=(--method "$METHOD")

echo "运行命令: ${cmd[*]}"
"${cmd[@]}"
