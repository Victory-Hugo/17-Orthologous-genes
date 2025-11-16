#!/bin/bash

# 根据 manifest 提取蛋白 + CDS 序列

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="${SCRIPT_DIR}/4-2-提取命中CDS蛋白序列.py"

# ==================== 用户可修改区域 ====================
MANIFEST="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/hits_manifest.tsv"
SAMPLE_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/data/3-样本序列"
OUTPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/sequences"
QUERIES=("lolD.trim")
# ==================== 用户可修改区域 ====================

cmd=(python3 "$PYTHON_SCRIPT" --manifest "$MANIFEST" --sample-dir "$SAMPLE_DIR" --output-dir "$OUTPUT_DIR")
if [[ ${#QUERIES[@]} -gt 0 ]]; then
    cmd+=(--queries "${QUERIES[@]}")
fi

echo "运行命令: ${cmd[*]}"
"${cmd[@]}"
