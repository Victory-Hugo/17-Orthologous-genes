#!/bin/bash

# 根据 manifest 提取蛋白 + CDS 序列

set -euo pipefail



# ==================== 用户可修改区域 ====================
SAMPLE_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/data/3-样本序列" #! 原始的序列文件目录
PIPEDIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/"
MANIFEST="${PIPEDIR}/hits_manifest.tsv" #! 由4-1-整理HMM命中列表.sh
OUTPUT_DIR="${PIPEDIR}/sequences"
QUERIES=("lolD.aln")
PYTHON_SCRIPT="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/script/4-2-提取命中CDS蛋白序列.py"

# ==================== 用户可修改区域 ====================

cmd=(python3 "$PYTHON_SCRIPT" --manifest "$MANIFEST" --sample-dir "$SAMPLE_DIR" --output-dir "$OUTPUT_DIR")
if [[ ${#QUERIES[@]} -gt 0 ]]; then
    cmd+=(--queries "${QUERIES[@]}")
fi

echo "运行命令: ${cmd[*]}"
"${cmd[@]}"
