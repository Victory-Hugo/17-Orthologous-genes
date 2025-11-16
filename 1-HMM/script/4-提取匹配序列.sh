#!/bin/bash

# 单个样本提取 LolD/LolF 等基因的 HMM 命中序列

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="${SCRIPT_DIR}/4-提取匹配序列.py"

# ==================== 用户可修改区域 ====================
FAA_FILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/data/3-样本序列/GCA_000351625.1_Esch_coli_KTE66_V1_genomic.faa"
TBL_FILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/tbl_merge.tsv"
OUTPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output"
MAX_EVALUE="1e-20"
MIN_COVERAGE="0.6"
MIN_LENGTH="180"
EXTRACT_MODE="only_one"   # 可选: all / only_one
MODEL_LENGTH="225"        # HMM LENG，用于覆盖度
# ==================== 用户可修改区域 ====================

cmd=(python3 "$PYTHON_SCRIPT" "$FAA_FILE" "$TBL_FILE" "$OUTPUT_DIR")
[[ -n "$MAX_EVALUE" ]] && cmd+=(--max-evalue "$MAX_EVALUE")
[[ -n "$MIN_COVERAGE" ]] && cmd+=(--min-coverage "$MIN_COVERAGE")
[[ -n "$MIN_LENGTH" ]] && cmd+=(--min-length "$MIN_LENGTH")
[[ -n "$EXTRACT_MODE" ]] && cmd+=(--mode "$EXTRACT_MODE")
[[ -n "$MODEL_LENGTH" ]] && cmd+=(--model-length "$MODEL_LENGTH")

echo "运行命令: ${cmd[*]}"
"${cmd[@]}"
