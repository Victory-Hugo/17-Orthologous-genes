#!/bin/bash

# 蛋白比对 + CDS -> 密码子比对 + 可选 AXT

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="${SCRIPT_DIR}/5-2-蛋白CDS转密码子比对.py"

# ==================== 用户可修改区域 ====================
ALIGNMENT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/alignments/protein"
CDS_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/sequences/cds"
OUTPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/alignments/codon"
AXT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/alignments/axt_all"
QUERIES=("lolD.trim")
PAL2NAL_BIN="pal2nal.pl"     # 如果不在 PATH，请写绝对路径
CODON_TABLE=11               # 细菌通常为 11
PAL2NAL_EXTRA=()             # 例如 ("-nogap")
# ==================== 用户可修改区域 ====================

cmd=(python3 "$PYTHON_SCRIPT" --alignment-dir "$ALIGNMENT_DIR" --cds-dir "$CDS_DIR" --output-dir "$OUTPUT_DIR")
[[ -n "$PAL2NAL_BIN" ]] && cmd+=(--pal2nal-bin "$PAL2NAL_BIN")
cmd+=(--codon-table "$CODON_TABLE")
if [[ ${#PAL2NAL_EXTRA[@]} -gt 0 ]]; then
    cmd+=(--pal2nal-extra "${PAL2NAL_EXTRA[@]}")
fi
[[ -n "$AXT_DIR" ]] && cmd+=(--axt-dir "$AXT_DIR")
if [[ ${#QUERIES[@]} -gt 0 ]]; then
    cmd+=(--queries "${QUERIES[@]}")
fi

echo "运行命令: ${cmd[*]}"
"${cmd[@]}"
