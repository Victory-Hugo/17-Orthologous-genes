#!/bin/bash

# 构建 KaKs_Calculator 需要的 AXT 文件

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="${SCRIPT_DIR}/6-1-构建KaKs参考AXT.py"

# ==================== 用户可修改区域 ====================
ALIGNMENT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/alignments/codon"
OUTPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/axt"
REFERENCE="GCA_000351625.1_Esch_coli_KTE66_V1_genomic|KB732762.1_1373|score=384.5|cov=1.04"
QUERIES=("lolD.trim")
# ==================== 用户可修改区域 ====================

cmd=(python3 "$PYTHON_SCRIPT" --alignment-dir "$ALIGNMENT_DIR" --output-dir "$OUTPUT_DIR")
[[ -n "$REFERENCE" ]] && cmd+=(--reference "$REFERENCE")
if [[ ${#QUERIES[@]} -gt 0 ]]; then
    cmd+=(--queries "${QUERIES[@]}")
fi

echo "运行命令: ${cmd[*]}"
"${cmd[@]}"
