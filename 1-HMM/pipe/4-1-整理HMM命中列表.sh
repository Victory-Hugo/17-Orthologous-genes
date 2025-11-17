#!/bin/bash

# 整理 *.hits.faa 为 TSV 清单

set -euo pipefail



# ==================== 用户可修改区域 ====================
HITS_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output"
OUTPUT_TSV="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/hits_manifest.tsv"
PYTHON_SCRIPT="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/script/4-1-整理HMM命中列表.py"

QUERIES=("lolD.aln")   # 留空 () 表示处理全部 query
# ==================== 用户可修改区域 ====================

cmd=(python3 "$PYTHON_SCRIPT" --hits-dir "$HITS_DIR" --output-tsv "$OUTPUT_TSV")
if [[ ${#QUERIES[@]} -gt 0 ]]; then
    cmd+=(--queries "${QUERIES[@]}")
fi

echo "运行命令: ${cmd[*]}"
"${cmd[@]}"
