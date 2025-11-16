#!/bin/bash

# 使用 KaKs_Calculator3 计算 Ka/Ks

set -euo pipefail

KAKS_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/4-dNdS-KaKs/1-KaKs_Calculator3/src"
KAKS_BIN="${KAKS_DIR}/KaKs"
AXT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/axt"
OUTPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/kaks"

if [[ ! -x "$KAKS_BIN" ]]; then
    echo "❌ 未找到 KaKs 可执行文件: $KAKS_BIN"
    exit 1
fi

if [[ ! -d "$AXT_DIR" ]]; then
    echo "❌ AXT 目录不存在: $AXT_DIR"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

for axt in "$AXT_DIR"/*.axt; do
    [[ -e "$axt" ]] || continue
    filename=$(basename "$axt")
    base="${filename%.axt}"
    output_file="$OUTPUT_DIR/${base}.kaks.tsv"
    echo "▶️  KaKs: $filename -> $(basename "$output_file")"
    "$KAKS_BIN" -i "$axt" -o "$output_file"
done

echo "✓ 所有 KaKs 计算完成，输出目录: $OUTPUT_DIR"
