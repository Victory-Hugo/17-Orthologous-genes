#!/bin/bash

# 使用 MAFFT 对提取的蛋白序列进行多序列比对

set -euo pipefail

INPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/sequences/protein"
OUTPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/alignments/protein"
THREADS=8
MAFFT_BIN="$(command -v mafft)"

if [[ -z "$MAFFT_BIN" ]]; then
    echo "❌ 未找到 mafft，请先安装"
    exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "❌ 输入目录不存在: $INPUT_DIR"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

for fasta in "$INPUT_DIR"/*.faa; do
    [[ -e "$fasta" ]] || continue
    filename=$(basename "$fasta")
    base="${filename%.faa}"
    output_file="$OUTPUT_DIR/${base}.aln.faa"
    echo "▶️  MAFFT 比对: $filename -> $(basename "$output_file")"
    "$MAFFT_BIN" --thread "$THREADS" --maxiterate 1000 --localpair "$fasta" > "$output_file"
done

echo "✓ 所有比对完成，输出目录: $OUTPUT_DIR"
