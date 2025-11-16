#!/bin/bash

# 从蛋白比对开始，一次性跑完 PAL2NAL -> 参考 AXT -> KaKs -> 统计 -> 绘图

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_PAL2NAL="${SCRIPT_DIR}/5-2-蛋白CDS转密码子比对.py"
PY_AXT="${SCRIPT_DIR}/6-1-构建KaKs参考AXT.py"
PY_STATS="${SCRIPT_DIR}/7-1-统计KaKs结果.py"
PY_PLOT="${SCRIPT_DIR}/7-2-绘制KaKs分布.py"

# ==================== 用户可修改区域 ====================
# 输入/输出目录
ALIGNMENT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/alignments/protein"
CDS_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/sequences/cds"
CODON_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/alignments/codon"
AXT_ALL_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/alignments/axt_all"
AXT_REF_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/axt"
KAKS_OUTPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/kaks"
KAKS_STATS_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/kaks_stats"
KAKS_PLOT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/LolD_pipeline/kaks_plots"

# Query & 参考
QUERIES=("lolD.aln")
REFERENCE="GCA_000351625.1_Esch_coli_KTE66_V1_genomic|KB732762.1_1373|score=384.5|cov=1.04"

# PAL2NAL
PAL2NAL_BIN="pal2nal.pl"   # 若不在 PATH，请填绝对路径
CODON_TABLE=11
PAL2NAL_EXTRA=()           # 示例: ("-nogap")

# KaKs_Calculator
KAKS_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/4-dNdS-KaKs/1-KaKs_Calculator3/src"
KAKS_BIN="${KAKS_DIR}/KaKs"

# KaKs 绘图
KAKS_METHOD=""  # 可选：指定 Method 名称
# ==================== 用户可修改区域 ====================

run_cmd() {
    echo "▶️  $*"
    "$@"
}

if [[ ! -x "$KAKS_BIN" ]]; then
    echo "❌ 未找到 KaKs 可执行文件: $KAKS_BIN"
    exit 1
fi

# Step 1: PAL2NAL
cmd_pal=(python3 "$PY_PAL2NAL" --alignment-dir "$ALIGNMENT_DIR" --cds-dir "$CDS_DIR" --output-dir "$CODON_DIR")
[[ -n "$PAL2NAL_BIN" ]] && cmd_pal+=(--pal2nal-bin "$PAL2NAL_BIN")
cmd_pal+=(--codon-table "$CODON_TABLE")
if [[ ${#PAL2NAL_EXTRA[@]} -gt 0 ]]; then
    cmd_pal+=(--pal2nal-extra "${PAL2NAL_EXTRA[@]}")
fi
[[ -n "$AXT_ALL_DIR" ]] && cmd_pal+=(--axt-dir "$AXT_ALL_DIR")
if [[ ${#QUERIES[@]} -gt 0 ]]; then
    cmd_pal+=(--queries "${QUERIES[@]}")
fi
run_cmd "${cmd_pal[@]}"

# Step 2: 构建参考 AXT
cmd_axt=(python3 "$PY_AXT" --alignment-dir "$CODON_DIR" --output-dir "$AXT_REF_DIR")
[[ -n "$REFERENCE" ]] && cmd_axt+=(--reference "$REFERENCE")
if [[ ${#QUERIES[@]} -gt 0 ]]; then
    cmd_axt+=(--queries "${QUERIES[@]}")
fi
run_cmd "${cmd_axt[@]}"

# Step 3: KaKs
mkdir -p "$KAKS_OUTPUT_DIR"
shopt -s nullglob
for axt in "$AXT_REF_DIR"/*.axt; do
    filename=$(basename "$axt")
    base="${filename%.axt}"
    output_file="$KAKS_OUTPUT_DIR/${base}.kaks.tsv"
    echo "▶️  KaKs_Calculator: $filename"
    "$KAKS_BIN" -i "$axt" -o "$output_file"
done
shopt -u nullglob

# Step 4: 统计
cmd_stats=(python3 "$PY_STATS" --kaks-dir "$KAKS_OUTPUT_DIR" --output-dir "$KAKS_STATS_DIR")
run_cmd "${cmd_stats[@]}"

# Step 5: 绘图
SUMMARY_FILE="$KAKS_STATS_DIR/KaKs_summary.tsv"
cmd_plot=(python3 "$PY_PLOT" --summary "$SUMMARY_FILE" --output-dir "$KAKS_PLOT_DIR")
[[ -n "$KAKS_METHOD" ]] && cmd_plot+=(--method "$KAKS_METHOD")
run_cmd "${cmd_plot[@]}"

echo "✓ 下游 KaKs 流水线完成"
