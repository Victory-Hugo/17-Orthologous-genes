#!/bin/bash

################################################################################
# 脚本名称: 5-下游KaKs流水线.sh
# 作者: BigLin
#
# 功能: 从蛋白 MSA 起步依次执行 PAL2NAL -> 参考 AXT -> KaKs -> 统计 -> 绘图。
# 依赖软件:
#   1. Python 3.8+（用于调用 5-2/6-1/7-1/7-2 等脚本）
#   2. PAL2NAL (pal2nal.pl) - https://www.bork.embl.de/pal2nal/
#   3. KaKs_Calculator 3.0 - 已编译的 KaKs 可执行文件
#   4. MAFFT（由 5-1-蛋白多序列比对.sh 负责）
#   5. Python 包 pandas / matplotlib（随 7-1、7-2 使用）
################################################################################

set -euo pipefail

log() {
    local level="$1"; shift
    local timestamp
    timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$timestamp][$level] $*"
}
log_info() { log "INFO" "$*"; }
log_warn() { log "WARN" "$*"; }
log_error() { log "ERROR" "$*"; }

trap 'log_error "流水线在行 $LINENO 被中断"; exit 1' ERR

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
    log_info "执行命令: $*"
    "$@"
}

log_info "检查输入目录与依赖"
for path in "$ALIGNMENT_DIR" "$CDS_DIR"; do
    if [[ ! -d "$path" ]]; then
        log_error "必需目录缺失: $path"
        exit 1
    fi
done

if [[ ! -x "$KAKS_BIN" ]]; then
    log_error "未找到 KaKs 可执行文件: $KAKS_BIN"
    exit 1
fi

log_info "Query: ${QUERIES[*]:-全部} | 参考: ${REFERENCE:-默认首个样本}"

# Step 1: PAL2NAL
log_info "步骤1: PAL2NAL 生成密码子对齐"
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
log_info "步骤2: 根据密码子比对构建参考 AXT"
cmd_axt=(python3 "$PY_AXT" --alignment-dir "$CODON_DIR" --output-dir "$AXT_REF_DIR")
[[ -n "$REFERENCE" ]] && cmd_axt+=(--reference "$REFERENCE")
if [[ ${#QUERIES[@]} -gt 0 ]]; then
    cmd_axt+=(--queries "${QUERIES[@]}")
fi
run_cmd "${cmd_axt[@]}"

# Step 3: KaKs
log_info "步骤3: 调用 KaKs_Calculator"
mkdir -p "$KAKS_OUTPUT_DIR"
shopt -s nullglob
count_kaks=0
for axt in "$AXT_REF_DIR"/*.axt; do
    filename=$(basename "$axt")
    base="${filename%.axt}"
    output_file="$KAKS_OUTPUT_DIR/${base}.kaks.tsv"
    log_info "KaKs_Calculator 输入: $filename"
    if ! "$KAKS_BIN" -i "$axt" -o "$output_file"; then
        log_error "KaKs_Calculator 处理 $filename 失败"
        exit 1
    fi
    count_kaks=$((count_kaks + 1))
done
shopt -u nullglob
if [[ $count_kaks -eq 0 ]]; then
    log_error "在 $AXT_REF_DIR 中未找到任何 .axt 文件，无法继续"
    exit 1
fi
log_info "KaKs_Calculator 完成：生成 ${count_kaks} 个结果文件"

# Step 4: 统计
log_info "步骤4: 统计 Ka/Ks 结果"
cmd_stats=(python3 "$PY_STATS" --kaks-dir "$KAKS_OUTPUT_DIR" --output-dir "$KAKS_STATS_DIR")
run_cmd "${cmd_stats[@]}"

# Step 5: 绘图
log_info "步骤5: 绘制 Ka/Ks 分布图"
SUMMARY_FILE="$KAKS_STATS_DIR/KaKs_summary.tsv"
if [[ ! -s "$SUMMARY_FILE" ]]; then
    log_error "未找到 KaKs_summary.tsv，统计步骤可能失败"
    exit 1
fi
cmd_plot=(python3 "$PY_PLOT" --summary "$SUMMARY_FILE" --output-dir "$KAKS_PLOT_DIR")
[[ -n "$KAKS_METHOD" ]] && cmd_plot+=(--method "$KAKS_METHOD")
run_cmd "${cmd_plot[@]}"

log_info "全部步骤完成。结果目录："
log_info "  密码子对齐: $CODON_DIR"
log_info "  all-pairs AXT: ${AXT_ALL_DIR:-已禁用}"
log_info "  参考 AXT: $AXT_REF_DIR"
log_info "  KaKs 表格: $KAKS_OUTPUT_DIR"
log_info "  统计: $KAKS_STATS_DIR"
log_info "  图像: $KAKS_PLOT_DIR"
log_info "✓ 下游 KaKs 流水线完成。"
