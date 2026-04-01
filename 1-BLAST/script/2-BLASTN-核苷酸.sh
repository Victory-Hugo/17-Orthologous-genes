#!/usr/bin/env bash
set -euo pipefail

# ===============================================
# 使用 BLASTN 搜索 Rv1819c 核苷酸同源序列
# 用法: ./2-BLASTN-核苷酸.sh [身份阈值] [覆盖率阈值] [最小比对长度]
# ===============================================

# ---------- 参数配置 ----------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONF_FILE="/mnt/f/onedrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-BLAST/conf/Rv1819c/2-BLASTN-核苷酸.conf"
[[ -f "$CONF_FILE" ]] || { echo "错误: 配置文件不存在: $CONF_FILE" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONF_FILE"

IDENTITY_THRESHOLD=${1:-$IDENTITY_THRESHOLD}
COVERAGE_THRESHOLD=${2:-$COVERAGE_THRESHOLD}
LENGTH_THRESHOLD=${3:-$LENGTH_THRESHOLD}

# ---------- 路径配置 ----------
# 所有路径变量由配置文件提供

resolve_target_fna() {
    local input_path="$1"
    local candidate=""

    if [[ -f "$input_path" && "$input_path" == *.fna ]]; then
        printf '%s\n' "$input_path"
        return 0
    fi

    if [[ "$input_path" == *.faa ]]; then
        candidate="${input_path%.faa}.fna"
    else
        candidate="${input_path%.*}.fna"
    fi

    if [[ -f "$candidate" ]]; then
        printf '%s\n' "$candidate"
        return 0
    fi

    return 1
}

# ---------- 前置检查 ----------
command -v blastn >/dev/null 2>&1 || { echo "错误: 未找到 blastn，请先安装 BLAST+。" >&2; exit 1; }
[[ -f "$QUERY_FASTA" ]] || { echo "错误: 查询序列不存在: $QUERY_FASTA" >&2; exit 1; }
[[ -f "$TARGET_LIST" ]] || { echo "错误: 目标列表不存在: $TARGET_LIST" >&2; exit 1; }

mkdir -p "$OUTPUT_DIR" "$SEQ_DIR"

# ---------- 断点续跑准备 ----------
declare -A DONE
if [[ -s "$RAW_RESULT" ]]; then
    echo "▶ 检测到已有结果，启用断点续跑..."
    while IFS=$'\t' read -r src _; do
        [[ -z "$src" || "$src" == "source_file" ]] && continue
        DONE["$src"]=1
    done < "$RAW_RESULT"
else
    echo -e "source_file\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen" > "$RAW_RESULT"
fi

echo "▶ 开始 BLASTN 搜索..."
while IFS= read -r target_path; do
    [[ -z "$target_path" ]] && continue

    if ! resolved_fna="$(resolve_target_fna "$target_path")"; then
        echo "⚠ 跳过，未找到对应核苷酸文件: $target_path" >&2
        continue
    fi

    filename="$(basename "$resolved_fna")"
    if [[ -n "${DONE[$filename]:-}" ]]; then
        echo "  ↩︎ 跳过（已完成）: $filename"
        continue
    fi

    blastn \
        -query "$QUERY_FASTA" \
        -subject "$resolved_fna" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
        -evalue "$EVALUE" \
        -num_threads "$THREADS" \
        | awk -v file="$filename" -v minlen="$LENGTH_THRESHOLD" 'BEGIN{OFS="\t"} $4>=minlen {print file, $0}' >> "$RAW_RESULT"

    echo "  ✓ 已处理: $filename"
done < "$TARGET_LIST"

echo "▶ 过滤结果..."
python3 "$PY_FILTER" "$RAW_RESULT" "$FILTERED_RESULT" "$IDENTITY_THRESHOLD" "$COVERAGE_THRESHOLD"
echo "  ✓ 过滤结果已保存: $FILTERED_RESULT"

echo "▶ 提取匹配序列..."
python3 "$PY_EXTRACT" "$FILTERED_RESULT" "$TARGET_LIST" "$SEQ_DIR"
echo "  ✓ 序列已保存到: $SEQ_DIR"

echo ""
echo "✓ 全部流程完成。"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "过滤条件:"
echo "  • 身份相似度: ≥ ${IDENTITY_THRESHOLD}%"
echo "  • 查询覆盖率: ≥ ${COVERAGE_THRESHOLD}%"
echo "  • 比对长度: ≥ ${LENGTH_THRESHOLD}"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "输出文件:"
echo "  • 原始结果: $RAW_RESULT"
echo "  • 过滤结果: $FILTERED_RESULT"
echo "  • 匹配序列: $SEQ_DIR"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
