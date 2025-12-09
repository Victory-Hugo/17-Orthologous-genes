#!/usr/bin/env bash
set -euo pipefail

# ===============================================
# 使用 DIAMOND 搜索 Rv0194（精确结构域）同源序列
# 用法: ./2-diamond.sh [身份阈值] [覆盖率阈值] [最小比对长度]
# ===============================================

# ---------- 参数配置 ----------
IDENTITY_THRESHOLD=${1:-40}    # 身份相似度阈值 (%)
COVERAGE_THRESHOLD=${2:-70}    # 查询覆盖率阈值 (%)
LENGTH_THRESHOLD=${3:-1000}    # 比对长度阈值 (alignment length)
EVALUE=1e-5                    # E-value 阈值
THREADS=${THREADS:-16}          # DIAMOND 线程数 (可通过环境变量覆盖)

# ---------- 路径配置 ----------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
QUERY_FASTA="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-BLAST/conf/Rv0194.aln.faa"
FAA_LIST="/mnt/l/18-Rv0194-Gene/1-BLAST/meta/Actinobacteria_Mycobacterium_faa_path.txt"
OUTPUT_DIR="/mnt/l/18-Rv0194-Gene/1-BLAST/data/search"
RAW_RESULT="$OUTPUT_DIR/Rv0194_diamond_raw.tsv"
FILTERED_RESULT="$OUTPUT_DIR/Rv0194_diamond.tsv"
SEQ_DIR="$OUTPUT_DIR/sequence"

PY_FILTER="$SCRIPT_DIR/filter_blastp_results.py"
PY_EXTRACT="$SCRIPT_DIR/extract_sequences.py"

# ---------- 前置检查 ----------
command -v diamond >/dev/null 2>&1 || { echo "错误: 未找到 diamond，请先安装。" >&2; exit 1; }
[[ -f "$QUERY_FASTA" ]] || { echo "错误: 查询序列不存在: $QUERY_FASTA" >&2; exit 1; }
[[ -f "$FAA_LIST" ]] || { echo "错误: FAA 列表不存在: $FAA_LIST" >&2; exit 1; }

mkdir -p "$OUTPUT_DIR" "$SEQ_DIR"

TEMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TEMP_DIR"' EXIT

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

echo "▶ 开始 DIAMOND 搜索..."
while IFS= read -r faa_file; do
    [[ -z "$faa_file" ]] && continue
    if [[ ! -f "$faa_file" ]]; then
        echo "⚠ 跳过，文件不存在: $faa_file" >&2
        continue
    fi

    filename="$(basename "$faa_file")"
    if [[ -n "${DONE[$filename]:-}" ]]; then
        echo "  ↩︎ 跳过（已完成）: $filename"
        continue
    fi

    db_prefix="$TEMP_DIR/${filename%.faa}"

    # 构建临时 DIAMOND 数据库
    diamond makedb --in "$faa_file" -d "$db_prefix" --quiet

    # 执行比对并在首列标记来源文件名
    diamond blastp \
        --query "$QUERY_FASTA" \
        --db "${db_prefix}.dmnd" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
        --evalue "$EVALUE" \
        --threads "$THREADS" \
        | awk -v file="$filename" -v minlen="$LENGTH_THRESHOLD" 'BEGIN{OFS="\t"} $4>=minlen {print file, $0}' >> "$RAW_RESULT"

    echo "  ✓ 已处理: $filename"
done < "$FAA_LIST"

echo "▶ 过滤结果..."
python3 "$PY_FILTER" "$RAW_RESULT" "$FILTERED_RESULT" "$IDENTITY_THRESHOLD" "$COVERAGE_THRESHOLD"
echo "  ✓ 过滤结果已保存: $FILTERED_RESULT"

echo "▶ 提取匹配序列..."
python3 "$PY_EXTRACT" "$FILTERED_RESULT" "$FAA_LIST" "$SEQ_DIR"
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
