#!/bin/bash
set -e
#===============================================
# BLASTP 同源基因搜索脚本
# 功能: 在多个蛋白质序列数据库中搜索目标基因
# 用法: ./1-BLASTP.sh [身份阈值] [覆盖率阈值]
#===============================================
# =============== 参数配置 ===============
# 接收命令行参数，使用默认值作为后备
IDENTITY_THRESHOLD=${1:-40}    # 最小身份相似度 (%)
COVERAGE_THRESHOLD=${2:-70}    # 最小查询覆盖率 (%)
EVALUE=1e-5                    # E值阈值
# =============== 文件路径配置 ===============
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FAA_PATH_LIST="/mnt/l/18-Rv0194-Gene/1-BLAST/meta/Actinobacteria_Mycobacterium_faa_path.txt"
QUERY_FASTA="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-BLAST/conf/Rv0194.aln.faa"
OUTPUT_DIR="/mnt/l/18-Rv0194-Gene/1-BLAST/data/search/"
OUT_FILE="$OUTPUT_DIR/Rv0194_blastp.tsv"
OUT_SEQ_DIR="$OUTPUT_DIR/sequence"
PY_FILTER="$SCRIPT_DIR/filter_blastp_results.py"
PY_EXTRACT="$SCRIPT_DIR/extract_sequences.py"
PYTHON="python3"
# =============== 初始化工作目录 ===============
mkdir -p "$OUTPUT_DIR"

TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

# =============== 运行 BLASTP ===============
echo "▶ 正在执行 BLASTP 搜索..."

all_results="$TEMP_DIR/all_blastp_results.tsv"
echo -e "source_file\traw_results" > "$all_results"

while IFS= read -r faa_file; do
    [ -z "$faa_file" ] && continue
    
    if [ ! -f "$faa_file" ]; then
        echo "⚠ 警告: 文件不存在: $faa_file" >&2
        continue
    fi
    
    filename=$(basename "$faa_file")
    blast_db="$TEMP_DIR/blast_db_$filename"
    
    # 创建 BLAST 数据库
    makeblastdb -in "$faa_file" \
               -dbtype prot \
               -out "$blast_db" \
               -logfile "$TEMP_DIR/makeblastdb_$filename.log" 2>/dev/null
    
    # 执行 BLASTP 搜索
    blastp -query "$QUERY_FASTA" \
           -db "$blast_db" \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
           -evalue "$EVALUE" 2>/dev/null | while read -r line; do
        echo "$filename	$line" >> "$all_results"
    done
    
    echo "  ✓ 已处理: $filename"
    
done < "$FAA_PATH_LIST"

# =============== 过滤结果 ===============
echo "▶ 正在过滤搜索结果..."

$PYTHON "$PY_FILTER" \
    "$all_results" \
    "$OUT_FILE" \
    "$IDENTITY_THRESHOLD" \
    "$COVERAGE_THRESHOLD"

echo "  ✓ 结果已保存: $OUT_FILE"

# =============== 提取序列 ===============
echo "▶ 正在提取匹配的序列..."

$PYTHON "$PY_EXTRACT" \
    "$OUT_FILE" \
    "$FAA_PATH_LIST" \
    "$OUT_SEQ_DIR"

echo "  ✓ 序列已保存: $OUT_SEQ_DIR"

# =============== 完成 ===============
echo ""
echo "✓ 流程完成!"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "过滤条件:"
echo "  • 身份相似度: ≥ ${IDENTITY_THRESHOLD}%"
echo "  • 查询覆盖率: ≥ ${COVERAGE_THRESHOLD}%"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "输出文件:"
echo "  • BLASTP 结果: $OUT_FILE"
echo "  • 提取序列: $OUT_SEQ_DIR"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
