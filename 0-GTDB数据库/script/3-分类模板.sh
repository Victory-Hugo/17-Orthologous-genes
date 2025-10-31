#!/bin/bash

# GTDB-Tk 分类运行模板
# 使用方法: bash 3-分类模板.sh <基因组目录> [输出目录] [CPU数量]

# 设置错误处理
set -e

# 脚本路径
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 参数设置
GENOME_DIR="${1:-.}"
OUT_DIR="${2:-./gtdbtk_output}"
CPUS="${3:-4}"

echo "========================================"
echo "GTDB-Tk 基因组分类"
echo "========================================"
echo ""
echo "配置信息:"
echo "  基因组目录: $GENOME_DIR"
echo "  输出目录: $OUT_DIR"
echo "  CPU数量: $CPUS"
echo ""

# 检查基因组目录
if [ ! -d "$GENOME_DIR" ]; then
    echo "❌ 错误: 基因组目录不存在: $GENOME_DIR"
    echo "用法: bash 3-分类模板.sh <基因组目录> [输出目录] [CPU数量]"
    echo "示例: bash 3-分类模板.sh ./genomes ./results 8"
    exit 1
fi

# 检查基因组文件
GENOME_COUNT=$(find "$GENOME_DIR" -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" \) | wc -l)

if [ $GENOME_COUNT -eq 0 ]; then
    echo "❌ 错误: 在 $GENOME_DIR 中未找到基因组文件"
    echo "支持的格式: *.fasta, *.fa, *.fna"
    exit 1
fi

echo "✅ 找到 $GENOME_COUNT 个基因组文件"
echo ""

# 检查 GTDBTK_DATA_PATH
if [ -z "$GTDBTK_DATA_PATH" ]; then
    echo "⚠️  环境变量 GTDBTK_DATA_PATH 未设置"
    echo "请先运行: source $SCRIPT_DIR/set_gtdbtk_env.sh"
    exit 1
fi

echo "✅ GTDBTK 数据路径: $GTDBTK_DATA_PATH"
echo ""

# 创建输出目录
mkdir -p "$OUT_DIR"

# 运行分类
echo "📌 开始分类流程..."
echo "---"
echo "时间: $(date)"
echo ""

gtdbtk classify_wf \
    --genome_dir "$GENOME_DIR" \
    --out_dir "$OUT_DIR" \
    --cpus "$CPUS" \
    --extension fasta \
    --force

echo ""
echo "---"
echo "✅ 分类完成！"
echo "时间: $(date)"
echo ""

# 显示结果摘要
echo "📊 结果文件:"
echo "---"

if [ -f "$OUT_DIR/classify/classify_wf.summary.tsv" ]; then
    echo ""
    echo "分类摘要 ($OUT_DIR/classify/classify_wf.summary.tsv):"
    head -n 5 "$OUT_DIR/classify/classify_wf.summary.tsv" | column -t -s $'\t'
    echo "..."
fi

if [ -f "$OUT_DIR/classify/classify_wf.markers_summary.tsv" ]; then
    echo ""
    echo "标记基因摘要 ($OUT_DIR/classify/classify_wf.markers_summary.tsv):"
    head -n 5 "$OUT_DIR/classify/classify_wf.markers_summary.tsv" | column -t -s $'\t'
    echo "..."
fi

echo ""
echo "📂 完整输出目录结构:"
tree "$OUT_DIR" -L 2 2>/dev/null || find "$OUT_DIR" -type f | head -20

echo ""
echo "========================================"
echo "✅ 所有操作完成"
echo "========================================"
