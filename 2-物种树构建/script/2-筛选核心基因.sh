#!/bin/bash

# 筛选核心基因


# 设置路径
INPUT_DIR="/mnt/f/15_Bam_Tam/2-物种树/output/all_single_copy"
OUTPUT_DIR="/mnt/f/15_Bam_Tam/2-物种树/output/all_core_genomes"
COVERAGE_FILE="/mnt/f/15_Bam_Tam/2-物种树/output/coverage_analysis/gene_species_count.csv"

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 覆盖率的阈值
THRESHOLD=XXX #? 这里填入的数字代表：在XXX个物种中都存在的基因，才被认为是核心基因，例如填入100，代表在本次研究中，100个物种中都存在的基因才被认为是核心基因

echo "开始筛选核心基因..."
echo "覆盖率阈值: $THRESHOLD 个物种"

# 筛选满足条件的基因
echo "正在筛选基因..."
CORE_GENES=$(awk -F',' -v threshold="$THRESHOLD" 'NR>1 && $2>=threshold {print $1}' "$COVERAGE_FILE")

# 统计核心基因数量
CORE_COUNT=$(echo "$CORE_GENES" | wc -l)
echo "找到 $CORE_COUNT 个核心基因（覆盖率≥70%）"

# 复制核心基因目录
echo "正在复制核心基因到目标目录..."
count=0
for gene in $CORE_GENES; do
    if [ -d "$INPUT_DIR/$gene" ]; then
        cp -r "$INPUT_DIR/$gene" "$OUTPUT_DIR/"
        count=$((count + 1))
        echo "已复制: $gene (进度: $count/$CORE_COUNT)"
    else
        echo "警告: 基因目录不存在 $INPUT_DIR/$gene"
    fi
done

echo "核心基因筛选完成！"
echo "源目录: $INPUT_DIR"
echo "目标目录: $OUTPUT_DIR"
echo "总基因数: $(ls -1 "$INPUT_DIR" | wc -l)"
echo "核心基因数: $count"
echo "覆盖率阈值: ≥$THRESHOLD 个物种 (≥70%)"

# 生成核心基因列表
echo "正在生成核心基因列表..."
ls -1 "$OUTPUT_DIR" > "$OUTPUT_DIR/core_genes_list.txt"
echo "核心基因列表已保存到: $OUTPUT_DIR/core_genes_list.txt"

# 生成统计报告
echo "正在生成统计报告..."
{
    echo "核心基因筛选统计报告"
    echo "====================="
    echo "筛选日期: $(date)"
    echo "总基因数: $(ls -1 "$INPUT_DIR" | wc -l)"
    echo "核心基因数: $count"
    echo "筛选比例: $(echo "scale=2; $count * 100 / $(ls -1 "$INPUT_DIR" | wc -l)" | bc)%"
    echo ""
    echo "核心基因列表:"
    echo "============"
    cat "$OUTPUT_DIR/core_genes_list.txt"
} > "$OUTPUT_DIR/core_genes_report.txt"

echo "统计报告已保存到: $OUTPUT_DIR/core_genes_report.txt"
echo "筛选完成！"