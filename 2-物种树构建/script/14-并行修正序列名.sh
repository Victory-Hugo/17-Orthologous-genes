#!/bin/bash

# 并行修正核心基因序列名脚本
# 使用GNU Parallel加速处理
# 日期: 2025-10-29

CORE_GENES_DIR="/mnt/f/15_Bam_Tam/2-物种树/output/all_core_genomes"
LOG_FILE="/mnt/f/15_Bam_Tam/2-物种树/output/sequence_rename_log.txt"
THREADS=$(nproc)  # 使用所有CPU核心

echo "=========================================="
echo "并行修正核心基因序列名脚本"
echo "日期: $(date)"
echo "使用CPU核心数: $THREADS"
echo "=========================================="

# 创建日志文件
> "$LOG_FILE"

# 统计总文件数
echo "正在统计文件数量..."
TOTAL_FILES=$(find "$CORE_GENES_DIR" -name "*.faa" | wc -l)
echo "需要处理的文件总数: $TOTAL_FILES"

# 定义处理单个文件的函数
process_faa_file() {
    faa_file="$1"
    log_file="$2"
    
    # 获取文件名（不包含.faa后缀）
    filename=$(basename "$faa_file" .faa)
    
    # 读取原始序列名
    original_seqname=$(head -1 "$faa_file" | sed 's/^>//')
    
    # 创建临时文件
    temp_file="${faa_file}.tmp.$$"
    
    # 修改序列名为文件名，保留原始注释
    if [[ "$original_seqname" =~ ^([A-Z0-9_]+\.[0-9]+)(.*)$ ]]; then
        # 提取原始注释部分
        annotation="${BASH_REMATCH[2]}"
        # 创建新的序列名：文件名 + 原始注释
        new_seqname="${filename}${annotation}"
    else
        # 如果格式不匹配，直接使用文件名
        new_seqname="$filename"
    fi
    
    # 写入新文件
    {
        echo ">$new_seqname"
        tail -n +2 "$faa_file"
    } > "$temp_file"
    
    # 替换原文件
    if [[ -s "$temp_file" ]]; then
        mv "$temp_file" "$faa_file"
        echo "$faa_file: $original_seqname -> $new_seqname" >> "$log_file"
        echo "SUCCESS: $filename"
    else
        echo "ERROR: $filename" >&2
        rm -f "$temp_file"
        return 1
    fi
}

# 导出函数和变量
export -f process_faa_file
export LOG_FILE

echo ""
echo "开始并行处理所有.faa文件..."
start_time=$(date +%s)

# 使用parallel处理所有.faa文件
find "$CORE_GENES_DIR" -name "*.faa" | \
parallel -j "$THREADS" --progress process_faa_file {} "$LOG_FILE"

end_time=$(date +%s)
duration=$((end_time - start_time))

echo ""
echo "=========================================="
echo "并行序列名修正完成！"
echo "总耗时: ${duration}秒"
echo "处理文件数: $TOTAL_FILES"
echo "平均速度: $(echo "scale=2; $TOTAL_FILES / $duration" | bc) 文件/秒"
echo "日志文件: $LOG_FILE"
echo "=========================================="

# 统计处理结果
success_count=$(grep -c "SUCCESS:" "$LOG_FILE" 2>/dev/null || echo "0")
error_count=$(grep -c "ERROR:" "$LOG_FILE" 2>/dev/null || echo "0")

echo ""
echo "处理统计:"
echo "- 成功: $success_count 个文件"
echo "- 失败: $error_count 个文件"

# 验证修正结果
echo ""
echo "验证修正结果（随机检查10个文件）:"
find "$CORE_GENES_DIR" -name "*.faa" | shuf | head -10 | while read file; do
    filename=$(basename "$file" .faa)
    seqname=$(head -1 "$file" | sed 's/^>//; s/ .*//')
    if [[ "$filename" == "$seqname" ]]; then
        echo "✓ $filename: 序列名匹配"
    else
        echo "✗ $filename: 序列名不匹配 ($seqname)"
    fi
done

echo ""
echo "序列名修正完成！可以继续执行并行多序列比对。"