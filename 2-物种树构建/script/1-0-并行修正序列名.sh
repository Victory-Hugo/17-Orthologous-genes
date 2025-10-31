#!/bin/bash
# 并行修正核心基因序列名脚本
# 使用GNU Parallel加速处理
# 支持多序列文件自动编号
# 作者: 你未来的自己
# 日期: 2025-10-31

# ===== 可配置区域 =====
EXT="${EXT:-.faa}"  # 文件后缀（可通过命令行 EXT=.fna ./script.sh 指定）
CORE_GENES_DIR="/mnt/f/15_Bam_Tam/2-物种树/output/all_core_genomes" #* 将该文件夹下（子文件夹）中的所有序列修改
LOG_FILE="/mnt/f/15_Bam_Tam/2-物种树/output/sequence_rename_log.txt"
THREADS=$(nproc)
# =====================

echo "=========================================="
echo "并行修正核心基因序列名脚本"
echo "日期: $(date)"
echo "CPU核心数: $THREADS"
echo "目标后缀: $EXT"
echo "=========================================="

> "$LOG_FILE"

echo "正在统计文件数量..."
TOTAL_FILES=$(find "$CORE_GENES_DIR" -name "*${EXT}" | wc -l)
echo "需要处理的文件总数: $TOTAL_FILES"

process_faa_file() {
    faa_file="$1"
    log_file="$2"
    filename=$(basename "$faa_file" "$EXT")
    temp_file="${faa_file}.tmp.$$"

    seq_count=$(grep -c '^>' "$faa_file")
    if [[ "$seq_count" -eq 0 ]]; then
        echo "ERROR: $filename (无序列)" >&2
        return 1
    fi

    # 修改点：删除注释，仅保留'>'后的第一个字段作为原序列名基础
    awk -v prefix="$filename" -v log="$log_file" -v file="$faa_file" '
        BEGIN { seqnum=0 }
        /^>/ {
            seqnum++
            # 去掉注释，只保留>后面的第一个字段
            split($0, arr, " ")
            oldname = arr[1]
            newname = ">" prefix "_" seqnum
            print newname
            print file ": " oldname " -> " newname >> log
            next
        }
        { print }
    ' "$faa_file" > "$temp_file"

    if [[ -s "$temp_file" ]]; then
        mv "$temp_file" "$faa_file"
        echo "SUCCESS: $filename ($seq_count 序列)"
    else
        echo "ERROR: $filename (生成失败)" >&2
        rm -f "$temp_file"
        return 1
    fi
}

export -f process_faa_file
export LOG_FILE
export EXT

echo ""
echo "开始并行处理所有${EXT}文件..."
start_time=$(date +%s)

find "$CORE_GENES_DIR" -name "*${EXT}" | \
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
