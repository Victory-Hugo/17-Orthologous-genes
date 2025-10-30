#!/bin/bash

# 脚本功能: 并发执行多个物种的基因组组装ID查询（完整合并版本）
# 用法: ./get_taxon_concurrent.sh [max_concurrent_jobs]
# 示例: ./get_taxon_concurrent.sh 5

# ===== API密钥配置 =====
export NCBI_API_KEY="29b326d54e7a21fc6c8b9afe7d71f441d809"

# 并发任务数，默认为4
MAX_CONCURRENT_JOBS="${1:-4}"

# 脚本路径配置
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TAXON_FILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-数据获取/2-下载特定物种序列/conf/taxon.txt"
OUTPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-数据获取/2-下载特定物种序列/download/"

# ===== 参数验证 =====
if [ ! -f "$TAXON_FILE" ]; then
    echo "错误: 物种文件不存在: $TAXON_FILE" >&2
    exit 1
fi

# 检查必要工具
command -v datasets >/dev/null 2>&1 || {
    echo "错误: 未找到 datasets 工具，请确保已安装 NCBI Datasets CLI" >&2
    exit 1
}

command -v jq >/dev/null 2>&1 || {
    echo "错误: 未找到 jq 工具，请安装 jq" >&2
    exit 1
}

# 创建输出目录
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/logs"

# ===== 日志文件 =====
MAIN_LOG="$OUTPUT_DIR/logs/concurrent_run_$(date +%Y%m%d_%H%M%S).log"
ERROR_LOG="$OUTPUT_DIR/logs/errors_$(date +%Y%m%d_%H%M%S).log"

# ===== 函数定义 =====

# 处理单个物种的函数（内置逻辑）
process_taxon() {
    local taxon="$1"
    local output_dir="$2"
    local log_file="$3"
    
    # 跳过空行和注释
    taxon=$(echo "$taxon" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    if [[ -z "$taxon" ]] || [[ "$taxon" =~ ^# ]]; then
        return 0
    fi
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] 开始处理物种: $taxon" | tee -a "$log_file"
    
    # ========== 内置获取物种基因组逻辑 ==========
    
    # 生成输出文件名（替换空格为下划线，移除特殊字符）
    local SAFE_NAME=$(echo "$taxon" | tr ' ' '_' | sed 's/[^a-zA-Z0-9_-]//g')
    local OUTPUT_FILE="$output_dir/${SAFE_NAME}_complete_genomes.txt"
    local SUMMARY_FILE="$output_dir/${SAFE_NAME}_summary.txt"
    
    echo "  输出文件: $OUTPUT_FILE" >> "$log_file"
    echo "  摘要文件: $SUMMARY_FILE" >> "$log_file"
    
    # 获取完整基因组数据
    local TEMP_JSON=$(mktemp)
    local TEMP_ERROR=$(mktemp)
    local MAX_RESULTS=1000
    
    # 执行 datasets 查询
    local DATASETS_SUCCESS=false
    
    echo "  执行命令: datasets summary genome taxon \"$taxon\" --limit $MAX_RESULTS" >> "$log_file"
    
    if timeout 300 datasets summary genome taxon "$taxon" --limit "$MAX_RESULTS" > "$TEMP_JSON" 2>"$TEMP_ERROR"; then
        DATASETS_SUCCESS=true
    fi
    
    if [ "$DATASETS_SUCCESS" = false ]; then
        local EXIT_CODE=$?
        if [ $EXIT_CODE -eq 124 ]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] ✗ 查询超时（5分钟）: $taxon" | tee -a "$log_file" "$ERROR_LOG"
        else
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] ✗ 无法获取数据: $taxon" | tee -a "$log_file" "$ERROR_LOG"
            cat "$TEMP_ERROR" >> "$log_file"
        fi
        rm -f "$TEMP_JSON" "$TEMP_ERROR"
        return 1
    fi
    
    rm -f "$TEMP_ERROR"
    
    # 检查是否有数据返回
    local TOTAL_COUNT=$(jq -r '.total_count // 0' "$TEMP_JSON")
    if [ "$TOTAL_COUNT" -eq 0 ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ✗ 未找到数据: $taxon" | tee -a "$log_file" "$ERROR_LOG"
        rm -f "$TEMP_JSON"
        return 1
    fi
    
    echo "  找到总计 $TOTAL_COUNT 个基因组" >> "$log_file"
    
    # 提取完整基因组的组装ID
    jq -r '.reports[] | select(.assembly_info.assembly_level == "Complete Genome") | .accession' "$TEMP_JSON" > "$OUTPUT_FILE"
    
    # 检查是否有完整基因组
    local COMPLETE_COUNT=$(wc -l < "$OUTPUT_FILE")
    if [ "$COMPLETE_COUNT" -eq 0 ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ✗ 未找到完整基因组: $taxon" | tee -a "$log_file" "$ERROR_LOG"
        rm -f "$TEMP_JSON" "$OUTPUT_FILE"
        return 1
    fi
    
    # 生成摘要信息
    {
        echo "=========================================="
        echo "基因组数据摘要 - $taxon"
        echo "查询日期: $(date '+%Y-%m-%d %H:%M:%S')"
        echo "结果限制: $MAX_RESULTS"
        echo "=========================================="
        echo ""
        echo "总基因组数量: $TOTAL_COUNT"
        echo "完整基因组数量: $COMPLETE_COUNT"
        echo ""
        echo "完整基因组菌株信息:"
        echo "------------------------------------------"
        jq -r '.reports[] | select(.assembly_info.assembly_level == "Complete Genome") | [.accession, .organism.organism_name] | @tsv' "$TEMP_JSON" | \
        sort | \
        awk -F'\t' '{printf "%-20s %s\n", $1, $2}'
        echo ""
        echo "独特菌株列表:"
        echo "------------------------------------------"
        jq -r '.reports[] | select(.assembly_info.assembly_level == "Complete Genome") | .organism.organism_name' "$TEMP_JSON" | \
        sort -u | \
        nl -w2 -s'. '
        echo ""
        echo "文件输出位置:"
        echo "- 组装ID列表: $OUTPUT_FILE"
        echo "- 摘要信息: $SUMMARY_FILE"
        echo "=========================================="
    } > "$SUMMARY_FILE"
    
    rm -f "$TEMP_JSON"
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ✓ 物种处理成功: $taxon (找到 $COMPLETE_COUNT 个完整基因组)" | tee -a "$log_file"
    return 0
}

# 等待后台任务的函数
wait_for_job() {
    local max_jobs="$1"
    while [ $(jobs -r | wc -l) -ge "$max_jobs" ]; do
        sleep 2
    done
}

# ===== 主程序 =====

echo "=========================================="
echo "开始并发查询物种基因组信息"
echo "=========================================="
echo "并发任务数: $MAX_CONCURRENT_JOBS"
echo "物种列表: $TAXON_FILE"
echo "输出目录: $OUTPUT_DIR"
echo "主日志: $MAIN_LOG"
echo "错误日志: $ERROR_LOG"
echo "=========================================="
echo ""

# 计数器
total_count=0
success_count=0
failed_count=0
running_jobs=0

# 记录开始时间
START_TIME=$(date +%s)

# 读取物种文件并并发处理
while IFS= read -r taxon; do
    # 跳过空行和注释
    taxon=$(echo "$taxon" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    if [[ -z "$taxon" ]] || [[ "$taxon" =~ ^# ]]; then
        continue
    fi
    
    total_count=$((total_count + 1))
    
    # 等待直到运行任务数少于 MAX_CONCURRENT_JOBS
    wait_for_job "$MAX_CONCURRENT_JOBS"
    
    # 在后台处理物种
    (
        process_taxon "$taxon" "$OUTPUT_DIR" "$MAIN_LOG"
        exit $?
    ) &
    
    running_jobs=$((running_jobs + 1))
    
done < "$TAXON_FILE"

# 等待所有后台任务完成
echo ""
echo "等待所有后台任务完成..."
echo ""

while [ $(jobs -r | wc -l) -gt 0 ]; do
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] 运行中的任务数: $(jobs -r | wc -l)"
    sleep 5
done

# 计算完成时间
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

# 统计结果
if [ -f "$ERROR_LOG" ]; then
    failed_count=$(grep -c "✗" "$MAIN_LOG" 2>/dev/null || echo 0)
fi
success_count=$((total_count - failed_count))

# ===== 打印总结 =====
echo ""
echo "=========================================="
echo "完成统计"
echo "=========================================="
echo "总物种数: $total_count"
echo "成功: $success_count"
echo "失败: $failed_count"
echo "总耗时: $(printf '%d分%d秒\n' $((DURATION/60)) $((DURATION%60)))"
echo "=========================================="
echo ""

# 生成结果摘要
echo "生成结果摘要..."
SUMMARY_FILE="$OUTPUT_DIR/summary_$(date +%Y%m%d_%H%M%S).txt"
{
    echo "=========================================="
    echo "并发查询物种基因组信息 - 执行摘要"
    echo "=========================================="
    echo "执行日期: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "并发任务数: $MAX_CONCURRENT_JOBS"
    echo ""
    echo "统计信息:"
    echo "- 总物种数: $total_count"
    echo "- 成功: $success_count"
    echo "- 失败: $failed_count"
    echo "- 总耗时: $(printf '%d分%d秒\n' $((DURATION/60)) $((DURATION%60)))"
    echo ""
    echo "输出文件位置:"
    echo "- 主日志: $MAIN_LOG"
    echo "- 错误日志: $ERROR_LOG"
    echo "- 数据输出: $OUTPUT_DIR"
    echo "=========================================="
    
    if [ "$failed_count" -gt 0 ]; then
        echo ""
        echo "失败物种列表:"
        echo "------------------------------------------"
        grep "✗" "$MAIN_LOG" 2>/dev/null || echo "无失败记录"
        echo "------------------------------------------"
    fi
} > "$SUMMARY_FILE"

echo "摘要已保存到: $SUMMARY_FILE"
echo ""
echo "查看日志:"
echo "  主日志: tail -f $MAIN_LOG"
echo "  错误日志: tail -f $ERROR_LOG"
echo "  完整摘要: cat $SUMMARY_FILE"

# 如果有失败，返回非零值
if [ "$failed_count" -gt 0 ]; then
    echo ""
    echo "⚠ 部分物种查询失败，请查看错误日志"
    exit 1
else
    echo ""
    echo "✓ 所有物种查询成功"
    exit 0
fi
