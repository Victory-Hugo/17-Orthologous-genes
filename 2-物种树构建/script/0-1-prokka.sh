#!/bin/bash

###############################################################################
#                        Prokka 并行处理脚本 (优化版)
# 功能: 并行运行 Prokka，支持断点续跑、彩色输出、优雅中
###############################################################################

# ==================== 变量写死配置 ====================
INPUT_LIST="/mnt/f/15_Bam_Tam/2-物种树/prokka/conf/correct_fna_input.txt"
OUTPUT_DIR="/mnt/f/15_Bam_Tam/2-物种树/prokka/output"
LOG_DIR="/mnt/f/15_Bam_Tam/2-物种树/prokka/logs"
STATUS_DIR="/mnt/f/15_Bam_Tam/2-物种树/prokka/status"
PARALLEL_JOBS=4
THREADS_PER_JOB=8
TEMP_FIFO="/tmp/prokka_pipe_$$"

# ==================== 颜色定义 ====================
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
NC='\033[0m' # No Color

# ==================== 函数定义 ====================

# 彩色打印函数
print_info() {
    echo -e "${BLUE}[INFO]${NC} $*"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $*"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $*"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $*" >&2
}

print_start() {
    echo -e "${CYAN}[START]${NC} $*"
}

print_done() {
    echo -e "${MAGENTA}[DONE]${NC} $*"
}

# 初始化函数
init() {
    print_info "初始化脚本环境..."
    
    # 创建必要目录
    mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$STATUS_DIR"
    
    # 检查 parallel 是否安装
    if ! command -v parallel &> /dev/null; then
        print_error "GNU Parallel 未安装。请运行: apt-get install parallel (或 brew install parallel)"
        exit 1
    fi
    
    # 检查 prokka 是否可用
    if ! command -v prokka &> /dev/null; then
        print_error "Prokka 不可用。请激活包含 Prokka 的环境。"
        exit 1
    fi
    
    # 检查输入文件
    if [ ! -f "$INPUT_LIST" ]; then
        print_error "输入文件不存在: $INPUT_LIST"
        exit 1
    fi
    
    print_success "初始化完成"
}

# 检查任务是否已完成
is_task_completed() {
    local sample="$1"
    [ -f "$STATUS_DIR/${sample}.completed" ]
}

# 标记任务完成
mark_task_completed() {
    local sample="$1"
    touch "$STATUS_DIR/${sample}.completed"
}

# 标记任务失败
mark_task_failed() {
    local sample="$1"
    touch "$STATUS_DIR/${sample}.failed"
}

# 清理失败标记
clear_task_failed() {
    local sample="$1"
    rm -f "$STATUS_DIR/${sample}.failed"
}

# 获取未完成的任务列表
get_pending_tasks() {
    local task_count=0
    while IFS= read -r fna; do
        [ -z "$fna" ] && continue
        local base=$(basename "$fna")
        # 处理 _CDS.fna 或其他后缀
        local sample=${base%.*}          # 移除 .fna
        sample=${sample%_CDS}            # 移除 _CDS 后缀（如果有）
        
        if ! is_task_completed "$sample"; then
            echo "$fna"
            ((task_count++))
        fi
    done < "$INPUT_LIST"
    
    return $task_count
}

# Prokka 处理函数
run_prokka() {
    local fna="$1"
    local base=$(basename "$fna")
    # 处理 _CDS.fna 或其他后缀
    local sample=${base%.*}              # 移除 .fna
    sample=${sample%_CDS}                # 移除 _CDS 后缀（如果有）
    local output_sample_dir="$OUTPUT_DIR/$sample"
    local log_file="$LOG_DIR/${sample}.log"
    
    # 检查是否已完成
    if is_task_completed "$sample"; then
        print_info "[$sample] 已完成，跳过"
        return 0
    fi
    
    # 清理之前的失败标记
    clear_task_failed "$sample"
    
    print_start "[$sample] 开始处理: $fna"
    
    # 创建输出目录
    mkdir -p "$output_sample_dir"
    
    # 运行 Prokka
    if prokka "$fna" \
        --outdir "$output_sample_dir" \
        --prefix "$sample" \
        --cpus "$THREADS_PER_JOB" \
        --force \
        --kingdom Bacteria \
        --genus Unknown \
        --usegenus \
        --addgenes \
        --locustag "$sample" \
        --compliant \
        > "$log_file" 2>&1; then
        
        mark_task_completed "$sample"
        print_done "[$sample] 处理完成"
        return 0
    else
        mark_task_failed "$sample"
        print_error "[$sample] 处理失败，详见: $log_file"
        return 1
    fi
}

export -f run_prokka is_task_completed mark_task_completed mark_task_failed clear_task_failed
export -f print_info print_success print_warning print_error print_start print_done
export OUTPUT_DIR LOG_DIR STATUS_DIR THREADS_PER_JOB
export RED GREEN YELLOW BLUE CYAN MAGENTA NC

# 中断处理函数
trap_exit() {
    print_warning "检测到用户中断信号..."
    sleep 1
    print_warning "等待进行中的任务完成..."
    sleep 2
    print_warning "清理未完成的临时文件..."
    
    # 列出未完成的任务
    local incomplete_count=0
    while IFS= read -r fna; do
        [ -z "$fna" ] && continue
        local base=$(basename "$fna")
        local sample=${base%.*}
        if [ -f "$STATUS_DIR/${sample}.failed" ]; then
            local incomplete_dir="$OUTPUT_DIR/$sample"
            if [ -d "$incomplete_dir" ]; then
                print_warning "删除未完成的: $incomplete_dir"
                rm -rf "$incomplete_dir"
                rm -f "$STATUS_DIR/${sample}.failed"
                ((incomplete_count++))
            fi
        fi
    done < "$INPUT_LIST"
    
    print_warning "已删除 $incomplete_count 个未完成的任务"
    print_info "脚本中断，已保留所有完成的任务"
    exit 0
}

# 设置信号捕获
trap trap_exit SIGINT SIGTERM

# ==================== 主程序 ====================

print_info "=========================================="
print_info "Prokka 并行处理脚本 (优化版)"
print_info "=========================================="

# 初始化
init

# 统计任务
total_tasks=0
while IFS= read -r fna; do
    [ -z "$fna" ] && continue
    ((total_tasks++))
done < "$INPUT_LIST"

# 检查待处理任务
pending_count=0
while IFS= read -r fna; do
    [ -z "$fna" ] && continue
    base=$(basename "$fna")
    sample=${base%.*}
    if ! is_task_completed "$sample"; then
        ((pending_count++))
    fi
done < "$INPUT_LIST"

completed_count=$((total_tasks - pending_count))

print_info "任务统计: 总计=$total_tasks, 已完成=$completed_count, 待处理=$pending_count"

if [ "$pending_count" -eq 0 ]; then
    print_success "所有任务已完成！"
    exit 0
fi

# 运行并行处理
print_info "开始并行处理 (并行数: $PARALLEL_JOBS, 每任务线程数: $THREADS_PER_JOB)..."

get_pending_tasks | parallel \
    --jobs "$PARALLEL_JOBS" \
    --line-buffer \
    --tag \
    "run_prokka {}"

# 检查最终状态
print_info "处理完成，检查最终状态..."

failed_count=0
success_count=0

while IFS= read -r fna; do
    [ -z "$fna" ] && continue
    base=$(basename "$fna")
    sample=${base%.*}
    
    if is_task_completed "$sample"; then
        ((success_count++))
    else
        ((failed_count++))
        print_error "失败: $sample"
    fi
done < "$INPUT_LIST"

print_info "=========================================="
print_success "成功: $success_count 个任务"
if [ "$failed_count" -gt 0 ]; then
    print_error "失败: $failed_count 个任务"
fi
print_info "日志位置: $LOG_DIR"
print_info "状态追踪: $STATUS_DIR"
print_info "=========================================="

if [ "$failed_count" -eq 0 ]; then
    print_success "所有任务处理成功！"
    exit 0
else
    print_warning "部分任务处理失败，请查看日志"
    exit 1
fi
