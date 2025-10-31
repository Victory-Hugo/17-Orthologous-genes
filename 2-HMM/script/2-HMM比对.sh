#!/bin/bash

################################################################################
# 脚本名称: HMM蛋白序列比对工具
# 作者: BigLin
# 功能: 使用HMMER进行大规模并行HMM比对，支持断点续跑和进度显示
################################################################################

################################################################################
# 依赖软件和包:
# 1. hmmer (>=3.3) - 用于HMM比对
#    安装: apt-get install hmmer (Ubuntu/Debian)
#          conda install -c bioconda hmmer (Conda)
#
# 2. GNU parallel (>=20210222) - 用于并行处理
#    安装: apt-get install parallel (Ubuntu/Debian)
#          conda install -c conda-forge parallel (Conda)
#          brew install parallel (macOS)
#
# 3. pv (Pipe Viewer) - 用于显示进度条
#    安装: apt-get install pv (Ubuntu/Debian)
#          brew install pv (macOS)
#
# 4. GNU coreutils - 基础工具
#    安装: 通常已预装
################################################################################

# 启用严格错误检查
set -o pipefail

# ============================================================================
# 配置变量（全部写死）
# ============================================================================
LIST_FILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/2-HMM/conf/鲍曼faa.list.txt"
HMM_FILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/2-HMM/data/LolD.hmm"
OUTPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/2-HMM/output/"
CHECKPOINT_DIR="${OUTPUT_DIR}.checkpoint"
LOG_FILE="${OUTPUT_DIR}processing.log"
FAILED_LOG_FILE="${OUTPUT_DIR}failed_tasks.log"
STATE_FILE="${CHECKPOINT_DIR}/state.txt"
COMPLETED_FILE="${CHECKPOINT_DIR}/completed.txt"
TEMP_DIR="${CHECKPOINT_DIR}/temp"
JOBS="4"
CPU_PER_JOB="2"
E_VALUE="1e-50"

# ============================================================================
# 颜色定义
# ============================================================================
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
NC='\033[0m' # No Color

# ============================================================================
# 函数定义
# ============================================================================

# 彩色打印函数
print_info() {
    echo -e "${BLUE}[INFO]${NC} $1" | tee -a "$LOG_FILE"
}

print_success() {
    echo -e "${GREEN}[✓ SUCCESS]${NC} $1" | tee -a "$LOG_FILE"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1" | tee -a "$LOG_FILE"
}

print_error() {
    echo -e "${RED}[✗ ERROR]${NC} $1" | tee -a "$LOG_FILE"
}

print_header() {
    echo ""
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}$1${NC}"
    echo -e "${CYAN}========================================${NC}"
    echo ""
}

# 初始化目录和检查点
initialize_dirs() {
    print_info "初始化目录结构..."
    mkdir -p "$OUTPUT_DIR" || { print_error "无法创建输出目录: $OUTPUT_DIR"; exit 1; }
    mkdir -p "$CHECKPOINT_DIR" || { print_error "无法创建检查点目录: $CHECKPOINT_DIR"; exit 1; }
    mkdir -p "$TEMP_DIR" || { print_error "无法创建临时目录: $TEMP_DIR"; exit 1; }
    
    # 初始化日志文件
    : > "$LOG_FILE"
    : > "$FAILED_LOG_FILE"
}

# 检查依赖工具
check_dependencies() {
    print_info "检查依赖软件..."
    local missing=0
    
    if ! command -v hmmsearch &> /dev/null; then
        print_error "未找到 hmmsearch，请安装 HMMER"
        missing=1
    fi
    
    if ! command -v parallel &> /dev/null; then
        print_error "未找到 parallel，请安装 GNU parallel"
        missing=1
    fi
    
    if ! command -v pv &> /dev/null; then
        print_warning "未找到 pv，进度条将不可用"
    fi
    
    if [[ $missing -eq 1 ]]; then
        exit 1
    fi
    
    print_success "所有依赖软件检查通过"
}

# 检查输入文件
check_input_files() {
    print_info "检查输入文件..."
    
    if [[ ! -f "$LIST_FILE" ]]; then
        print_error "列表文件不存在: $LIST_FILE"
        exit 1
    fi
    
    if [[ ! -f "$HMM_FILE" ]]; then
        print_error "HMM文件不存在: $HMM_FILE"
        exit 1
    fi
    
    local empty_lines=0
    local total_lines=0
    while IFS= read -r line; do
        ((total_lines++))
        if [[ -z "$line" ]] || [[ "$line" == " "* ]]; then
            ((empty_lines++))
        elif [[ ! -f "$line" ]]; then
            print_warning "FAA文件不存在: $line"
        fi
    done < "$LIST_FILE"
    
    print_info "列表文件总行数: $total_lines (空行: $empty_lines)"
    print_success "输入文件检查完成"
}

# 加载已完成的任务
load_completed_tasks() {
    if [[ -f "$COMPLETED_FILE" ]]; then
        mapfile -t COMPLETED_TASKS < "$COMPLETED_FILE"
    else
        COMPLETED_TASKS=()
    fi
}

# 检查任务是否已完成
is_task_completed() {
    local filename="$1"
    for completed in "${COMPLETED_TASKS[@]}"; do
        if [[ "$completed" == "$filename" ]]; then
            return 0
        fi
    done
    return 1
}

# 记录已完成的任务
mark_task_completed() {
    local filename="$1"
    echo "$filename" >> "$COMPLETED_FILE"
}

# 处理单个FAA文件的函数
process_faa_file() {
    local faa_file="$1"
    local output_dir="$2"
    local hmm_file="$3"
    local checkpoint_dir="$4"
    local temp_dir="$5"
    local cpu_per_job="$6"
    local e_value="$7"
    
    # 跳过空行和注释行
    [[ -z "$faa_file" ]] || [[ "$faa_file" == " "* ]] && return 0
    
    # 检查文件是否存在
    if [[ ! -f "$faa_file" ]]; then
        echo "[ERROR] 文件不存在: $faa_file" >> "${temp_dir}/errors.log"
        return 1
    fi
    
    # 提取文件名（不含路径和扩展名）
    local filename=$(basename "$faa_file" .faa)
    
    # 检查是否已完成
    if [[ -f "${checkpoint_dir}/completed.txt" ]] && grep -q "^${filename}$" "${checkpoint_dir}/completed.txt"; then
        echo "[SKIP] $filename - 已完成" >> "${temp_dir}/process.log"
        return 0
    fi
    
    # 设置输出文件名
    local output_file="${output_dir}/hits_${filename}.tbl"
    local temp_output="${temp_dir}/${filename}.tbl.tmp"
    
    # 确保临时目录存在
    mkdir -p "$temp_dir" 2>/dev/null
    
    # 执行 hmmsearch，捕获完整错误信息
    local error_log="${temp_dir}/${filename}.err"
    if hmmsearch --cpu "$cpu_per_job" \
        -E "$e_value" \
        --tblout "$temp_output" \
        "$hmm_file" \
        "$faa_file" 2> "$error_log"; then
        
        # 移动临时文件到最终位置
        if [[ -f "$temp_output" ]]; then
            mv "$temp_output" "$output_file"
            echo "[SUCCESS] $filename" >> "${temp_dir}/success.log"
            rm -f "$error_log"
            return 0
        else
            echo "[FAILED] $filename - 输出文件未生成" >> "${temp_dir}/failed.log"
            cat "$error_log" >> "${temp_dir}/failed.log"
            return 1
        fi
    else
        echo "[FAILED] $filename - hmmsearch返回错误" >> "${temp_dir}/failed.log"
        cat "$error_log" >> "${temp_dir}/failed.log"
        rm -f "$temp_output" "$error_log"
        return 1
    fi
}

export -f process_faa_file

# 清理临时文件和未完成的任务
cleanup_and_exit() {
    local signal="$1"
    print_warning "收到中止信号 ($signal)，进行清理..."
    
    # 等待所有后台任务完成
    jobs -p | xargs -r kill 2>/dev/null
    wait
    
    # 清理临时文件
    print_info "清理临时文件..."
    rm -f "$TEMP_DIR"/*.tbl.tmp
    
    # 合并日志
    if [[ -f "$TEMP_DIR/success.log" ]]; then
        cat "$TEMP_DIR/success.log" >> "$OUTPUT_DIR/success.log"
    fi
    
    if [[ -f "$TEMP_DIR/failed.log" ]]; then
        cat "$TEMP_DIR/failed.log" >> "$FAILED_LOG_FILE"
    fi
    
    print_warning "脚本已中止，已完成的任务已保留"
    exit 130
}

# 处理进度显示
show_progress() {
    local total="$1"
    local completed="$2"
    local percent=$((completed * 100 / total))
    local bar_length=40
    local filled=$((percent * bar_length / 100))
    local empty=$((bar_length - filled))
    
    printf "\r${CYAN}进度: [${GREEN}"
    printf "%${filled}s" | tr ' ' '='
    printf "${NC}%${empty}s${CYAN}] %3d%% (%d/%d)${NC}" | tr ' ' '-' | sed "s/-//g" >&1
    printf "%${empty}s${CYAN}] %3d%% (%d/%d)${NC}" | sed "s/ //g"
    printf "\r${CYAN}进度: [${GREEN}"
    for ((i=0; i<filled; i++)); do printf "="; done
    for ((i=0; i<empty; i++)); do printf "-"; done
    printf "${NC}] ${percent}%% (${completed}/${total})${NC}"
}

# ============================================================================
# 主程序
# ============================================================================

# 设置中止信号处理
trap 'cleanup_and_exit SIGINT' SIGINT SIGTERM

# 打印头部信息
print_header "HMM蛋白序列比对工具 v1.0 - 作者: BigLin"

# 初始化
initialize_dirs
check_dependencies
check_input_files
load_completed_tasks

print_info "开始处理..."
print_info "输出目录: $OUTPUT_DIR"
print_info "并行任务数: $JOBS"
print_info "每任务CPU线程数: $CPU_PER_JOB"

# 统计总任务数
local_total_tasks=$(grep -cv '^\s*$' "$LIST_FILE")
print_info "总任务数: $local_total_tasks"

# 使用 parallel 进行并行处理，带进度条
local_completed_before=$((${#COMPLETED_TASKS[@]}))
print_info "已完成任务数: $local_completed_before"

# 初始化临时日志文件
: > "$TEMP_DIR/process.log"
: > "$TEMP_DIR/success.log"
: > "$TEMP_DIR/failed.log"
: > "$TEMP_DIR/errors.log"

cat "$LIST_FILE" | \
    grep -v '^\s*$' | \
    parallel --jobs "$JOBS" \
        --tag \
        --line-buffer \
        --joblog "${CHECKPOINT_DIR}/joblog.txt" \
        process_faa_file {} "$OUTPUT_DIR" "$HMM_FILE" "$CHECKPOINT_DIR" "$TEMP_DIR" "$CPU_PER_JOB" "$E_VALUE"

# 等待所有后台任务完成
wait

# 收集结果
print_info "收集处理结果..."

local_success_count=0
local_failed_count=0

if [[ -f "$TEMP_DIR/success.log" ]]; then
    local_success_count=$(wc -l < "$TEMP_DIR/success.log")
    while IFS= read -r line; do
        # 格式: [SUCCESS] filename
        filename=$(echo "$line" | sed 's/\[SUCCESS\] //')
        [[ -n "$filename" ]] && mark_task_completed "$filename"
    done < "$TEMP_DIR/success.log"
    print_success "成功完成任务数: $local_success_count"
fi

if [[ -f "$TEMP_DIR/failed.log" ]]; then
    local_failed_count=$(wc -l < "$TEMP_DIR/failed.log")
    cat "$TEMP_DIR/failed.log" >> "$FAILED_LOG_FILE"
    if [[ $local_failed_count -gt 0 ]]; then
        print_error "失败任务数: $local_failed_count"
        print_warning "失败详情请查看: $FAILED_LOG_FILE"
    fi
fi

# 打印完成信息
print_header "处理完成"
print_success "所有任务处理完成！"
print_info "输出目录: $OUTPUT_DIR"
print_info "日志文件: $LOG_FILE"
[[ -s "$FAILED_LOG_FILE" ]] && print_warning "失败日志: $FAILED_LOG_FILE"

echo ""
print_success "脚本执行成功！" 