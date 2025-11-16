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
LIST_FILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/conf/1-prodigal_faa_list.txt"
# TXT文件，每行一个HMM文件的绝对路径
HMM_FILE_TXT="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/conf/2-hmm库_list.txt"
OUTPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/"
CHECKPOINT_DIR="${OUTPUT_DIR}.checkpoint"
LOG_FILE="${OUTPUT_DIR}processing.log"
FAILED_LOG_FILE="${OUTPUT_DIR}failed_tasks.log"
STATE_FILE="${CHECKPOINT_DIR}/state.txt"
COMPLETED_FILE="${CHECKPOINT_DIR}/completed.txt"
TEMP_DIR="${CHECKPOINT_DIR}/temp"
JOBS="8"
CPU_PER_JOB="2"
E_VALUE="1e-1"

declare -a FAA_FILES=()
declare -a HMM_FILES=()
TOTAL_FAA_FILES=0
TOTAL_HMM_FILES=0
TOTAL_TASKS=0

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
    
    if [[ ! -f "$HMM_FILE_TXT" ]]; then
        print_error "HMM列表文件不存在: $HMM_FILE_TXT"
        exit 1
    fi
    
    FAA_FILES=()
    HMM_FILES=()

    local faa_empty_lines=0
    local faa_total_entries=0
    while IFS= read -r line || [[ -n "$line" ]]; do
        line="${line%$'\r'}"
        if [[ "$line" =~ ^[[:space:]]*$ ]]; then
            ((faa_empty_lines++))
            continue
        fi

        local trimmed
    trimmed=$(printf '%s\n' "$line" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
        [[ "$trimmed" == \#* ]] && continue

        ((faa_total_entries++))
        if [[ ! -f "$trimmed" ]]; then
            print_warning "FAA文件不存在: $trimmed"
            continue
        fi
        FAA_FILES+=("$trimmed")
    done < "$LIST_FILE"

    TOTAL_FAA_FILES=${#FAA_FILES[@]}
    print_info "FAA列表有效文件数: $TOTAL_FAA_FILES (总行数: $faa_total_entries, 空行: $faa_empty_lines)"

    local hmm_empty_lines=0
    local hmm_total_entries=0
    while IFS= read -r line || [[ -n "$line" ]]; do
        line="${line%$'\r'}"
        if [[ "$line" =~ ^[[:space:]]*$ ]]; then
            ((hmm_empty_lines++))
            continue
        fi

    local trimmed
    trimmed=$(printf '%s\n' "$line" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
        [[ "$trimmed" == \#* ]] && continue

        ((hmm_total_entries++))
        if [[ ! -f "$trimmed" ]]; then
            print_warning "HMM文件不存在: $trimmed"
            continue
        fi
        HMM_FILES+=("$trimmed")
    done < "$HMM_FILE_TXT"

    TOTAL_HMM_FILES=${#HMM_FILES[@]}
    print_info "HMM列表有效文件数: $TOTAL_HMM_FILES (总行数: $hmm_total_entries, 空行: $hmm_empty_lines)"

    TOTAL_TASKS=$((TOTAL_FAA_FILES * TOTAL_HMM_FILES))
    if [[ $TOTAL_TASKS -eq 0 ]]; then
        print_error "没有可执行的任务，请检查FAA列表和HMM列表"
        exit 1
    fi

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
    local task_key="$1"
    for completed in "${COMPLETED_TASKS[@]}"; do
        if [[ "$completed" == "$task_key" ]]; then
            return 0
        fi
    done
    return 1
}

# 记录已完成的任务
mark_task_completed() {
    local task_key="$1"
    if [[ -z "$task_key" ]]; then
        return
    fi

    if [[ -f "$COMPLETED_FILE" ]] && grep -qxF "$task_key" "$COMPLETED_FILE"; then
        return
    fi

    echo "$task_key" >> "$COMPLETED_FILE"
    COMPLETED_TASKS+=("$task_key")
}

# 处理 HMM 与 FAA 组合的函数
process_hmm_faa_combo() {
    local hmm_file="$1"
    local faa_file="$2"
    local output_dir="$3"
    local checkpoint_dir="$4"
    local temp_dir="$5"
    local cpu_per_job="$6"
    local e_value="$7"

    [[ -z "$hmm_file" ]] && return 0
    [[ -z "$faa_file" ]] && return 0

    if [[ ! -f "$hmm_file" ]]; then
        echo "[ERROR] HMM不存在: $hmm_file" >> "${temp_dir}/errors.log"
        return 1
    fi

    if [[ ! -f "$faa_file" ]]; then
        echo "[ERROR] FAA不存在: $faa_file" >> "${temp_dir}/errors.log"
        return 1
    fi

    local hmm_filename=$(basename "$hmm_file")
    local hmm_basename="${hmm_filename%.*}"
    local faa_filename=$(basename "$faa_file")
    local faa_basename="${faa_filename%.*}"
    local task_key="${hmm_basename}|${faa_basename}"

    if [[ -f "${checkpoint_dir}/completed.txt" ]] && grep -qxF "$task_key" "${checkpoint_dir}/completed.txt"; then
        echo "[SKIP] ${task_key} - 已完成" >> "${temp_dir}/process.log"
        return 0
    fi

    local output_subdir="${output_dir}/${hmm_basename}"
    local output_file="${output_subdir}/${faa_basename}.tbl"
    local temp_output="${temp_dir}/${hmm_basename}_${faa_basename}.tbl.tmp"

    mkdir -p "$output_subdir" "$temp_dir" 2>/dev/null

    local error_log="${temp_dir}/${hmm_basename}_${faa_basename}.err"
    if hmmsearch --cpu "$cpu_per_job" \
        -E "$e_value" \
        --tblout "$temp_output" \
        "$hmm_file" \
        "$faa_file" 2> "$error_log"; then

        if [[ -f "$temp_output" ]]; then
            mv "$temp_output" "$output_file"
            echo "[SUCCESS] $task_key" >> "${temp_dir}/success.log"
            rm -f "$error_log"
            return 0
        else
            echo "[FAILED] $task_key - 输出文件未生成" >> "${temp_dir}/failed.log"
            cat "$error_log" >> "${temp_dir}/failed.log"
            return 1
        fi
    else
        echo "[FAILED] $task_key - hmmsearch返回错误" >> "${temp_dir}/failed.log"
        cat "$error_log" >> "${temp_dir}/failed.log"
        rm -f "$temp_output" "$error_log"
        return 1
    fi
}

export -f process_hmm_faa_combo
export -f is_task_completed

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
print_header "HMM序列比对工具 v1.0 - 作者: BigLin"

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
local_total_tasks=$TOTAL_TASKS
print_info "总任务数: $local_total_tasks"

# 使用 parallel 进行并行处理，带进度条
local_completed_before=$((${#COMPLETED_TASKS[@]}))
print_info "已完成任务数: $local_completed_before"

# 初始化临时日志文件
: > "$TEMP_DIR/process.log"
: > "$TEMP_DIR/success.log"
: > "$TEMP_DIR/failed.log"
: > "$TEMP_DIR/errors.log"

tasks_file="${TEMP_DIR}/tasks.tsv"
: > "$tasks_file"
for hmm_file in "${HMM_FILES[@]}"; do
    for faa_file in "${FAA_FILES[@]}"; do
        printf '%s\t%s\n' "$hmm_file" "$faa_file" >> "$tasks_file"
    done
done

if [[ ! -s "$tasks_file" ]]; then
    print_warning "任务列表为空，未执行任何操作"
else
    parallel --jobs "$JOBS" \
        --colsep '\t' \
        --tag \
        --line-buffer \
        --joblog "${CHECKPOINT_DIR}/joblog.txt" \
        process_hmm_faa_combo {1} {2} "$OUTPUT_DIR" "$CHECKPOINT_DIR" "$TEMP_DIR" "$CPU_PER_JOB" "$E_VALUE" \
        :::: "$tasks_file"
fi

# 打印完成信息
print_header "处理完成"
print_success "所有任务处理完成！"
print_info "输出目录: $OUTPUT_DIR"
print_info "日志文件: $LOG_FILE"
[[ -s "$FAILED_LOG_FILE" ]] && print_warning "失败日志: $FAILED_LOG_FILE"

echo ""
print_success "脚本执行成功！" 
python3 /home/luolintao/test_mail.py "2-HMM/script/2-HMM比对.sh任务完成通知" "<p>2-HMM/script/2-HMM比对.sh分析已完成，请查看结果目录。</p>"