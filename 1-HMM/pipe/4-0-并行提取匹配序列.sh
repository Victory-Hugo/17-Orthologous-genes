#!/bin/bash

################################################################################
# 脚本名称: 并行提取最佳匹配序列
# 作者: BigLin
# 描述: 使用 GNU Parallel 并行处理 HMM 搜索结果，提取最佳匹配序列
#
# 依赖软件和包:
#   - GNU Parallel (brew install parallel 或 apt-get install parallel)
#   - Python 3.6+
#   - 依赖的 Python 模块: (在下方指定)
#
# 使用方法:
#   ./4-并行提取最佳匹配序列.sh
#
# 中断处理:
#   按 Ctrl+C 中断脚本，已完成的任务会被保留，未完成的任务会被删除
#
################################################################################

set -o pipefail

# ==================== 颜色定义 ====================
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

#* 输出目录结构示例:
#*├── 基因名字1
#*├── 基因名字2
#*├── 基因名字3

# ==================== 全部变量写死 ====================
BASEDIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM"
SCRIPT_DIR="${BASEDIR}/script"
LIST_FILE="${BASEDIR}/conf/1-prodigal_faa_list.txt"

TBL_FILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/tbl_merge.tsv" #! 这里写上一步合成的tbl文件
OUTPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/" #! 这里写提取序列的输出目录

JOBS=4
STATE_DIR="${OUTPUT_DIR}/.job_state"
LOG_FILE="${OUTPUT_DIR}/parallel_extract.log"
TEMP_FILE_LIST="/tmp/extract_file_list_$$.txt"
CHECKPOINT_FILE="${STATE_DIR}/checkpoint_$$.json"

# ==================== 命中过滤参数（可以按需修改） ====================
MAX_EVALUE="1e-20"        # 允许的最大 E-value，留空表示不过滤
MIN_COVERAGE="0.6"        # 命中序列长度 / HMM 模型长度，留空表示不过滤
MIN_LENGTH="180"          # 命中序列的最小氨基酸长度，留空表示不过滤
EXTRACT_MODE="only_one"   # all: 输出全部命中；only_one: 每株只保留最高分
HMM_MODEL_LENGTH="225"    # LolD.hmm 的 LENG 数值（用于覆盖度计算）

# ==================== 函数定义 ====================

# 彩色打印函数
print_color() {
    local color=$1
    shift
    echo -e "${color}$*${NC}"
}

# 打印信息
print_info() {
    print_color "$BLUE" "ℹ️  $*"
}

# 打印成功
print_success() {
    print_color "$GREEN" "✓ $*"
}

# 打印警告
print_warning() {
    print_color "$YELLOW" "⚠️  $*"
}

# 打印错误
print_error() {
    print_color "$RED" "✗ $*"
}

# 检查依赖
check_dependencies() {
    print_info "检查依赖软件..."
    
    local missing=0
    
    if ! command -v parallel &> /dev/null; then
        print_error "GNU Parallel 未安装"
        print_info "安装方法: sudo apt-get install parallel (Linux) 或 brew install parallel (macOS)"
        missing=1
    else
        print_success "GNU Parallel 已安装: $(parallel --version | head -n1)"
    fi
    
    if ! command -v python3 &> /dev/null; then
        print_error "Python3 未安装"
        missing=1
    else
        print_success "Python3 已安装: $(python3 --version)"
    fi
    
    if [ $missing -eq 1 ]; then
        print_error "缺少必要的依赖软件，请先安装"
        exit 1
    fi
    
    echo ""
}

# 初始化环境
init_env() {
    print_info "初始化环境..."
    
    # 创建必要的目录
    mkdir -p "$OUTPUT_DIR" "$STATE_DIR" "$(dirname "$LOG_FILE")"
    
    # 清空日志文件
    > "$LOG_FILE"
    
    print_success "环境初始化完成"
    echo ""
}

# 验证输入文件
validate_inputs() {
    print_info "验证输入文件..."
    
    if [ ! -f "$LIST_FILE" ]; then
        print_error "列表文件不存在: $LIST_FILE"
        exit 1
    fi
    
    
    if [ ! -d "$SCRIPT_DIR" ]; then
        print_error "脚本目录不存在: $SCRIPT_DIR"
        exit 1
    fi
    
    if [ ! -f "$SCRIPT_DIR/4-提取匹配序列.py" ]; then
        print_error "Python 脚本不存在: $SCRIPT_DIR/4-提取匹配序列.py"
        exit 1
    fi
    
    print_success "所有输入文件验证通过"
    echo ""
}

# 获取需要处理的文件列表
get_file_list() {
    cat "$LIST_FILE" > "$TEMP_FILE_LIST"
}

# 获取待处理文件数
get_total_count() {
    wc -l < "$TEMP_FILE_LIST"
}

build_python_command() {
    local cmd="python3 '$SCRIPT_DIR/4-提取匹配序列.py' {} '$TBL_FILE' '$OUTPUT_DIR'"
    if [[ -n "$MAX_EVALUE" ]]; then
        cmd+=" --max-evalue '$MAX_EVALUE'"
    fi
    if [[ -n "$MIN_COVERAGE" ]]; then
        cmd+=" --min-coverage '$MIN_COVERAGE'"
    fi
    if [[ -n "$MIN_LENGTH" ]]; then
        cmd+=" --min-length '$MIN_LENGTH'"
    fi
    if [[ -n "$EXTRACT_MODE" ]]; then
        cmd+=" --mode '$EXTRACT_MODE'"
    fi
    if [[ -n "$HMM_MODEL_LENGTH" ]]; then
        cmd+=" --model-length '$HMM_MODEL_LENGTH'"
    fi
    echo "$cmd"
}


# 中断处理函数
cleanup_on_interrupt() {
    print_error "\n检测到用户中断 (Ctrl+C)..."
    print_warning "正在清理未完成的任务..."
    
    # 等待 parallel 进程自行清理
    wait 2>/dev/null || true
    
    print_info "删除临时文件..."
    rm -f "$TEMP_FILE_LIST" 2>/dev/null || true
    
    print_warning "脚本已中断，已完成的任务已保留"
    print_info "日志文件: $LOG_FILE"
    exit 130  # 标准中断退出码
}

# 捕获 SIGINT 和 SIGTERM 信号
trap cleanup_on_interrupt SIGINT SIGTERM

# 运行并行处理
run_parallel_processing() {
    print_info "开始并行处理..."
    echo ""
    
    local total=$(get_total_count)
    local python_cmd
    python_cmd=$(build_python_command)
    print_info "待处理文件总数: $total"
    print_info "并行任务数: $JOBS"
    print_info "过滤条件: E-value<=${MAX_EVALUE:-不限}, 覆盖度>=${MIN_COVERAGE:-不限}, 长度>=${MIN_LENGTH:-不限}, 模式=${EXTRACT_MODE:-all}"
    echo ""
    
    # 使用 parallel 并行处理，添加进度条
    cat "$TEMP_FILE_LIST" | \
    parallel \
        --jobs "$JOBS" \
        --tag \
        --bar \
        --line-buffer \
        --joblog "$STATE_DIR/job_log_$$.txt" \
        "$python_cmd 2>&1" >> "$LOG_FILE" 2>&1
    local exit_code=$?
    
    if [ $exit_code -eq 0 ]; then
        print_success "所有任务处理完成"
    elif [ $exit_code -eq 130 ]; then
        print_error "任务被中断"
        return 130
    else
        print_error "部分任务处理失败 (退出码: $exit_code)"
        print_info "请查看日志: $LOG_FILE"
        return $exit_code
    fi
    
    echo ""
}

# 显示处理结果统计
show_statistics() {
    print_info "处理结果统计..."
    echo ""
    
    if [ -f "$STATE_DIR/job_log_$$.txt" ]; then
        local total=$(tail -n +2 "$STATE_DIR/job_log_$$.txt" 2>/dev/null | wc -l)
        local succeeded=$(awk 'NR>1 && $7=="0" {count++} END {print count+0}' "$STATE_DIR/job_log_$$.txt" 2>/dev/null)
        local failed=$((total - succeeded))
        
        print_info "总任务数: $total"
        print_success "成功: $succeeded"
        if [ "${failed:-0}" -gt 0 ]; then
            print_error "失败: $failed"
        fi
    else
        print_warning "未找到任务日志文件"
    fi
    
    print_info "输出目录: $OUTPUT_DIR"
    print_info "日志文件: $LOG_FILE"
    echo ""
}

# 清理临时文件
cleanup() {
    rm -f "$TEMP_FILE_LIST" 2>/dev/null || true
}

# ==================== 主程序执行 ====================

main() {
    print_color "$CYAN" "
╔═══════════════════════════════════════════════════════════╗
║     并行提取最佳匹配序列 - 作者: BigLin              ║
╚═══════════════════════════════════════════════════════════╝
    "
    
    # 检查依赖
    check_dependencies
    
    # 初始化环境
    init_env
    
    # 验证输入
    validate_inputs
    
    # 获取文件列表
    get_file_list
    
    # 运行并行处理
    run_parallel_processing
    local exit_code=$?
    
    if [ $exit_code -ne 130 ]; then
        # 显示统计信息
        show_statistics
    fi
    
    # 清理临时文件
    cleanup
    
    if [ $exit_code -eq 0 ]; then
        print_color "$GREEN" "
╔═══════════════════════════════════════════════════════════╗
║                     ✓ 脚本执行成功                    ║
╚═══════════════════════════════════════════════════════════╝
        "
    else
        print_color "$RED" "
╔═══════════════════════════════════════════════════════════╗
║                     ✗ 脚本执行失败                    ║
╚═══════════════════════════════════════════════════════════╝
        "
    fi
    
    exit $exit_code
}

# 执行主程序
main "$@"
python3 /home/luolintao/test_mail.py "2-HMM/script/4-并行提取最佳匹配序列.sh任务完成通知" "<p>2-HMM/script/4-并行提取最佳匹配序列.sh分析已完成，请查看结果目录。</p>"
