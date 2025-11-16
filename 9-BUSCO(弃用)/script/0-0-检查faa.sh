#!/bin/bash

################################################################################
# 检查样本中是否存在 .faa 文件
# 功能：
#   1. 并行检查所有样本目录
#   2. 支持断点续跑
#   3. 支持彩色输出
#   4. 支持优雅停止（Ctrl+C）
#   5. 生成详细报告
################################################################################

################################################################################
# 配置部分 - 所有变量写死
################################################################################
INPUT_DIR="/mnt/f/15_Bam_Tam/5-补齐更多物种/download"
OUTPUT_DIR="/mnt/f/15_Bam_Tam/5-补齐更多物种/output"
REPORT_FILE="${OUTPUT_DIR}/faa_check_report.txt"
JOBS=10

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# 进度追踪相关
PROGRESS_DIR="${OUTPUT_DIR}/.progress_faa_check"
LOCK_FILE="${PROGRESS_DIR}/lock"
CHECKED_FILE="${PROGRESS_DIR}/checked.txt"
HAS_FAA_FILE="${PROGRESS_DIR}/has_faa.txt"
NO_FAA_FILE="${PROGRESS_DIR}/no_faa.txt"

################################################################################
# 颜色打印函数
################################################################################
log_info() {
    echo -e "${BLUE}[INFO]${NC} $*"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $*"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $*"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $* " >&2
}

log_running() {
    echo -e "${CYAN}[RUNNING]${NC} $*"
}

################################################################################
# 初始化函数
################################################################################
init_progress_tracking() {
    mkdir -p "${PROGRESS_DIR}"
    mkdir -p "${OUTPUT_DIR}"
    
    # 如果进度文件不存在，则创建
    if [ ! -f "${CHECKED_FILE}" ]; then
        touch "${CHECKED_FILE}"
    fi
    
    if [ ! -f "${HAS_FAA_FILE}" ]; then
        touch "${HAS_FAA_FILE}"
    fi
    
    if [ ! -f "${NO_FAA_FILE}" ]; then
        touch "${NO_FAA_FILE}"
    fi
}

################################################################################
# 清理资源函数
################################################################################
cleanup() {
    log_warning "捕获到中断信号，开始清理..."
    
    # 等待所有后台进程完成
    wait 2>/dev/null || true
    
    # 清理锁文件
    rm -f "${LOCK_FILE}"
    
    log_warning "清理完成。已完成的检查已保留在进度文件中"
}

# 设置信号处理
trap cleanup SIGINT SIGTERM

################################################################################
# 主要检查函数
################################################################################
check_faa() {
    local dir="$1"
    local sp
    sp=$(basename "$dir")
    
    # 检查是否已检查过
    if grep -q "^${sp}$" "${CHECKED_FILE}" 2>/dev/null; then
        return 0
    fi
    
    log_running "检查样本: $sp"
    
    # 查找 .faa 文件
    local faa_file
    faa_file=$(find "$dir" -maxdepth 1 -type f -name "*.faa" 2>/dev/null | head -n 1)
    
    # 线程安全地记录检查状态
    if [ -n "$faa_file" ]; then
        {
            flock 200
            echo "$sp" >> "${HAS_FAA_FILE}"
            echo "$sp" >> "${CHECKED_FILE}"
        } 200>"${LOCK_FILE}"
        log_success "✓ $sp"
    else
        {
            flock 200
            echo "$sp" >> "${NO_FAA_FILE}"
            echo "$sp" >> "${CHECKED_FILE}"
        } 200>"${LOCK_FILE}"
        log_error "✗ $sp (未找到.faa文件)"
    fi
}

################################################################################
# 生成报告函数
################################################################################
generate_report() {
    local total_samples=0
    local has_faa=0
    local no_faa=0
    
    # 统计总数
    total_samples=$(find "${INPUT_DIR}" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l)
    
    # 统计有faa文件的数量
    if [ -f "${HAS_FAA_FILE}" ]; then
        has_faa=$(grep -c . "${HAS_FAA_FILE}" 2>/dev/null)
        if [ -z "$has_faa" ] || [ "$has_faa" -eq 0 ]; then
            has_faa=0
        fi
    else
        has_faa=0
    fi
    
    # 统计没有faa文件的数量
    if [ -f "${NO_FAA_FILE}" ]; then
        no_faa=$(grep -c . "${NO_FAA_FILE}" 2>/dev/null)
        if [ -z "$no_faa" ] || [ "$no_faa" -eq 0 ]; then
            no_faa=0
        fi
    else
        no_faa=0
    fi
    
    # 生成报告文件
    {
        echo "================================================================================"
        echo "                         .faa 文件检查报告"
        echo "================================================================================"
        echo ""
        echo "检查时间: $(date '+%Y-%m-%d %H:%M:%S')"
        echo "检查目录: ${INPUT_DIR}"
        echo ""
        echo "================================================================================"
        echo "                              统计信息"
        echo "================================================================================"
        echo "总样本数:            $total_samples"
        echo "存在.faa文件:        $has_faa"
        echo "不存在.faa文件:      $no_faa"
        echo "完成率:              $(awk "BEGIN {printf \"%.2f%%\", ($has_faa + $no_faa) / $total_samples * 100}")"
        echo ""
        
        if [ "$no_faa" -gt 0 ]; then
            echo "================================================================================"
            echo "                    以下样本不存在 .faa 文件"
            echo "================================================================================"
            sort "${NO_FAA_FILE}"
            echo ""
        fi
        
        if [ "$has_faa" -gt 0 ]; then
            echo "================================================================================"
            echo "                    以下样本存在 .faa 文件 ($has_faa 个)"
            echo "================================================================================"
            echo "(样本列表已省略，共 $has_faa 个样本)"
            echo ""
        fi
        
        echo "================================================================================"
        echo "报告生成时间: $(date '+%Y-%m-%d %H:%M:%S')"
        echo "================================================================================"
    } > "${REPORT_FILE}"
    
    # 打印到控制台
    cat "${REPORT_FILE}"
}

################################################################################
# 统计并显示汇总信息
################################################################################
print_summary() {
    local total_samples=0
    local has_faa=0
    local no_faa=0
    
    # 统计总数
    total_samples=$(find "${INPUT_DIR}" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l)
    
    # 统计完成数
    if [ -f "${HAS_FAA_FILE}" ]; then
        has_faa=$(grep -c . "${HAS_FAA_FILE}" 2>/dev/null)
        if [ -z "$has_faa" ] || [ "$has_faa" -eq 0 ]; then
            has_faa=0
        fi
    else
        has_faa=0
    fi
    
    # 统计失败数
    if [ -f "${NO_FAA_FILE}" ]; then
        no_faa=$(grep -c . "${NO_FAA_FILE}" 2>/dev/null)
        if [ -z "$no_faa" ] || [ "$no_faa" -eq 0 ]; then
            no_faa=0
        fi
    else
        no_faa=0
    fi
    
    echo ""
    echo -e "${CYAN}════════════════════════════════════════${NC}"
    echo -e "${CYAN}           检查结果统计${NC}"
    echo -e "${CYAN}════════════════════════════════════════${NC}"
    echo -e "总样本数:              ${BLUE}$total_samples${NC}"
    echo -e "存在.faa文件:          ${GREEN}$has_faa${NC}"
    echo -e "不存在.faa文件:        ${RED}$no_faa${NC}"
    echo -e "${CYAN}════════════════════════════════════════${NC}"
    echo ""
}

################################################################################
# 主程序
################################################################################
main() {
    log_info "========== .faa 文件检查开始 =========="
    log_info "检查目录: $INPUT_DIR"
    log_info "并行任务数: $JOBS"
    echo ""
    
    # 初始化进度追踪
    init_progress_tracking
    
    # 验证目录
    if [ ! -d "${INPUT_DIR}" ]; then
        log_error "输入目录不存在: $INPUT_DIR"
        exit 1
    fi
    
    mkdir -p "${OUTPUT_DIR}"
    
    # 导出函数和变量供parallel使用
    export INPUT_DIR OUTPUT_DIR PROGRESS_DIR
    export CHECKED_FILE HAS_FAA_FILE NO_FAA_FILE LOCK_FILE
    export RED GREEN YELLOW BLUE CYAN NC
    export -f check_faa log_info log_success log_warning log_error log_running
    
    # 使用parallel并行检查
    find "${INPUT_DIR}" -mindepth 1 -maxdepth 1 -type d -print0 2>/dev/null | \
        parallel --null \
                 --jobs "${JOBS}" \
                 --line-buffer \
                 --tag \
                 check_faa
    
    local parallel_exit_code=$?
    
    if [ $parallel_exit_code -eq 0 ]; then
        log_success "所有检查完成"
    elif [ $parallel_exit_code -eq 101 ]; then
        log_warning "检查被中断或出错，部分任务可能未完成"
    fi
    
    # 生成报告
    log_info "生成检查报告..."
    generate_report
    
    # 打印统计信息
    print_summary
    
    # 显示报告文件位置
    log_success "详细报告已保存到: ${REPORT_FILE}"
    
    log_info "========== .faa 文件检查结束 =========="
}

# 执行主程序
main "$@"
