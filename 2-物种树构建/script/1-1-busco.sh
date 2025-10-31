#!/bin/bash

################################################################################
# 使用BUSCO评估蛋白质序列的完整性
# 功能：
#   1. 批量运行BUSCO分析
#   2. 支持断点续跑
#   3. 支持彩色输出
#   4. 支持优雅停止（Ctrl+C）
#################################################################################

#*===========准备文件===========*#
#* download文件夹下：推荐使用prokka注释得到
# ├── GCA_000006905.1
# │   ├── GCA_000006905.1.faa
# │   ├── GCA_000006905.1.fasta
# │   ├── GCA_000006905.1.gff
# │   └── GCA_000006905.1_CDS.fna
# ├── GCA_000007925.1
# │   ├── GCA_000007925.1.faa
# │   ├── GCA_000007925.1.fasta
# │   ├── GCA_000007925.1.gff
# │   └── GCA_000007925.1_CDS.fna
#*===========DATABASE_DIR===========*#
#* busco_downloads文件夹下：下载得到
# /mnt/d/4-BUSCO-Database/busco_downloads
# ├── file_versions.tsv
# └── lineages
#     └── bacteria_odb10
#         ├── ancestral
#         ├── ancestral_variants
#         ├── dataset.cfg
#         ├── hmms
#         └── …
#! 先激活 conda activate busco

################################################################################
# 配置部分 - 所有变量写死
################################################################################
INPUT_DIR="/mnt/f/15_Bam_Tam/5-补齐更多物种/download"
OUTPUT_DIR="/mnt/f/15_Bam_Tam/5-补齐更多物种/output/busco_out"
DATABASE_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/2-物种树构建/database/busco_downloads/"
JOBS=7
LINEAGE="bacteria_odb10"
MODE="proteins"

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# 进度追踪相关
PROGRESS_DIR="${OUTPUT_DIR}/.progress"
LOCK_FILE="${PROGRESS_DIR}/lock"
COMPLETED_FILE="${PROGRESS_DIR}/completed.txt"
FAILED_FILE="${PROGRESS_DIR}/failed.txt"
RUNNING_PIDS_FILE="${PROGRESS_DIR}/running_pids.txt"

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
    echo -e "${RED}[ERROR]${NC} $*" >&2
}

log_running() {
    echo -e "${CYAN}[RUNNING]${NC} $*"
}

################################################################################
# 初始化函数
################################################################################
init_progress_tracking() {
    mkdir -p "${PROGRESS_DIR}"
    
    # 如果进度文件不存在，则创建
    if [ ! -f "${COMPLETED_FILE}" ]; then
        touch "${COMPLETED_FILE}"
    fi
    
    if [ ! -f "${FAILED_FILE}" ]; then
        touch "${FAILED_FILE}"
    fi
    
    if [ ! -f "${RUNNING_PIDS_FILE}" ]; then
        touch "${RUNNING_PIDS_FILE}"
    fi
}

################################################################################
# 清理资源函数
################################################################################
cleanup() {
    log_warning "捕获到中断信号，开始清理..."
    
    # 读取所有运行中的PID并杀死
    if [ -f "${RUNNING_PIDS_FILE}" ]; then
        while IFS= read -r pid; do
            if [ -n "$pid" ] && kill -0 "$pid" 2>/dev/null; then
                log_warning "杀死运行中的进程: PID=$pid"
                kill "$pid" 2>/dev/null || true
            fi
        done < "${RUNNING_PIDS_FILE}"
    fi
    
    # 等待所有后台进程完成
    wait 2>/dev/null || true
    
    # 清理锁文件
    rm -f "${LOCK_FILE}"
    
    log_warning "清理完成。已完成的任务已保留在 ${COMPLETED_FILE}"
    log_warning "未完成的任务已删除。"
}

# 设置信号处理
trap cleanup SIGINT SIGTERM

################################################################################
# 主要处理函数
################################################################################
run_busco() {
    local dir="$1"
    local sp
    sp=$(basename "$dir")
    
    # 检查是否已完成
    if grep -q "^${sp}$" "${COMPLETED_FILE}" 2>/dev/null; then
        log_success "跳过已完成的任务: $sp"
        return 0
    fi
    
    local faa
    faa=$(find "$dir" -type f -name "*.faa" | head -n 1)
    
    if [ -z "$faa" ]; then
        log_error "未找到 .faa 文件: $sp"
        echo "$sp" >> "${FAILED_FILE}"
        return 1
    fi
    
    log_running "开始处理: $sp (文件: $faa)"
    
    # 运行BUSCO（-f 强制覆盖已存在的结果）
    if busco -i "$faa" \
             -l "${LINEAGE}" \
             -o "$sp" \
             -m "${MODE}" \
             -f \
             --download_path "${DATABASE_DIR}" \
             --offline \
             --out_path "${OUTPUT_DIR}"; then
        
        log_success "完成: $sp"
        
        # 线程安全地记录完成状态
        {
            flock 200
            echo "$sp" >> "${COMPLETED_FILE}"
        } 200>"${LOCK_FILE}"
        
        return 0
    else
        log_error "失败: $sp"
        
        # 线程安全地记录失败状态
        {
            flock 200
            echo "$sp" >> "${FAILED_FILE}"
        } 200>"${LOCK_FILE}"
        
        return 1
    fi
}

################################################################################
# 统计结果
################################################################################
print_summary() {
    local total_samples=0
    local completed=0
    local failed=0
    local remaining=0
    
    # 统计总数
    total_samples=$(find "${INPUT_DIR}" -mindepth 1 -maxdepth 1 -type d | wc -l)
    
    # 统计完成数
    if [ -f "${COMPLETED_FILE}" ]; then
        completed=$(grep -c . "${COMPLETED_FILE}" 2>/dev/null)
        if [ -z "$completed" ] || [ "$completed" -eq 0 ]; then
            completed=0
        fi
    else
        completed=0
    fi
    
    # 统计失败数
    if [ -f "${FAILED_FILE}" ]; then
        failed=$(grep -c . "${FAILED_FILE}" 2>/dev/null)
        if [ -z "$failed" ] || [ "$failed" -eq 0 ]; then
            failed=0
        fi
    else
        failed=0
    fi
    
    # 计算未处理
    remaining=$((total_samples - completed - failed))
    
    echo ""
    echo -e "${CYAN}════════════════════════════════════════${NC}"
    echo -e "${CYAN}           任务统计总结${NC}"
    echo -e "${CYAN}════════════════════════════════════════${NC}"
    echo -e "总数:        ${BLUE}$total_samples${NC}"
    echo -e "已完成:      ${GREEN}$completed${NC}"
    echo -e "已失败:      ${RED}$failed${NC}"
    echo -e "未处理:      ${YELLOW}$remaining${NC}"
    echo -e "${CYAN}════════════════════════════════════════${NC}"
    echo ""
}

################################################################################
# 主程序
################################################################################
main() {
    log_info "========== BUSCO 批量分析开始 =========="
    log_info "输入目录: $INPUT_DIR"
    log_info "输出目录: $OUTPUT_DIR"
    log_info "数据库目录: $DATABASE_DIR"
    log_info "并行任务数: $JOBS"
    log_info "系统发育谱系: $LINEAGE"
    echo ""
    
    # 初始化进度追踪
    init_progress_tracking
    
    # 验证目录
    if [ ! -d "${INPUT_DIR}" ]; then
        log_error "输入目录不存在: $INPUT_DIR"
        exit 1
    fi
    
    if [ ! -d "${DATABASE_DIR}" ]; then
        log_error "数据库目录不存在: $DATABASE_DIR"
        exit 1
    fi
    
    cd "${INPUT_DIR}" || exit 1
    mkdir -p "${OUTPUT_DIR}"
    
    # 导出函数和变量供parallel使用
    export INPUT_DIR OUTPUT_DIR DATABASE_DIR JOBS LINEAGE MODE PROGRESS_DIR
    export COMPLETED_FILE FAILED_FILE LOCK_FILE RUNNING_PIDS_FILE
    export RED GREEN YELLOW BLUE CYAN NC
    export -f run_busco log_info log_success log_warning log_error log_running
    
    # 使用parallel并行处理，并记录PID便于后续清理
    find . -mindepth 1 -maxdepth 1 -type d -print0 | \
        parallel --null \
                 --jobs "${JOBS}" \
                 --line-buffer \
                 --tag \
                 --joblog "${PROGRESS_DIR}/parallel.log" \
                 run_busco
    
    local parallel_exit_code=$?
    
    if [ $parallel_exit_code -eq 0 ]; then
        log_success "所有任务完成"
    elif [ $parallel_exit_code -eq 101 ]; then
        log_error "任务被中断或出错，部分任务可能未完成"
    fi
    
    # 打印统计
    print_summary
    
    log_info "========== BUSCO 批量分析结束 =========="
}

# 执行主程序
main "$@"
