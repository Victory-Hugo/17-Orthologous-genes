#!/bin/bash

################################################################################
# 脚本名称: 分析基因覆盖率和物种分布
# 作者: BigLin
# 功能: 分析单拷贝基因的物种覆盖率和物种的基因覆盖率
################################################################################

# ============================================================================
# 软件和依赖要求
# ============================================================================
# 必需软件:
#   - bash 4.0+
#   - GNU coreutils (find, wc, sort, uniq, awk, sed, grep)
#   - GNU parallel (用于并行处理)
#   - bc (用于数学计算)
# 
# 安装方法(Linux/Debian系统):
#   sudo apt-get install -y parallel bc
#
# macOS:
#   brew install parallel
#
# ============================================================================

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# ============================================================================
# 硬编码的变量配置
# ============================================================================
SINGLE_COPY_DIR="/mnt/f/15_Bam_Tam/2-物种树/output/all_single_copy"
OUTPUT_DIR="/mnt/f/15_Bam_Tam/2-物种树/output/coverage_analysis"
LOG_FILE="${OUTPUT_DIR}/analysis.log"
CHECKPOINT_DIR="${OUTPUT_DIR}/.checkpoint"
TEMP_DIR="${OUTPUT_DIR}/.temp"
NUM_JOBS=4  # 并行任务数
MAX_RETRIES=3  # 最大重试次数


mkdir -p "$OUTPUT_DIR"
mkdir -p "$CHECKPOINT_DIR"
mkdir -p "$TEMP_DIR"

# ============================================================================
# 工具函数
# ============================================================================

# 彩色打印函数
print_info() {
    echo -e "${BLUE}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$LOG_FILE"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$LOG_FILE"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$LOG_FILE"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$LOG_FILE"
}

# 进度条函数
show_progress() {
    local current=$1
    local total=$2
    local width=50
    local percentage=$((current * 100 / total))
    local completed=$((width * current / total))
    
    printf "\r${CYAN}["
    printf "%${completed}s" | tr ' ' '='
    printf "%$((width - completed))s" | tr ' ' '-'
    printf "]${NC} ${percentage}%% (${current}/${total})"
}

# 检查依赖
check_dependencies() {
    print_info "检查依赖软件..."
    
    local missing_deps=()
    
    for cmd in find wc sort uniq awk sed grep bc parallel; do
        if ! command -v "$cmd" &> /dev/null; then
            missing_deps+=("$cmd")
        fi
    done
    
    if [ ${#missing_deps[@]} -gt 0 ]; then
        print_error "缺少以下依赖: ${missing_deps[*]}"
        print_error "请安装依赖后重新运行"
        exit 1
    fi
    
    print_success "所有依赖已安装"
}

# 保存检查点
save_checkpoint() {
    local checkpoint_name=$1
    local data=$2
    echo "$data" > "${CHECKPOINT_DIR}/${checkpoint_name}"
    print_info "检查点已保存: ${checkpoint_name}"
}

# 恢复检查点
restore_checkpoint() {
    local checkpoint_name=$1
    if [ -f "${CHECKPOINT_DIR}/${checkpoint_name}" ]; then
        cat "${CHECKPOINT_DIR}/${checkpoint_name}"
        return 0
    fi
    return 1
}

# 检查检查点是否存在
has_checkpoint() {
    local checkpoint_name=$1
    [ -f "${CHECKPOINT_DIR}/${checkpoint_name}" ]
}

# 清理临时文件
cleanup_temp() {
    if [ -d "$TEMP_DIR" ]; then
        rm -rf "$TEMP_DIR"
    fi
}

# 清理未完成的任务
cleanup_on_exit() {
    print_warning "捕捉到中断信号，开始清理..."
    cleanup_temp
    print_info "清理完成，保留的检查点信息已保存"
    exit 130
}

# 设置信号处理
trap cleanup_on_exit SIGINT SIGTERM

# ============================================================================
# 主程序
# ============================================================================

main() {
    print_info "=========================================="
    print_info "开始分析基因覆盖率和物种分布"
    print_info "==========================================="
    
    check_dependencies
    
    # 验证输入目录
    if [ ! -d "$SINGLE_COPY_DIR" ]; then
        print_error "输入目录不存在: $SINGLE_COPY_DIR"
        exit 1
    fi
    
    # 第一步: 统计每个基因的物种数量
    step1_gene_species_count
    
    print_success "=========================================="
    print_success "分析完成！所有任务已完成"
    print_success "==========================================="
}

# 第一步: 统计每个基因的物种数量
step1_gene_species_count() {
    print_info "开始步骤1: 统计每个基因的物种数量..."
    
    local gene_count_file="${OUTPUT_DIR}/gene_species_count.csv"
    
    # 检查检查点
    if has_checkpoint "step1_done"; then
        print_warning "步骤1已完成，跳过此步"
        return 0
    fi
    
    echo "基因ID,物种数量" > "$gene_count_file"
    
    local total_genes=$(find "$SINGLE_COPY_DIR" -mindepth 1 -maxdepth 1 -type d | wc -l)
    local processed=0
    
    for gene_dir in "$SINGLE_COPY_DIR"/*; do
        if [ -d "$gene_dir" ]; then
            gene_id=$(basename "$gene_dir")
            if [ "$gene_id" != "README.txt" ]; then
                species_count=$(find "$gene_dir" -name "*.faa" -type f | wc -l)
                echo "$gene_id,$species_count" >> "$gene_count_file"
                
                processed=$((processed + 1))
                show_progress "$processed" "$total_genes"
            fi
        fi
    done
    
    echo ""
    save_checkpoint "step1_done" "true"
    print_success "步骤1完成: 已处理 $total_genes 个基因"
}

# ============================================================================
# 清理资源并启动主程序
# ============================================================================

# 清理前次运行的临时文件(可选)
# 如果需要强制重新计算，取消注释下一行
# rm -rf "$CHECKPOINT_DIR" "$TEMP_DIR"

main

# 最后清理临时目录
cleanup_temp

print_info "脚本执行完成"