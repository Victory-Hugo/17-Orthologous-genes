#!/bin/bash

################################################################################
# 脚本名称: 筛选核心基因 (Core Genes Filter)
# 作者: BigLin
# 
# 脚本依赖项:
# - bash (>= 4.0)
# - awk (GNU awk 或 mawk)
# - bc (任意版本)
# - grep (GNU grep 或 BSD grep)
# - parallel (optional, 用于并行处理，推荐安装)
# - tput (用于彩色输出)
#
# 安装依赖项:
# Ubuntu/Debian: sudo apt-get install -y gawk bc parallel
# CentOS/RHEL: sudo yum install -y gawk bc parallel
# macOS: brew install awk bc parallel
#
################################################################################

set -euo pipefail

# ==================== 颜色和格式定义 ====================
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'  # No Color

# ==================== 变量定义 ====================
INPUT_DIR="/mnt/f/15_Bam_Tam/5-补齐更多物种/output/all_single_copy"
OUTPUT_DIR="/mnt/f/15_Bam_Tam/5-补齐更多物种/output/all_core_genomes"
COVERAGE_FILE="/mnt/f/15_Bam_Tam/5-补齐更多物种/output/coverage_analysis/gene_species_count.csv"
THRESHOLD=1200  #? 这个值代表在1200个物种中存在该基因
TEMP_DIR="/tmp/core_genes_$$"
MAX_PARALLEL_JOBS=8
PROCESSED_DIR="${OUTPUT_DIR}/.processed"
FAILED_DIR="${OUTPUT_DIR}/.failed"

# ==================== 函数定义 ====================

# 彩色打印函数
print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[✓]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# 进度条函数
show_progress_bar() {
    local current=$1
    local total=$2
    local width=50
    local percent=$((current * 100 / total))
    local filled=$((width * current / total))
    
    printf "\r${CYAN}[进度]${NC} ["
    printf "%${filled}s" | tr ' ' '='
    printf "%$((width - filled))s" | tr ' ' '-'
    printf "] ${percent}%% (${current}/${total})"
}

# 清理函数 - 处理用户中断
cleanup() {
    print_error "脚本被中断..."
    
    # 保存已完成的任务
    if [ -d "$PROCESSED_DIR" ] && [ "$(ls -A "$PROCESSED_DIR" 2>/dev/null)" ]; then
        print_info "正在保存已完成的任务..."
        find "$PROCESSED_DIR" -maxdepth 1 -type f -name "*.done" | while read -r done_file; do
            gene=$(basename "$done_file" .done)
            if [ -d "$TEMP_DIR/$gene" ]; then
                mv "$TEMP_DIR/$gene" "$OUTPUT_DIR/" 2>/dev/null || true
            fi
        done
        print_success "已完成的任务已保存到 $OUTPUT_DIR"
    fi
    
    # 删除临时目录
    rm -rf "$TEMP_DIR"
    
    # 清理失败的任务
    if [ -d "$FAILED_DIR" ]; then
        print_warning "以下基因复制失败，已保存到 $FAILED_DIR/failed_genes.txt"
        ls "$FAILED_DIR" 2>/dev/null | head -10
    fi
    
    exit 1
}

# 设置中断信号处理
trap cleanup SIGINT SIGTERM

# 并行复制函数
copy_gene_parallel() {
    local gene="$1"
    local src="$INPUT_DIR/$gene"
    local dst="$TEMP_DIR/$gene"
    local done_file="$PROCESSED_DIR/$gene.done"
    
    if [ -d "$src" ]; then
        mkdir -p "$dst"
        if cp -r "$src"/* "$dst/" 2>/dev/null; then
            touch "$done_file"
            echo "done"
        else
            mkdir -p "$FAILED_DIR"
            echo "$gene" >> "$FAILED_DIR/failed_genes.txt" 2>/dev/null || true
            echo "fail"
        fi
    else
        mkdir -p "$FAILED_DIR"
        echo "$gene" >> "$FAILED_DIR/failed_genes.txt" 2>/dev/null || true
        echo "notfound"
    fi
}

export -f copy_gene_parallel
export INPUT_DIR TEMP_DIR PROCESSED_DIR FAILED_DIR

# ==================== 主程序 ====================

print_info "=========================================="
print_info "核心基因筛选工具"
print_info "作者: BigLin"
print_info "=========================================="

# 创建输出目录
mkdir -p "$OUTPUT_DIR"
mkdir -p "$TEMP_DIR"
mkdir -p "$PROCESSED_DIR"

# 检查输入文件
if [ ! -f "$COVERAGE_FILE" ]; then
    print_error "覆盖率文件不存在: $COVERAGE_FILE"
    exit 1
fi

if [ ! -d "$INPUT_DIR" ]; then
    print_error "输入目录不存在: $INPUT_DIR"
    exit 1
fi

print_info "开始筛选核心基因..."
print_info "覆盖率阈值: $THRESHOLD 个物种"
print_info "输入目录: $INPUT_DIR"
print_info "输出目录: $OUTPUT_DIR"

# 筛选满足条件的基因
print_info "正在筛选基因..."
CORE_GENES=$(awk -F',' -v threshold="$THRESHOLD" 'NR>1 && $2>=threshold {print $1}' "$COVERAGE_FILE" 2>/dev/null | sort)

# 统计核心基因数量
CORE_COUNT=$(echo "$CORE_GENES" | grep -c . || echo "0")

if [ "$CORE_COUNT" -eq 0 ]; then
    print_warning "未找到满足条件的核心基因"
    exit 0
fi

print_success "找到 $CORE_COUNT 个核心基因（覆盖率 ≥ $THRESHOLD 个物种）"

# 复制核心基因目录 - 使用并行处理
print_info "正在使用并行处理复制核心基因 (最多 $MAX_PARALLEL_JOBS 个并发任务)..."

# 检查是否安装了 GNU parallel
if command -v parallel &> /dev/null; then
    # 使用 GNU parallel
    count=0
    echo "$CORE_GENES" | parallel --no-notice -j "$MAX_PARALLEL_JOBS" \
        "copy_gene_parallel {} && echo -ne '\\r已处理: $((++count))/$CORE_COUNT' || true" 2>/dev/null || true
    
    count=$(find "$PROCESSED_DIR" -name "*.done" | wc -l)
else
    # 回退到顺序处理，但仍显示进度条
    print_warning "未安装 GNU parallel，使用顺序处理（推荐安装 parallel 以加快速度）"
    count=0
    total=$(echo "$CORE_GENES" | wc -l)
    
    while IFS= read -r gene; do
        if [ -n "$gene" ]; then
            show_progress_bar "$count" "$total"
            copy_gene_parallel "$gene" > /dev/null 2>&1 || true
            count=$((count + 1))
        fi
    done <<< "$CORE_GENES"
    
    echo ""  # 换行
fi

# 将临时目录的文件移动到最终输出目录
print_info "正在转移完成的基因到输出目录..."
mv "$TEMP_DIR"/* "$OUTPUT_DIR/" 2>/dev/null || true
count=$(ls -1 "$OUTPUT_DIR" | grep -v "^\." | wc -l)

print_success "核心基因复制完成！"
print_info "源目录: $INPUT_DIR"
print_info "目标目录: $OUTPUT_DIR"
print_info "总基因数: $(ls -1 "$INPUT_DIR" | wc -l)"
print_info "核心基因数: $count"
print_info "覆盖率阈值: ≥ $THRESHOLD 个物种"

# 生成核心基因列表
print_info "正在生成核心基因列表..."
ls -1 "$OUTPUT_DIR" | grep -v "^\." > "$OUTPUT_DIR/core_genes_list.txt"
print_success "核心基因列表已保存到: $OUTPUT_DIR/core_genes_list.txt"

# 生成统计报告
print_info "正在生成统计报告..."
{
    echo "=============================================="
    echo "核心基因筛选统计报告"
    echo "=============================================="
    echo "脚本作者: BigLin"
    echo "筛选日期: $(date '+%Y-%m-%d %H:%M:%S')"
    echo ""
    echo "配置参数:"
    echo "  输入目录: $INPUT_DIR"
    echo "  输出目录: $OUTPUT_DIR"
    echo "  覆盖率阈值: ≥ $THRESHOLD 个物种"
    echo ""
    echo "统计结果:"
    echo "  总基因数: $(ls -1 "$INPUT_DIR" | wc -l)"
    echo "  核心基因数: $count"
    echo "  筛选比例: $(echo "scale=2; $count * 100 / $(ls -1 "$INPUT_DIR" | wc -l)" | bc)%"
    echo ""
    echo "核心基因列表:"
    echo "=============================================="
    cat "$OUTPUT_DIR/core_genes_list.txt"
    echo "=============================================="
} > "$OUTPUT_DIR/core_genes_report.txt"

print_success "统计报告已保存到: $OUTPUT_DIR/core_genes_report.txt"

# 清理临时文件
rm -rf "$TEMP_DIR"

print_success "筛选完成！"