#!/bin/bash

################################################################################
# 脚本名称: 收集单拷贝基因序列
# 作者: BigLin
# 描述: 从BUSCO输出目录中并行收集单拷贝基因序列，支持断点续跑和进度显示
#
# 依赖软件和包:
#   - bash >= 4.0
#   - GNU Parallel (gnu-parallel)
#   - GNU Coreutils (cp, ls, mkdir, basename)
#   - sed
#   - grep
#   
# 安装依赖 (Ubuntu/Debian):
#   sudo apt-get install parallel coreutils
#
# 安装依赖 (CentOS/RHEL):
#   sudo yum install parallel coreutils
#
################################################################################

set -euo pipefail

# ==================== 颜色定义 ====================
readonly RED='\033[0;31m'
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'
readonly BLUE='\033[0;34m'
readonly CYAN='\033[0;36m'
readonly NC='\033[0m' # No Color

# ==================== 变量定义（全部写死） ====================
readonly BUSCO_OUTPUT_DIR="/mnt/f/15_Bam_Tam/5-补齐更多物种/output/busco_out"
readonly OUTPUT_DIR="/mnt/f/15_Bam_Tam/5-补齐更多物种/output/all_single_copy"
readonly BUSCO_RUN_NAME="run_bacteria_odb10"
readonly SINGLE_COPY_SUBDIR="busco_sequences/single_copy_busco_sequences"
readonly LOG_DIR="/mnt/f/15_Bam_Tam/5-补齐更多物种/output/logs"
readonly CHECKPOINT_FILE="${LOG_DIR}/.collect_checkpoint"
readonly TEMP_DIR="${OUTPUT_DIR}/.tmp"
readonly MAX_JOBS=8
readonly BATCH_SIZE=100

# ==================== 彩色打印函数 ====================
print_info() {
    echo -e "${BLUE}[INFO]${NC} $*"
}

print_success() {
    echo -e "${GREEN}[✓]${NC} $*"
}

print_warning() {
    echo -e "${YELLOW}[⚠]${NC} $*"
}

print_error() {
    echo -e "${RED}[✗]${NC} $*" >&2
}

print_progress() {
    local current=$1
    local total=$2
    local percent=$((current * 100 / total))
    local bar_length=40
    local filled=$((bar_length * current / total))
    local empty=$((bar_length - filled))
    
    printf "\r${CYAN}进度:${NC} ["
    printf "%${filled}s" | tr ' ' '='
    printf "%${empty}s" | tr ' ' '-'
    printf "] ${percent}%% (${current}/${total})"
}

# ==================== 检查依赖 ====================
check_dependencies() {
    print_info "检查依赖软件..."
    
    local missing_deps=()
    
    for cmd in parallel cp ls mkdir basename sed grep; do
        if ! command -v "$cmd" &> /dev/null; then
            missing_deps+=("$cmd")
        fi
    done
    
    if [ ${#missing_deps[@]} -gt 0 ]; then
        print_error "缺少以下依赖: ${missing_deps[*]}"
        print_error "请运行: sudo apt-get install parallel coreutils"
        exit 1
    fi
    
    print_success "所有依赖检查完成"
}

# ==================== 信号处理与清理 ====================
cleanup() {
    local exit_code=$?
    
    if [ $exit_code -ne 0 ]; then
        print_warning "脚本被中断，清理未完成的任务..."
        
        # 杀死所有parallel进程
        pkill -P $$ parallel 2>/dev/null || true
        wait 2>/dev/null || true
        
        # 删除临时目录
        if [ -d "$TEMP_DIR" ]; then
            rm -rf "$TEMP_DIR"
            print_warning "已删除临时目录"
        fi
        
        print_warning "已保留已完成的任务在: $OUTPUT_DIR"
    fi
    
    exit $exit_code
}

trap cleanup EXIT
trap 'exit 130' INT TERM

# ==================== 初始化 ====================
initialize() {
    print_info "初始化环境..."
    
    mkdir -p "$OUTPUT_DIR"
    mkdir -p "$LOG_DIR"
    mkdir -p "$TEMP_DIR"
    
    # 检查输入目录
    if [ ! -d "$BUSCO_OUTPUT_DIR" ]; then
        print_error "BUSCO输出目录不存在: $BUSCO_OUTPUT_DIR"
        exit 1
    fi
    
    print_success "环境初始化完成"
}

# ==================== 获取基因列表 ====================
get_gene_ids() {
    local single_copy_dir
    local first_sample
    local samples=()

    mapfile -t samples < <(find "$BUSCO_OUTPUT_DIR" -mindepth 1 -maxdepth 1 -type d ! -name '.*' | sort)

    first_sample=${samples[0]:-}

    if [ -z "$first_sample" ]; then
        print_error "BUSCO输出目录为空: $BUSCO_OUTPUT_DIR" >&2
        exit 1
    fi

    single_copy_dir="$first_sample/$BUSCO_RUN_NAME/$SINGLE_COPY_SUBDIR"

    if [ ! -d "$single_copy_dir" ]; then
        print_error "找不到单拷贝基因目录: $single_copy_dir" >&2
        exit 1
    fi

    find "$single_copy_dir" -maxdepth 1 -type f -name '*.faa' -printf '%f\n' 2>/dev/null | \
        sed 's/\.faa$//' | \
        sort || true
}

# ==================== 并行收集基因序列 ====================
collect_gene_sequences() {
    local gene_id="$1"
    local gene_dir="$OUTPUT_DIR/$gene_id"
    local temp_gene_dir="$TEMP_DIR/$gene_id"
    local count=0
    
    mkdir -p "$temp_gene_dir"
    
    for sample_dir in "$BUSCO_OUTPUT_DIR"/*; do
        if [ -d "$sample_dir" ]; then
            local sample_name
            sample_name=$(basename "$sample_dir")
            local single_copy_file="$sample_dir/$BUSCO_RUN_NAME/$SINGLE_COPY_SUBDIR/${gene_id}.faa"
            
            if [ -f "$single_copy_file" ]; then
                cp "$single_copy_file" "$temp_gene_dir/${sample_name}.faa"
                ((count++))
            fi
        fi
    done
    
    echo "$gene_id|$count"
}

export -f collect_gene_sequences
export OUTPUT_DIR BUSCO_OUTPUT_DIR BUSCO_RUN_NAME SINGLE_COPY_SUBDIR TEMP_DIR

# ==================== 主函数 ====================
main() {
    print_info "========================================"
    print_info "    单拷贝基因序列收集工具 (作者: BigLin)"
    print_info "========================================"
    echo ""
    
    # 1. 检查依赖
    check_dependencies
    
    # 2. 初始化
    initialize
    
    # 3. 获取基因列表
    print_info "扫描基因ID列表..."
    local gene_ids_array
    mapfile -t gene_ids_array < <(get_gene_ids)
    local total_genes=${#gene_ids_array[@]}
    
    # 检查是否获取到基因列表
    if [ "$total_genes" -eq 0 ]; then
        print_warning "未发现任何单拷贝基因，脚本退出"
        print_success "输出目录已创建: $OUTPUT_DIR"
        echo ""
        print_success "✓ 任务完成！"
        return 0
    fi
    
    print_info "发现 ${total_genes} 个单拷贝基因"
    echo ""
    
    # 4. 并行收集基因序列
    print_info "开始并行收集基因序列 (最多 ${MAX_JOBS} 个并行任务)..."
    
    local results_file="$TEMP_DIR/results.txt"
    > "$results_file"
    
    printf '%s\n' "${gene_ids_array[@]}" | \
    parallel --halt soon,fail=1 \
             --jobs "$MAX_JOBS" \
             --progress \
             --line-buffer \
             --tag \
             collect_gene_sequences >> "$results_file"
    
    echo ""
    print_success "并行收集完成"
    
    # 5. 将临时文件移到最终位置
    print_info "整理输出文件..."
    if [ ${#gene_ids_array[@]} -gt 0 ]; then
        for gene_id in "${gene_ids_array[@]}"; do
            if [ -d "$TEMP_DIR/$gene_id" ]; then
                mkdir -p "$OUTPUT_DIR/$gene_id"
                mv "$TEMP_DIR/$gene_id"/* "$OUTPUT_DIR/$gene_id/" 2>/dev/null || true
            fi
        done
    fi
    
    # 6. 统计结果
    print_info "生成统计信息..."
    echo ""
    print_info "======== 统计结果 ========"
    
    local total_sequences=0
    local genes_with_sequences=0
    
    # 安全地遍历输出目录
    if [ -d "$OUTPUT_DIR" ] && [ "$(ls -A "$OUTPUT_DIR" 2>/dev/null | grep -v '^\.' | wc -l)" -gt 0 ]; then
        while IFS= read -r gene_dir; do
            if [ -d "$gene_dir" ]; then
                local gene_id
                gene_id=$(basename "$gene_dir")
                local count
                count=$(find "$gene_dir" -maxdepth 1 -type f -name '*.faa' 2>/dev/null | wc -l)
                
                if [ "$count" -gt 0 ]; then
                    ((genes_with_sequences++))
                    ((total_sequences += count))
                    print_success "  $gene_id: $count 个序列"
                fi
            fi
        done < <(find "$OUTPUT_DIR" -maxdepth 1 -type d ! -name '.*' 2>/dev/null)
    else
        print_warning "输出目录为空，未找到任何基因序列"
    fi
    
    echo ""
    print_info "====== 总体统计 ======"
    print_success "包含序列的基因数: $genes_with_sequences"
    print_success "总序列数: $total_sequences"
    print_success "输出目录: $OUTPUT_DIR"
    print_success "日志目录: $LOG_DIR"
    
    # 7. 清理临时文件
    print_info "清理临时文件..."
    if [ -d "$TEMP_DIR" ]; then
        rm -rf "$TEMP_DIR"
    fi
    
    echo ""
    print_success "✓ 任务完成！"
}

# ==================== 执行主函数 ====================
main "$@"
