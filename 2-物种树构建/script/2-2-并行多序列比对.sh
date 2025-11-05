#!/bin/bash

################################################################################
# 并行多序列比对脚本 - 使用GNU Parallel加速处理
################################################################################
# 作者: BigLin
# 日期: 2025-10-31
#
# ===================== 脚本运行需要的软件和packages =====================
# 1. GNU Parallel - 用于并行处理
#    安装: sudo apt-get install parallel
#    或: brew install parallel
#
# 2. MAFFT - 多序列比对工具
#    安装: sudo apt-get install mafft
#    或: brew install mafft
#
# 3. standard Unix tools: find, sed, grep, cat, rm
#
# =========================================================================

# 设置路径 - 所有变量写死
CORE_GENES_DIR="/mnt/f/15_Bam_Tam/5-补齐更多物种/output/all_core_genomes"
OUTPUT_DIR="/mnt/f/15_Bam_Tam/5-补齐更多物种/output/alignments"
ALL_SPECIES_FILE="/mnt/f/15_Bam_Tam/2-物种树/conf/all_species.txt"
TEMP_GENES_LIST="/tmp/core_genes_list_$$.txt"
CHECKPOINT_DIR="/mnt/f/15_Bam_Tam/2-物种树/output/.checkpoint"
COMPLETED_GENES_FILE="$CHECKPOINT_DIR/completed_genes.txt"
PROGRESS_FILE="$CHECKPOINT_DIR/progress.log"

# 并行参数
THREADS=8

# 彩色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# 创建输出目录和checkpoint目录
mkdir -p "$OUTPUT_DIR" "$CHECKPOINT_DIR"

# 设置trap处理Ctrl+C中断
trap cleanup_on_interrupt SIGINT SIGTERM

# 清理函数 - 删除未完成的任务，保留完成的任务
cleanup_on_interrupt() {
    echo -e "\n${RED}[ERROR]${NC} 脚本被中断，正在清理未完成的任务..."
    
    # 获取已完成的基因列表
    local completed_genes=()
    if [[ -f "$COMPLETED_GENES_FILE" ]]; then
        mapfile -t completed_genes < "$COMPLETED_GENES_FILE"
    fi
    
    # 找出所有已生成的对齐文件
    local all_aligned_files=$(find "$OUTPUT_DIR" -name "*_aligned.faa" -type f)
    
    # 删除未完成的临时文件
    for file in $all_aligned_files; do
        local gene_name=$(basename "$file" | sed 's/_aligned.faa//')
        local is_completed=0
        
        for completed_gene in "${completed_genes[@]}"; do
            if [[ "$completed_gene" == "$gene_name" ]]; then
                is_completed=1
                break
            fi
        done
        
        if [[ $is_completed -eq 0 ]]; then
            echo -e "${YELLOW}[DELETE]${NC} 删除未完成的文件: $file"
            rm -f "$file"
        fi
    done
    
    # 删除临时文件
    rm -f "$TEMP_GENES_LIST"
    
    echo -e "${GREEN}[INFO]${NC} 清理完成。已完成的任务已保留。"
    exit 130
}

# 打印带颜色的消息
print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# 打印分隔线
print_separator() {
    echo -e "${CYAN}===========================================${NC}"
}

print_separator
echo -e "${CYAN}并行核心基因多序列比对流程${NC}"
print_separator
print_info "开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
print_info "使用CPU核心数: $THREADS"
print_info "输出目录: $OUTPUT_DIR"
print_separator

print_info "Step 1/5: 获取完整物种列表..."
find /mnt/f/15_Bam_Tam/2-物种树/output/all_single_copy -name "*.faa" | sed 's|.*/||; s|\.faa$||' | sort -u > "$ALL_SPECIES_FILE"
if [[ -f "$ALL_SPECIES_FILE" ]]; then
    local species_count=$(wc -l < "$ALL_SPECIES_FILE")
    print_success "已获取 $species_count 个物种"
else
    print_error "无法生成物种列表"
    exit 1
fi

print_info "Step 2/5: 获取核心基因列表..."
ls -1 "$CORE_GENES_DIR" | grep -E "^[0-9].*at2$" > "$TEMP_GENES_LIST"
local gene_count=$(wc -l < "$TEMP_GENES_LIST")
print_success "已找到 $gene_count 个核心基因"

# 检查checkpoint，识别已完成的基因
print_info "Step 3/5: 检查断点续跑信息..."
if [[ -f "$COMPLETED_GENES_FILE" ]]; then
    local completed_count=$(wc -l < "$COMPLETED_GENES_FILE")
    print_success "发现 $completed_count 个已完成的基因，将跳过这些基因"
    # 生成待处理的基因列表（排除已完成的）
    grep -vFf "$COMPLETED_GENES_FILE" "$TEMP_GENES_LIST" > "${TEMP_GENES_LIST}.todo" || true
    mv "${TEMP_GENES_LIST}.todo" "$TEMP_GENES_LIST"
    local remaining_count=$(wc -l < "$TEMP_GENES_LIST")
    print_info "还需处理 $remaining_count 个基因"
else
    print_info "首次运行，将处理所有基因"
fi

echo ""
print_info "Step 4/5: 定义处理函数..."

# 导出函数和变量供parallel使用
export CORE_GENES_DIR OUTPUT_DIR ALL_SPECIES_FILE CHECKPOINT_DIR COMPLETED_GENES_FILE
export RED GREEN YELLOW BLUE CYAN NC

# 定义处理单个基因的函数
process_gene() {
    local gene=$1
    
    local gene_dir="$CORE_GENES_DIR/$gene"
    local output_file="$OUTPUT_DIR/${gene}.faa"
    local aligned_file="$OUTPUT_DIR/${gene}_aligned.faa"
    
    # 创建包含所有物种的完整序列文件
    > "$output_file"  # 清空文件
    
    while read -r species; do
        local species_file="$gene_dir/${species}.faa"
        if [[ -f "$species_file" ]]; then
            # 物种有该基因，直接添加
            cat "$species_file" >> "$output_file"
        else
            # 物种缺失该基因，添加空序列
            echo ">${species}" >> "$output_file"
            echo "X" >> "$output_file"  # 使用X作为占位符
        fi
    done < "$ALL_SPECIES_FILE"
    
    # 使用MAFFT进行多序列比对
    local start_time=$(date +%s)
    mafft --auto --quiet "$output_file" > "$aligned_file" 2>/dev/null
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    # 标记基因为已完成，添加到checkpoint文件
    echo "$gene" >> "$COMPLETED_GENES_FILE"
    
    # 清理临时文件
    rm -f "$output_file"
    
    # 输出处理结果（支持并行的非阻塞输出）
    echo "[DONE] $gene (${duration}s)"
}

# 导出函数
export -f process_gene

print_info "Step 5/5: 使用GNU Parallel并行处理所有基因 (共 $(wc -l < "$TEMP_GENES_LIST") 个)..."
echo ""

# 使用parallel处理所有基因，显示彩色进度条
parallel -j "$THREADS" \
    --progress \
    --bar \
    --eta \
    --line-buffer \
    process_gene :::: "$TEMP_GENES_LIST" 2>&1 | tee -a "$PROGRESS_FILE"

echo ""
print_separator
print_success "并行多序列比对流程完成！"
print_info "完成时间: $(date '+%Y-%m-%d %H:%M:%S')"
print_info "输出目录: $OUTPUT_DIR"
print_info "已完成基因数: $(wc -l < "$COMPLETED_GENES_FILE" 2>/dev/null || echo 0)"

# 统计输出文件
local aligned_count=$(find "$OUTPUT_DIR" -name "*_aligned.faa" -type f | wc -l)
print_success "已生成 $aligned_count 个比对文件"

print_separator

# 清理临时文件
rm -f "$TEMP_GENES_LIST"

print_info "临时文件已清理"
print_separator