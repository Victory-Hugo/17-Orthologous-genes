#!/bin/bash

################################################################################
# 脚本名称: HMM蛋白序列比对工具（极简for循环版）
# 作者: BigLin
# 功能: 使用HMMER进行HMM比对，双重for循环直接跑全部组合，不做输入检查
################################################################################

set -o pipefail

# ============================================================================
# 配置变量（全部写死）
# ============================================================================
LIST_FILE="/data_raid/7_luolintao/2-NCBI-2024/1-All-data/1-faa_list.txt" # 每行一个FAA文件绝对路径
HMM_FILE_TXT="/home/luolintao/0_Github/17-Orthologous-genes/1-HMM/conf/2-hmm库_list.txt" # 每行一个HMM文件绝对路径
OUTPUT_DIR="/home/luolintao/7-Rv0194/output/"
LOG_FILE="${OUTPUT_DIR}processing.log"
FAILED_LOG_FILE="${OUTPUT_DIR}failed_tasks.log"
CPU_PER_JOB="2"
E_VALUE="1e-5"

# 存在缺失矩阵脚本及输出
MATRIX_SCRIPT="/home/luolintao/0_Github/17-Orthologous-genes/1-HMM/script/3-生成存在缺失矩阵.py"
PRESENCE_MATRIX_TSV="${OUTPUT_DIR}/presence_absence_matrix.tsv"
MERGED_TBL_TSV="${OUTPUT_DIR}/tbl_merge.tsv"
MATRIX_EVALUE_CUTOFF="1e-10"

declare -a FAA_FILES=()
declare -a HMM_FILES=()
TOTAL_FAA_FILES=0
TOTAL_HMM_FILES=0
TOTAL_TASKS=0

# ============================================================================
# 颜色 & 打印函数
# ============================================================================
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
NC='\033[0m'

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

# 初始化目录
initialize_dirs() {
    mkdir -p "$OUTPUT_DIR" || { echo "[ERROR] 无法创建输出目录: $OUTPUT_DIR"; exit 1; }
    : > "$LOG_FILE"
    : > "$FAILED_LOG_FILE"
}

# 检查依赖工具（只保留 hmmsearch）
check_dependencies() {
    print_info "检查依赖软件..."
    local missing=0

    if ! command -v hmmsearch &> /dev/null; then
        print_error "未找到 hmmsearch，请安装 HMMER"
        missing=1
    fi

    if [[ $missing -eq 1 ]]; then
        exit 1
    fi

    print_success "HMMER 检查通过"
}

# ============================================================================
# 主程序
# ============================================================================

print_header "HMM序列比对工具 v1.0（极简for循环版） - 作者: BigLin"

initialize_dirs
check_dependencies

# 不“检查输入文件”，直接读列表并开跑
print_info "读取 FAA 列表: $LIST_FILE"
while IFS= read -r line || [[ -n "$line" ]]; do
    line="${line%$'\r'}"
    # 去掉首尾空白 & 注释、空行
    [[ "$line" =~ ^[[:space:]]*$ ]] && continue
    trimmed=$(printf '%s\n' "$line" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
    [[ "$trimmed" == \#* ]] && continue
    FAA_FILES+=("$trimmed")
done < "$LIST_FILE"

print_info "读取 HMM 列表: $HMM_FILE_TXT"
while IFS= read -r line || [[ -n "$line" ]]; do
    line="${line%$'\r'}"
    [[ "$line" =~ ^[[:space:]]*$ ]] && continue
    trimmed=$(printf '%s\n' "$line" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
    [[ "$trimmed" == \#* ]] && continue
    HMM_FILES+=("$trimmed")
done < "$HMM_FILE_TXT"

TOTAL_FAA_FILES=${#FAA_FILES[@]}
TOTAL_HMM_FILES=${#HMM_FILES[@]}
TOTAL_TASKS=$((TOTAL_FAA_FILES * TOTAL_HMM_FILES))

print_info "FAA 文件数: $TOTAL_FAA_FILES"
print_info "HMM 文件数: $TOTAL_HMM_FILES"
print_info "总任务数: $TOTAL_TASKS"

task_idx=0

for hmm_file in "${HMM_FILES[@]}"; do
    hmm_filename=$(basename "$hmm_file")
    hmm_basename="${hmm_filename%.*}"
    output_subdir="${OUTPUT_DIR}/${hmm_basename}"
    mkdir -p "$output_subdir"

    for faa_file in "${FAA_FILES[@]}"; do
        faa_filename=$(basename "$faa_file")
        faa_basename="${faa_filename%.*}"
        output_file="${output_subdir}/${faa_basename}.tbl"
        err_log="${output_subdir}/${faa_basename}.err"

        ((task_idx++))
        print_info "任务 ${task_idx}/${TOTAL_TASKS}: HMM=${hmm_basename}, FAA=${faa_basename}"
        print_info "  hmmsearch 输入: $hmm_file  vs  $faa_file"

        if hmmsearch --cpu "$CPU_PER_JOB" \
            -E "$E_VALUE" \
            --tblout "$output_file" \
            "$hmm_file" \
            "$faa_file" 2> "$err_log"; then

            if [[ -s "$output_file" ]]; then
                print_success "完成: ${hmm_basename} vs ${faa_basename} -> $output_file"
                rm -f "$err_log"
            else
                print_warning "任务失败(无输出): ${hmm_basename} vs ${faa_basename}"
                echo "[FAILED] ${hmm_basename}|${faa_basename} - 输出文件未生成" >> "$FAILED_LOG_FILE"
                [[ -s "$err_log" ]] && cat "$err_log" >> "$FAILED_LOG_FILE"
            fi
        else
            print_warning "任务失败(hmmsearch错误): ${hmm_basename} vs ${faa_basename}"
            echo "[FAILED] ${hmm_basename}|${faa_basename} - hmmsearch返回错误" >> "$FAILED_LOG_FILE"
            [[ -s "$err_log" ]] && cat "$err_log" >> "$FAILED_LOG_FILE"
            rm -f "$output_file"
        fi
    done
done

# 生成存在/缺失矩阵
if [[ -f "$MATRIX_SCRIPT" ]]; then
    print_info "开始生成存在/缺失矩阵和合并表..."
    if python3 "$MATRIX_SCRIPT" \
        --base-dir "$OUTPUT_DIR" \
        --output-tsv "$PRESENCE_MATRIX_TSV" \
        --merged-tsv "$MERGED_TBL_TSV" \
        --evalue-cutoff "$MATRIX_EVALUE_CUTOFF"; then
        print_success "存在/缺失矩阵生成完成: $PRESENCE_MATRIX_TSV"
        print_success "合并tbl表生成完成: $MERGED_TBL_TSV"
    else
        print_warning "存在/缺失矩阵生成失败，请检查脚本和日志。"
    fi
else
    print_warning "未找到存在缺失矩阵脚本: $MATRIX_SCRIPT"
fi

print_header "处理完成"
print_success "所有任务处理完成！"
print_info "输出目录: $OUTPUT_DIR"
print_info "日志文件: $LOG_FILE"
[[ -s "$FAILED_LOG_FILE" ]] && print_warning "失败日志: $FAILED_LOG_FILE"

echo ""
print_success "脚本执行成功！"
python3 /home/luolintao/test_mail.py "2-HMM/script/2-HMM比对.sh任务完成通知" "<p>2-HMM/script/2-HMM比对.sh分析已完成，请查看结果目录。</p>"
