#!/bin/bash

# ============================================================================
# 脚本说明: 整理xml文件到指定文件夹
# 功能: 将meta目录中的所有.xml文件移动到meta/xml子文件夹中
# 日期: 2025-10-30
# ============================================================================
# ============================================================================
# 变量定义（全部写死）
# ============================================================================
META_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-数据获取/2-下载特定物种序列/meta"
XML_SUBDIR="${META_DIR}/xml"
PROGRESS_FILE="/tmp/organize_xml_progress.txt"
TEMP_FILE="/tmp/organize_xml_temp.txt"
MAX_PARALLEL=4

# ============================================================================
# 颜色定义
# ============================================================================
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color


# ============================================================================
# 清理函数
# ============================================================================
cleanup() {
    echo -e "\n${YELLOW}[警告]${NC} 脚本被中断，清理未完成任务..."
    
    # 杀死所有后台进程
    jobs -p | xargs -r kill 2>/dev/null
    
    # 等待所有后台进程完成
    wait 2>/dev/null
    
    # 清理临时文件
    rm -f "${TEMP_FILE}"
    
    echo -e "${YELLOW}[完成]${NC} 已清理所有未完成任务，已完成的任务已保留"
    exit 0
}

# ============================================================================
# 信号处理
# ============================================================================
trap cleanup SIGINT SIGTERM

# ============================================================================
# 彩色打印函数
# ============================================================================
print_success() {
    echo -e "${GREEN}✓ [成功]${NC} $1"
}

print_error() {
    echo -e "${RED}✗ [错误]${NC} $1"
}

print_info() {
    echo -e "${BLUE}ℹ [信息]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}⚠ [警告]${NC} $1"
}

# ============================================================================
# 断点续跑机制
# ============================================================================
load_progress() {
    if [[ -f "${PROGRESS_FILE}" ]]; then
        print_info "检测到上次进度，正在恢复..."
        cat "${PROGRESS_FILE}"
    fi
}

save_progress() {
    local file="$1"
    echo "$file" >> "${PROGRESS_FILE}"
}

is_file_processed() {
    local file="$1"
    if [[ -f "${PROGRESS_FILE}" ]]; then
        grep -q "^${file}$" "${PROGRESS_FILE}"
        return $?
    fi
    return 1
}

# ============================================================================
# 移动XML文件函数
# ============================================================================
move_xml_file() {
    local xml_file="$1"
    local xml_name=$(basename "$xml_file")
    
    # 检查是否已处理过
    if is_file_processed "$xml_name"; then
        print_info "已跳过（已处理）: $xml_name"
        return 0
    fi
    
    # 执行移动操作
    if mv "$xml_file" "${XML_SUBDIR}/" 2>/dev/null; then
        print_success "已移动: $xml_name"
        save_progress "$xml_name"
        return 0
    else
        print_error "移动失败: $xml_name"
        return 1
    fi
}

# ============================================================================
# 主函数
# ============================================================================
main() {
    echo -e "${BLUE}╔════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║${NC}  XML文件整理脚本 v1.0                ${BLUE}║${NC}"
    echo -e "${BLUE}╚════════════════════════════════════════╝${NC}"
    echo ""
    
    # 检查meta目录是否存在
    if [[ ! -d "${META_DIR}" ]]; then
        print_error "meta目录不存在: ${META_DIR}"
        exit 1
    fi
    
    print_info "源目录: ${META_DIR}"
    print_info "目标目录: ${XML_SUBDIR}"
    echo ""
    
    # 创建xml子文件夹
    if mkdir -p "${XML_SUBDIR}" 2>/dev/null; then
        print_success "已创建或确认xml文件夹"
    else
        print_error "创建xml文件夹失败"
        exit 1
    fi
    
    # 加载上次进度
    load_progress
    echo ""
    
    # 获取所有xml文件
    local xml_files=()
    while IFS= read -r -d '' file; do
        xml_files+=("$file")
    done < <(find "${META_DIR}" -maxdepth 1 -name "*.xml" -type f -print0)
    
    # 检查是否找到xml文件
    if [[ ${#xml_files[@]} -eq 0 ]]; then
        print_warning "在meta目录中未找到xml文件"
        exit 0
    fi
    
    print_info "找到 ${#xml_files[@]} 个xml文件，开始处理..."
    echo ""
    
    # 使用parallel进行并行处理
    local processed=0
    local failed=0
    
    printf '%s\n' "${xml_files[@]}" | parallel -j "${MAX_PARALLEL}" \
        "bash -c 'source \"$0\" && move_xml_file \"{}\"'" && {
            processed=$((${#xml_files[@]}))
        } || {
            # 计算失败的文件数
            for file in "${xml_files[@]}"; do
                if ! is_file_processed "$(basename "$file")"; then
                    failed=$((failed + 1))
                fi
            done
            processed=$((${#xml_files[@]} - failed))
        }
    
    echo ""
    echo -e "${BLUE}╔════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║${NC}  处理完成                              ${BLUE}║${NC}"
    echo -e "${BLUE}╚════════════════════════════════════════╝${NC}"
    
    print_info "已处理: $processed/${#xml_files[@]} 个文件"
    
    if [[ $failed -gt 0 ]]; then
        print_warning "失败: $failed 个文件"
    fi
    
    echo ""
    
    # 列出meta目录中保留的其他文件
    print_info "meta目录中保留的其他文件:"
    while IFS= read -r file; do
        echo "  - $(basename "$file")"
    done < <(find "${META_DIR}" -maxdepth 1 -type f ! -name "*.xml")
    
    # 列出xml子文件夹中的文件
    echo ""
    print_info "xml文件夹中的xml文件数量: $(find "${XML_SUBDIR}" -maxdepth 1 -name "*.xml" -type f | wc -l)"
    
    # 清理进度文件
    rm -f "${PROGRESS_FILE}"
    
    if [[ $failed -eq 0 ]]; then
        print_success "所有文件整理完成！"
        exit 0
    else
        print_error "部分文件整理失败"
        exit 1
    fi
}

# ============================================================================
# 脚本入口
# ============================================================================
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
