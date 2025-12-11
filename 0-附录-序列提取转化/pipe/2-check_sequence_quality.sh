#!/bin/bash

# 序列质量检查工具 - Shell 调用脚本
# 用法: ./check_sequence_quality.sh <输入FASTA文件> <输出报告路径>
# 示例: ./check_sequence_quality.sh input.fna output/report.json

set -e

# 获取脚本所在目录
SCRIPT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/0-附录-序列提取转化/script"
PYTHON_SCRIPT="${SCRIPT_DIR}/check_sequence_quality.py"

# 参数检查
if [ $# -lt 2 ]; then
    echo "用法: $(basename $0) <输入FASTA文件> <输出报告路径>"
    echo ""
    echo "参数说明:"
    echo "  <输入FASTA文件>   - FASTA 格式的核酸序列文件"
    echo "  <输出报告路径>    - 输出报告文件路径 (建议使用 .json 扩展名)"
    echo ""
    echo "示例:"
    echo "  $(basename $0) sequences.fna reports/quality_report.json"
    echo "  $(basename $0) merge_pick.fna ./output/report.json"
    exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="$2"

# 检查输入文件是否存在
if [ ! -f "$INPUT_FILE" ]; then
    echo "错误: 输入文件不存在 - $INPUT_FILE"
    exit 1
fi

# 检查 Python 脚本是否存在
if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "错误: Python 脚本不存在 - $PYTHON_SCRIPT"
    exit 1
fi

# 检查 Python 是否可用
if ! command -v python3 &> /dev/null; then
    echo "错误: Python3 未安装或不在 PATH 中"
    exit 1
fi

# 创建输出目录
OUTPUT_DIR="$(dirname "$OUTPUT_FILE")"
if [ ! -z "$OUTPUT_DIR" ] && [ "$OUTPUT_DIR" != "." ]; then
    mkdir -p "$OUTPUT_DIR"
fi

# 运行 Python 脚本
echo "开始序列质量检查..."
python3 "$PYTHON_SCRIPT" "$INPUT_FILE" "$OUTPUT_FILE"

exit_code=$?

if [ $exit_code -eq 0 ]; then
    echo ""
    echo "✓ 质量检查完成！"
    echo "✓ 详细报告位置:"
    echo "  - JSON 报告: $OUTPUT_FILE"
    echo "  - 文本报告: ${OUTPUT_FILE%.json}.txt"
else
    echo ""
    echo "✗ 检查过程出现错误，请检查输入文件格式"
    exit $exit_code
fi
