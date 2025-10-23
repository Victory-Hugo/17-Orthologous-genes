#!/bin/bash

# 设置变量
LIST_FILE="/mnt/c/Users/Administrator/Desktop/Lol-家族系统发育树/conf/鲍曼faa.list.txt"
RESULTS_DIR="/mnt/c/Users/Administrator/Desktop/Lol-家族系统发育树/data/鲍曼HMM结果/"
OUTPUT_DIR="/mnt/c/Users/Administrator/Desktop/Lol-家族系统发育树/data/鲍曼HMM结果/"
SCRIPT_DIR="/mnt/c/Users/Administrator/Desktop/Lol-家族系统发育树/script"
JOBS=4  # 并行任务数

# 确保输出目录存在
mkdir -p "$OUTPUT_DIR"

echo "=" * 80
echo "开始使用 parallel 提取最佳匹配序列"
echo "=" * 80
echo "并行任务数: $JOBS"
echo ""

# 使用 parallel 并行处理
cat "$LIST_FILE" | parallel --jobs $JOBS --tag \
    python3 "$SCRIPT_DIR/4-提取最佳匹配序列.py" {} "$RESULTS_DIR" "$OUTPUT_DIR"

echo ""
echo "=" * 80
echo "所有文件处理完成！"
echo "输出目录: $OUTPUT_DIR"
echo "=" * 80
