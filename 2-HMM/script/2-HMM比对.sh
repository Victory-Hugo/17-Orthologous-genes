#!/bin/bash

# 设置变量
LIST_FILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/2-HMM/conf/鲍曼faa.list.txt"
HMM_FILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/2-HMM/data/LolD.hmm"
OUTPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/2-HMM/data/鲍曼HMM结果/"
JOBS=4  # 并行任务数

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 定义处理单个文件的函数
process_faa_file() {
    local faa_file="$1"
    local output_dir="$2"
    local hmm_file="$3"
    
    # 跳过空行
    [[ -z "$faa_file" ]] && return
    
    # 检查文件是否存在
    if [[ ! -f "$faa_file" ]]; then
        echo "错误: 文件不存在 - $faa_file"
        return 1
    fi
    
    # 提取文件名（不含路径和扩展名）
    filename=$(basename "$faa_file" .faa)
    
    # 设置输出文件名
    output_file="${output_dir}/hits_${filename}.tbl"
    
    echo "正在处理: $faa_file"
    
    # 执行 hmmsearch
    hmmsearch --cpu 2 \
        -E 1e-50 \
        --tblout "$output_file" \
        "$hmm_file" \
        "$faa_file"
    
    if [[ $? -eq 0 ]]; then
        echo "✓ 完成: $filename"
    else
        echo "✗ 失败: $filename"
        return 1
    fi
}

export -f process_faa_file

# 使用 parallel 并行处理
echo "开始使用 parallel 进行 HMM 比对..."
cat "$LIST_FILE" | parallel --jobs $JOBS --tag process_faa_file {} "$OUTPUT_DIR" "$HMM_FILE"

echo ""
echo "所有文件处理完成！" 