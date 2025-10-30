#!/bin/bash

# 脚本功能：从 assembly_biosample_map.csv 文件提取第1、2、4列
# 输入文件：/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-数据获取/2-下载特定物种序列/meta/assembly_biosample_map.csv
# 输出文件：/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-数据获取/2-下载特定物种序列/meta/species_name.txt

INPUT_FILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-数据获取/2-下载特定物种序列/meta/assembly_biosample_map.csv"
OUTPUT_FILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-数据获取/2-下载特定物种序列/meta/species_name.txt"

# 提取第1、2、4列（排除标题行）
tail -n +2 "$INPUT_FILE" | cut -d',' -f1,2,4 > "$OUTPUT_FILE"

echo "物种名称提取完成！"
echo "输入文件：$INPUT_FILE"
echo "输出文件：$OUTPUT_FILE"
echo "提取的物种数量（去重后）：$(wc -l < "$OUTPUT_FILE")"
