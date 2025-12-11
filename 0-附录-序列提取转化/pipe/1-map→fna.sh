#!/usr/bin/env bash
set -euo pipefail

#*==============================
#* 脚本功能：根据映射文件提取指定基因的核酸序列，生成fna文件
#*==============================
base_dir="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/0-附录-序列提取转化"
core_script="${base_dir}/script/map_to_fna_core.sh"
threads="16"
#*==============================
map_file="${base_dir}/conf/map.csv"                                 #? 输入文件
out_dir="/mnt/l/18-Rv0194-Gene/3-SingleGene-tree/data/Rv0194"       #? 输出目录


if [[ ! -x "$core_script" ]]; then
    echo "未找到可执行核心脚本: $core_script" >&2
    exit 1
fi

exec "$core_script" \
    "$map_file" \
    "$out_dir" \
    "$threads"
