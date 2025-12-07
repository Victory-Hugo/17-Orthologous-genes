#!/usr/bin/env bash
#*=========================
#*======建立 HMM 库=======
#*=========================
set -Eeuo pipefail
shopt -s nullglob

INPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/data/1-HMM序列"
OUTPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/data/2-HMM库"
FILE_GLOB_PATTERN="*.aln.faa"         # 匹配输入 MSA 的文件名模式
HMMPROG="hmmbuild"                     # hmmbuild 程序名（在 PATH 中）


need_bin() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "[ERR ] 缺少可执行文件：$1 ；请安装后重试。"
    exit 127
  fi
}
need_bin "$HMMPROG"

mkdir -p "$OUTPUT_DIR"

############################
# 执行 hmmbuild             
############################
files=("$INPUT_DIR"/$FILE_GLOB_PATTERN)

if [[ ${#files[@]} -eq 0 ]]; then
  echo "[INFO] 输入目录中未找到匹配的文件：$FILE_GLOB_PATTERN"
  exit 0
fi

total=${#files[@]}
index=0

for in_file in "${files[@]}"; do
  index=$((index + 1))
  base="$(basename "$in_file")"
  stem="${base%.aln.faa}"
  out_file="${OUTPUT_DIR}/${stem}.hmm"

  echo "[INFO] (${index}/${total}) 开始构建：$in_file -> $out_file"
  "$HMMPROG" "$out_file" "$in_file"
echo "[DONE] 完成：$out_file"
done

echo "[ALL DONE] 所有任务顺序完成。"
