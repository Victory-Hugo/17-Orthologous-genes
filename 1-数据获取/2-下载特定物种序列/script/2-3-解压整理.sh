#!/usr/bin/env bash
set -euo pipefail
#! 使用之前请先自己解压到单独的文件夹
#*==========示例如下==========
# ├── GCF_003063885.1_downloaded
# │   ├── README.md
# │   ├── md5sum.txt
# │   └── ncbi_dataset
# │       └── data
# │           ├── GCF_003063885.1
# │           │   ├── GCF_003063885.1_ASM306388v1_genomic.fna
# │           │   ├── cds_from_genomic.fna
# │           │   └── genomic.gff
# │           ├── assembly_data_report.jsonl
# │           └── dataset_catalog.json
#*============================
ROOT="/mnt/f/15_Bam_Tam/5-补齐更多物种/download"

############################
# 1. 重命名 *_downloaded 目录
############################
for dir in "$ROOT"/*_downloaded; do
  [ -d "$dir" ] || continue
  new_name="${dir/_downloaded/}"
  mv "$dir" "$new_name"
done

############################
# 2. 删除无用文件 (并行)
############################
cd "$ROOT"
for pat in "README.md" "md5sum.txt" "assembly_data_report.jsonl" "dataset_catalog.json"; do
  find "$ROOT" -type f -name "$pat" | parallel --bar -j 8 rm {}
done


##################################################
# 3. 统一重命名工具函数
##################################################
rename_with_policy() {
  local src="$1" dst="$2"
  if [ -e "$dst" ]; then
    if [ "$OVERWRITE" -eq 1 ]; then
      echo "OVERWRITE: '$src' -> '$dst'"
      mv -f -- "$src" "$dst"
    elif [ "$BACKUP" -eq 1 ]; then
      bak="${dst}.bak.$(date +%s)"
      echo "BACKUP: '$dst' -> '$bak'"
      mv -- "$dst" "$bak"
      mv -v -- "$src" "$dst"
    else
      echo "SKIP (目标已存在): '$src' -> '$dst'"
    fi
  else
    mv -v -- "$src" "$dst"
  fi
}
##################################################


############################
# 4. 重命名 cds_from_genomic.fna -> <parent>_CDS.fna
############################
OVERWRITE=0; BACKUP=0; ROOT_SET=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --overwrite) OVERWRITE=1; shift ;;
    --backup) BACKUP=1; shift ;;
    --help|-h)
      echo "Usage: $0 [--overwrite] [--backup] [ROOT_DIRECTORY]"
      exit 0 ;;
    *) if [ $ROOT_SET -eq 0 ]; then ROOT="$1"; ROOT_SET=1; fi; shift ;;
  esac
done

find "$ROOT" -type f -name 'cds_from_genomic.fna' -print0 |
while IFS= read -r -d '' file; do
  dir="$(dirname "$file")"
  parent="$(basename "$dir")"
  newpath="$dir/${parent}_CDS.fna"
  rename_with_policy "$file" "$newpath"
done


############################
# 5. 重命名 genomic.gff -> <parent>.gff
############################
find "$ROOT" -type f -name 'genomic.gff' -print0 |
while IFS= read -r -d '' file; do
  dir="$(dirname "$file")"
  parent="$(basename "$dir")"
  newpath="$dir/${parent}.gff"
  rename_with_policy "$file" "$newpath"
done


############################
# 6. 重命名 *.fna -> <parent>.fasta (排除 *_CDS.fna)
############################
find "$ROOT" -type f -name '*.fna' ! -name '*_CDS.fna' -print0 |
while IFS= read -r -d '' file; do
  dir="$(dirname "$file")"
  parent="$(basename "$dir")"
  newpath="$dir/${parent}.fasta"
  rename_with_policy "$file" "$newpath"
done


############################
# 7. 重命名 protein.faa -> <parent>.faa
############################
find "$ROOT" -type f -name 'protein.faa' -print0 |
while IFS= read -r -d '' file; do
  dir="$(dirname "$file")"
  parent="$(basename "$dir")"
  newpath="$dir/${parent}.faa"
  rename_with_policy "$file" "$newpath"
done


############################
# 8. 展开 ncbi_dataset/data 目录结构
############################
for sample in "$ROOT"/*; do
  [ -d "$sample" ] || continue
  data_base="$sample/ncbi_dataset/data"
  [ -d "$data_base" ] || continue

  for inner in "$data_base"/*; do
    [ -d "$inner" ] || continue
    echo "处理 inner='$inner'"
    find "$inner" -mindepth 1 -maxdepth 1 -print0 |
    while IFS= read -r -d '' item; do
      name="$(basename "$item")"
      dest="$sample/$name"
      rename_with_policy "$item" "$dest"
    done
    rmdir "$inner" 2>/dev/null || echo "保留非空目录: $inner"
  done

  rmdir "$data_base" 2>/dev/null || echo "保留非空目录: $data_base"
  rmdir "$sample/ncbi_dataset" 2>/dev/null || echo "保留非空目录: $sample/ncbi_dataset"
done

echo "全部完成。"
