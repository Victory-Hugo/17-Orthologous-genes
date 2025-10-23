#!/bin/bash
set -euo pipefail

INDIR="/mnt/f/11_钟杰_cydCD/1-数据下载/download/cydD_蛋白"
OUT="/mnt/f/11_钟杰_cydCD/1-数据下载/download/cydD_merge.faa"
ID_TXT="/mnt/f/11_钟杰_cydCD/1-数据下载/download/cydD_merge.ID.txt"
mkdir -p "$(dirname "$OUT")"
: > "$OUT"

# 每个 .faa 文件抽一条：NP_ > YP_ > WP_ > 其他(则取第一条)，并在头行前加 GeneID
find "$INDIR" -type f -name 'Gene_ID_*.faa' -print0 | while IFS= read -r -d '' f; do
  gid="$(basename "$f" | sed -n 's/.*Gene_ID_\([0-9][0-9]*\)\.faa.*/\1/p')"
  [[ -z "$gid" ]] && { echo "跳过：无法解析GeneID -> $f" >&2; continue; }

  awk -v GID="$gid" '
    function pick() {
      if (cur_h!="") {
        p = (cur_h ~ />.*NP_/)?3:((cur_h ~ />.*YP_/)?2:((cur_h ~ />.*WP_/)?1:0))
        if (!seen || p > bestp) { besth=cur_h; bests=cur_s; bestp=p; seen=1 }
      }
      cur_h=""; cur_s=""
    }
    BEGIN { bestp=-1; seen=0; cur_h=""; cur_s="" }
    {
      sub(/\r$/,"")                                  # 去掉CR
      if ($0 ~ /^>/) { pick(); cur_h=$0; cur_s="" }  # 新条目
      else           { cur_s=cur_s $0 "\n" }         # 序列
    }
    END {
      pick();
      if (seen) {
        gsub(/^>/, ">" GID "|", besth);
        print besth; printf "%s", bests;
      }
    }
  ' "$f" >> "$OUT"
done

grep '>' ${OUT} \
    |awk -v FS='|' -v OFS='\t' '{print $1,$2}'\
    |awk -v FS='[' -v OFS='\t' '{print $1,$2}'\
    |awk -v FS=']' '{print $1,$2}'\
    |awk -v FS='>' '{print $2}' >  ${ID_TXT}