#!/bin/bash
GENE_NAME="lolD"
INPUT_TSV="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-数据获取/conf/gene_result.txt"
OUTPUT_TXT="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-数据获取/conf/${GENE_NAME}.txt"



csvgrep -t -c 6 \
    -m "$GENE_NAME" \
    "$INPUT_TSV" |\
    csvcut -c 3 -K 1  \
    > "$OUTPUT_TXT"