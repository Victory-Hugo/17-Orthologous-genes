#! /bin/bash

INPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-数据获取/2-下载特定物种序列/download"
OUTPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-数据获取/2-下载特定物种序列/conf/"

cd "${INPUT_DIR}" || { echo "无法进入目录: ${INPUT_DIR}"; exit 1; }

cat ${INPUT_DIR}/*complete_genomes.txt | sort | uniq >> "${OUTPUT_DIR}/balanced_accessions_ID.txt"