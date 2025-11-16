#!/usr/bin/env bash

set -euo pipefail
#*=======注意激活conda========
# conda activate gtotree
#*=======注意激活conda========


SCG_SET="Bacteria"    #? 设定这些物种所在的域
PARALLEL_JOBS="16"    #? GToTree 的 -j 选项
HMM_CPUS="16"         #? GToTree 的 -n 选项
MUSCLE_THREADS="16"   #? GToTree 的 -M 选项

#* -------- 目录设定 --------
PROJECT_ROOT="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/2-GToTree/example/1-多物种序列→直接建树"
DATA_DIR="${PROJECT_ROOT}/data"       #? 存放输入数据的目录,默认是每个物种一个fasta文件
OUTPUT_DIR="${PROJECT_ROOT}/output"   #? 存放输出结果的目录,自动生成
LOG_DIR="${OUTPUT_DIR}/logs"          #? 存放日志文件的目录,自动生成
#* -------- 目录设定 --------
CONF_DIR="${PROJECT_ROOT}/conf"       #? 存放配置文件的目录，其中给出fasta文件列表和基因组ID映射文件
FASTA_LIST="${CONF_DIR}/fasta_files.txt" #? 包含所有fasta文件路径的文本文件
MAP_FILE="${CONF_DIR}/genome_to_id_map.tsv" #? 基因组文件名到物种ID的映射文件
RUN_ID="run_$(date +%Y%m%d_%H%M%S)"
RUN_DIR="${OUTPUT_DIR}/${RUN_ID}"
LOG_FILE="${LOG_DIR}/GToTree.log"

# 若时间戳重复，自动等待并生成新的运行目录
while [[ -d "${RUN_DIR}" ]]; do
  sleep 1
  RUN_ID="run_$(date +%Y%m%d_%H%M%S)"
  RUN_DIR="${OUTPUT_DIR}/${RUN_ID}"
done


mkdir -p "${DATA_DIR}" "${OUTPUT_DIR}" "${LOG_DIR}"

timestamp() {
  date +"%Y-%m-%d %H:%M:%S"
}

log() {
  printf "[%s] %s\n" "$(timestamp)" "$*" | tee -a "${LOG_FILE}"
}

log "GToTree 开始运行，使用基因集：${SCG_SET}，输出目录：${RUN_DIR}"
#? -k：保留每个目标基因的独立比对（默认：false）。
#? -T <字符串>：指定构树程序，可选 "FastTree"（默认）或 "IQ-TREE"。
GToTree \
  -f "${FASTA_LIST}" \
  -H "${SCG_SET}" \
  -m "${MAP_FILE}" \
  -o "${RUN_DIR}" \
  -j "${PARALLEL_JOBS}" \
  -n "${HMM_CPUS}" \
  -M "${MUSCLE_THREADS}" \
  -k \
  -T "IQ-TREE"

log "GToTree 运行完成，结果位于：${RUN_DIR}"
