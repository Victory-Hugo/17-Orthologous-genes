#!/usr/bin/env bash

set -euo pipefail
#*=======注意激活conda========
# conda activate gtotree
#*=======注意激活conda========

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONF_FILE="/mnt/f/onedrive/文档（科研）/脚本/Download/17-Orthologous-genes/2-GToTree/example/1-多物种序列→直接建树/conf/多物种序列建树.conf"
if [[ $# -ge 1 ]]; then
  CONF_FILE="$1"
fi

if [[ ! -f "${CONF_FILE}" ]]; then
  printf "错误：找不到配置文件：%s\n" "${CONF_FILE}" >&2
  exit 1
fi

# shellcheck disable=SC1090
source "${CONF_FILE}"

required_vars=(
  SCG_SET
  PARALLEL_JOBS
  HMM_CPUS
  MUSCLE_THREADS
  PROJECT_ROOT
  DATA_DIR
  OUTPUT_DIR
  LOG_DIR
  CONF_DIR
  FASTA_LIST
  TREE_METHOD
  KEEP_ALIGN
  RUN_ID_PREFIX
  RUN_ID_TIME_FORMAT
  LOG_TIME_FORMAT
  LOG_FILE
)

for var_name in "${required_vars[@]}"; do
  if [[ -z "${!var_name:-}" ]]; then
    printf "错误：配置文件缺少变量：%s\n" "${var_name}" >&2
    exit 1
  fi
done

RUN_ID="${RUN_ID_PREFIX}_$(date +"${RUN_ID_TIME_FORMAT}")"
RUN_DIR="${OUTPUT_DIR}/${RUN_ID}"

# 若时间戳重复，自动等待并生成新的运行目录
while [[ -d "${RUN_DIR}" ]]; do
  sleep 1
  RUN_ID="${RUN_ID_PREFIX}_$(date +"${RUN_ID_TIME_FORMAT}")"
  RUN_DIR="${OUTPUT_DIR}/${RUN_ID}"
done


mkdir -p "${DATA_DIR}" "${OUTPUT_DIR}" "${LOG_DIR}"

timestamp() {
  date +"${LOG_TIME_FORMAT}"
}

log() {
  printf "[%s] %s\n" "$(timestamp)" "$*" | tee -a "${LOG_FILE}"
}

log "GToTree 开始运行，使用基因集：${SCG_SET}，输出目录：${RUN_DIR}"
#? -k：保留每个目标基因的独立比对（默认：false）。
#? -T <字符串>：指定构树程序，可选 "FastTree"（默认）或 "IQ-TREE"。

gtotree_args=(
  -f "${FASTA_LIST}"
  -H "${SCG_SET}"
  -o "${RUN_DIR}"
  -j "${PARALLEL_JOBS}"
  -n "${HMM_CPUS}"
  -M "${MUSCLE_THREADS}"
  -T "${TREE_METHOD}"
)

if [[ "${KEEP_ALIGN}" == "true" ]]; then
  gtotree_args+=(-k)
fi

if [[ -n "${MAP_FILE:-}" ]]; then
  if [[ ! -f "${MAP_FILE}" ]]; then
    log "错误：映射文件不存在：${MAP_FILE}"
    exit 1
  fi
  gtotree_args+=(-m "${MAP_FILE}")
fi

GToTree "${gtotree_args[@]}"
log "GToTree 运行完成，结果位于：${RUN_DIR}"
