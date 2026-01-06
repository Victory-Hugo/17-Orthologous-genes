#!/usr/bin/env bash
# 依赖列表:
# - clinker (CLI)
# - clustermap.js (可选)
# - python3
# - biopython
# - gffutils (可选)

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
CONF_FILE="/mnt/f/onedrive/文档（科研）/脚本/Download/17-Orthologous-genes/8-clinker/2-运行/conf/clinker.conf"

if [[ ! -f "${CONF_FILE}" ]]; then
  echo "配置文件不存在: ${CONF_FILE}" >&2
  exit 1
fi

# shellcheck disable=SC1090
source "${CONF_FILE}"

PY_SCRIPT="${SCRIPT_DIR}/run_clinker.py"
if [[ ! -f "${PY_SCRIPT}" ]]; then
  echo "脚本不存在: ${PY_SCRIPT}" >&2
  exit 1
fi

ARGS=(
  --conf "${CONF_TSV}"
  --outdir "${OUTDIR}"
  --logdir "${LOGDIR}"
  --flank "${FLANK}"
  --threads "${THREADS}"
)

if [[ "${KEEP_TMP}" == "1" ]]; then
  ARGS+=(--keep_tmp)
fi

python3 "${PY_SCRIPT}" "${ARGS[@]}" "$@"
