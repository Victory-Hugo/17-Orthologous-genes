#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG_DIR="${BASE_DIR}/log"
TMP_DIR="${BASE_DIR}/tmp"
PYTHON_BIN="${PYTHON_BIN:-python3}"
CONFIG_PATH="${BASE_DIR}/conf/hyphy_visualization.json"
JOBS="${JOBS:-1}"
FORCE="false"

usage() {
  cat <<USAGE
Usage: $0 [--config PATH] [--jobs N] [--force]
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG_PATH="$2"
      shift 2
      ;;
    --jobs)
      JOBS="$2"
      shift 2
      ;;
    --force)
      FORCE="true"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" 1>&2
      usage
      exit 1
      ;;
  esac
done

mkdir -p "$LOG_DIR" "$TMP_DIR"
SCRIPT_LOG="${LOG_DIR}/$(basename "$0").log"

exec > >(tee -a "$SCRIPT_LOG") 2>&1

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

trap 'status=$?; if [[ $status -ne 0 ]]; then log "STATUS=FAIL exit_code=$status"; else log "STATUS=SUCCESS"; fi' EXIT

log "START config=${CONFIG_PATH} jobs=${JOBS} force=${FORCE}"

if [[ ! -f "$CONFIG_PATH" ]]; then
  log "ERROR config not found: $CONFIG_PATH"
  exit 1
fi

DATASET_CMD=("$PYTHON_BIN" "${BASE_DIR}/python/hyphy_visualization.py" --config "$CONFIG_PATH" --list-datasets)
RUN_CMD=("$PYTHON_BIN" "${BASE_DIR}/python/hyphy_visualization.py" --config "$CONFIG_PATH" --log-dir "$LOG_DIR" --tmp-dir "$TMP_DIR")
if [[ "$FORCE" == "true" ]]; then
  RUN_CMD+=(--force)
fi
RUN_CMD+=(--dataset)

if command -v parallel >/dev/null 2>&1; then
  "${DATASET_CMD[@]}" | parallel -j "$JOBS" --halt now,fail=1 "${RUN_CMD[@]}" {}
else
  "${DATASET_CMD[@]}" | xargs -r -n 1 -P "$JOBS" -I {} "${RUN_CMD[@]}" {}
fi

log "DONE"
