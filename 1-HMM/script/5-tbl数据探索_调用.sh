#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/2-HMM/script/5-tbl数据探索.py"
R_SCRIPT="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/2-HMM/script/5-tbl数据探索.R"

INPUT=""
OUTPUT=""

usage() {
  cat <<'EOF'
Usage: 5-tbl数据探索_调用.sh --input <tbl.tsv> --output <dir>

Runs the Python analysis first and then triggers the R-based visualization.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)
      INPUT="${2:-}"
      shift 2
      ;;
    -o|--output)
      OUTPUT="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "[ERROR] Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ -z "${INPUT}" || -z "${OUTPUT}" ]]; then
  echo "[ERROR] Missing required arguments." >&2
  usage
  exit 1
fi

if [[ ! -f "${INPUT}" ]]; then
  echo "[ERROR] Input file not found: ${INPUT}" >&2
  exit 1
fi

mkdir -p "${OUTPUT}"

echo "[INFO] Running Python data analysis..."
python3 "${PY_SCRIPT}" --input "${INPUT}" --output "${OUTPUT}"

echo "[INFO] Running R visualization..."
Rscript "${R_SCRIPT}" --input "${INPUT}" --output "${OUTPUT}"

echo "[DONE] All tasks finished."

