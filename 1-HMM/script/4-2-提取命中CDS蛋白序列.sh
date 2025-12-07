#!/usr/bin/env bash
################################################################################
# 根据 manifest 提取蛋白 + CDS 序列（顺序版）
# 变量需外部传入：SAMPLE_DIR, MANIFEST, OUTPUT_DIR, PYTHON_SCRIPT, QUERIES (可选)
################################################################################

set -Eeuo pipefail

require_var() {
  local name="$1"
  if [[ -z "${!name:-}" ]]; then
    echo "[ERROR] 变量未设置: $name"
    exit 1
  fi
}

require_var SAMPLE_DIR
require_var MANIFEST
require_var OUTPUT_DIR
require_var PYTHON_SCRIPT

read -r -a queries_array <<< "${QUERIES:-}"

cmd=(python3 "$PYTHON_SCRIPT" --manifest "$MANIFEST" --sample-dir "$SAMPLE_DIR" --output-dir "$OUTPUT_DIR")
if [[ ${#queries_array[@]} -gt 0 ]]; then
  cmd+=(--queries "${queries_array[@]}")
fi

echo "运行命令: ${cmd[*]}"
"${cmd[@]}"
