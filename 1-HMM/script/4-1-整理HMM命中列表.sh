#!/usr/bin/env bash
################################################################################
# 整理 *.hits.faa 为 TSV 清单（顺序版）
# 变量需外部传入：HITS_DIR, OUTPUT_TSV, PYTHON_SCRIPT, QUERIES (可选，空表示全部)
################################################################################

set -Eeuo pipefail

require_var() {
  local name="$1"
  if [[ -z "${!name:-}" ]]; then
    echo "[ERROR] 变量未设置: $name"
    exit 1
  fi
}

require_var HITS_DIR
require_var OUTPUT_TSV
require_var PYTHON_SCRIPT

read -r -a queries_array <<< "${QUERIES:-}"

cmd=(python3 "$PYTHON_SCRIPT" --hits-dir "$HITS_DIR" --output-tsv "$OUTPUT_TSV")
if [[ ${#queries_array[@]} -gt 0 ]]; then
  cmd+=(--queries "${queries_array[@]}")
fi

echo "运行命令: ${cmd[*]}"
"${cmd[@]}"
