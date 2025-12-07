#!/usr/bin/env bash
################################################################################
# HMM 蛋白序列比对（简化顺序版）
# 作者: BigLin
# 功能: 按列表逐一执行 hmmsearch，无彩色输出、无断点续跑，最简单的双层 for 循环
################################################################################

set -Eeuo pipefail
shopt -s nullglob

################################
# 配置（由外部传入）
################################
# 必传：LIST_FILE, HMM_FILE_TXT, OUTPUT_DIR, CPU_PER_JOB, E_VALUE, MATRIX_SCRIPT,
#      PRESENCE_MATRIX_TSV, MERGED_TBL_TSV, MATRIX_EVALUE_CUTOFF
require_var() {
  local name="$1"
  if [[ -z "${!name:-}" ]]; then
    echo "[ERROR] 变量未设置: $name"
    exit 1
  fi
}

require_var LIST_FILE
require_var HMM_FILE_TXT
require_var OUTPUT_DIR
require_var CPU_PER_JOB
require_var E_VALUE
require_var MATRIX_SCRIPT
require_var PRESENCE_MATRIX_TSV
require_var MERGED_TBL_TSV
require_var MATRIX_EVALUE_CUTOFF

################################
# 基础函数
################################
need_bin() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "[ERROR] 未找到命令: $1"
    exit 127
  fi
}

trim_line() {
  # 去首尾空白
  printf '%s' "$1" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//'
}

################################
# 准备
################################
need_bin hmmsearch
mkdir -p "$OUTPUT_DIR"

echo "[INFO] 读取 FAA 列表: $LIST_FILE"
FAA_FILES=()
while IFS= read -r line || [[ -n "$line" ]]; do
  line="${line%$'\r'}"
  line="$(trim_line "$line")"
  [[ -z "$line" || "$line" == \#* ]] && continue
  FAA_FILES+=("$line")
done < "$LIST_FILE"

echo "[INFO] 读取 HMM 列表: $HMM_FILE_TXT"
HMM_FILES=()
while IFS= read -r line || [[ -n "$line" ]]; do
  line="${line%$'\r'}"
  line="$(trim_line "$line")"
  [[ -z "$line" || "$line" == \#* ]] && continue
  HMM_FILES+=("$line")
done < "$HMM_FILE_TXT"

TOTAL_FAA=${#FAA_FILES[@]}
TOTAL_HMM=${#HMM_FILES[@]}
TOTAL_TASKS=$((TOTAL_FAA * TOTAL_HMM))

echo "[INFO] FAA 文件数: $TOTAL_FAA"
echo "[INFO] HMM 文件数: $TOTAL_HMM"
echo "[INFO] 总任务数: $TOTAL_TASKS"

################################
# 主循环：顺序 hmmsearch
################################
task_idx=0
for hmm in "${HMM_FILES[@]}"; do
  hmm_name="$(basename "$hmm")"
  hmm_base="${hmm_name%.*}"
  out_dir="${OUTPUT_DIR}/${hmm_base}"
  mkdir -p "$out_dir"

  for faa in "${FAA_FILES[@]}"; do
    faa_name="$(basename "$faa")"
    faa_base="${faa_name%.*}"
    tbl="${out_dir}/${faa_base}.tbl"
    err="${out_dir}/${faa_base}.err"

    task_idx=$((task_idx + 1))
    echo "[TASK ${task_idx}/${TOTAL_TASKS}] $hmm_name vs $faa_name"

    if hmmsearch --cpu "$CPU_PER_JOB" -E "$E_VALUE" --tblout "$tbl" "$hmm" "$faa" 2> "$err"; then
      echo "[OK] 生成: $tbl"
      rm -f "$err"
    else
      echo "[FAIL] hmmsearch 失败: $hmm_name vs $faa_name；详情见 $err"
    fi
  done
done

################################
# 生成存在/缺失矩阵（如脚本存在）
################################
if [[ -f "$MATRIX_SCRIPT" ]]; then
  echo "[INFO] 生成存在/缺失矩阵..."
  if python3 "$MATRIX_SCRIPT" \
    --base-dir "$OUTPUT_DIR" \
    --output-tsv "$PRESENCE_MATRIX_TSV" \
    --merged-tsv "$MERGED_TBL_TSV" \
    --evalue-cutoff "$MATRIX_EVALUE_CUTOFF"; then
    echo "[OK] 矩阵: $PRESENCE_MATRIX_TSV"
    echo "[OK] 合并表: $MERGED_TBL_TSV"
  else
    echo "[WARN] 矩阵生成失败，请检查。"
  fi
else
  echo "[WARN] 未找到矩阵脚本: $MATRIX_SCRIPT"
fi

echo "[ALL DONE] 全部 hmmsearch 顺序完成。"
