#!/usr/bin/env bash

# ------------------------------------------------------------
# 使用说明：
#   1. 默认读取 conf/1-fna_path.txt 列出的 FASTA 文件。
#   2. 可通过环境变量 THREADS/OUTPUT_DIR/PRODIGAL_MODE 覆盖默认设置。
#   3. 支持断点续跑，会自动跳过已经生成完 gff/faa/ffn 的样本。
#   4. 借助 GNU parallel --bar 输出进度条，并将日志写入 output/logs。
# ------------------------------------------------------------

set -euo pipefail

BASE_DIR="/mnt/c/Users/Administrator/Desktop"
DEFAULT_LIST="$BASE_DIR/conf/1-fna_path.txt"
DEFAULT_OUTPUT="$BASE_DIR/output"

# 根据硬件情况推断默认线程数，优先使用 nproc。
if command -v nproc >/dev/null 2>&1; then
  DEFAULT_THREADS="$(nproc)"
else
  DEFAULT_THREADS=1
fi

OUTPUT_DIR="${OUTPUT_DIR:-$DEFAULT_OUTPUT}"
PRODIGAL_MODE="${PRODIGAL_MODE:-single}"
THREADS="${THREADS:-$DEFAULT_THREADS}"
LOG_DIR="$OUTPUT_DIR/logs"

# 清理临时文件用的 trap。
PENDING_LIST=""
cleanup() {
  if [[ -n "$PENDING_LIST" && -f "$PENDING_LIST" ]]; then
    rm -f "$PENDING_LIST"
  fi
}
trap cleanup EXIT

# 提取样本名（去掉扩展名）。
sample_name() {
  local fasta_path="$1"
  local file_name
  file_name="$(basename -- "$fasta_path")"
  echo "${file_name%.*}"
}

# 工作进程：真正调用 prodigal 的地方。
run_worker() {
  local fasta_path="${1:-}"

  if [[ -z "$fasta_path" ]]; then
    echo "[ERROR] Worker received an empty FASTA path." >&2
    exit 1
  fi

  if [[ ! -s "$fasta_path" ]]; then
    echo "[ERROR] Input FASTA not found or empty: $fasta_path" >&2
    exit 1
  fi

  mkdir -p "$OUTPUT_DIR"

  local sample
  sample="$(sample_name "$fasta_path")"

  local gff_file="$OUTPUT_DIR/${sample}.gff"
  local protein_file="$OUTPUT_DIR/${sample}.faa"
  local nucleotide_file="$OUTPUT_DIR/${sample}.ffn"

  prodigal \
    -i "$fasta_path" \
    -o "$gff_file" \
    -a "$protein_file" \
    -d "$nucleotide_file" \
    -f gff \
    -p "$PRODIGAL_MODE" \
    -q
}

# 被 GNU parallel 以 --worker 方式调用。
if [[ "${1:-}" == "--worker" ]]; then
  shift
  run_worker "$@"
  exit 0
fi

INPUT_LIST="${1:-$DEFAULT_LIST}"

if [[ ! -f "$INPUT_LIST" ]]; then
  echo "[ERROR] Input list not found: $INPUT_LIST" >&2
  exit 1
fi

command -v prodigal >/dev/null 2>&1 || {
  echo "[ERROR] prodigal is not available in PATH." >&2
  exit 1
}

command -v parallel >/dev/null 2>&1 || {
  echo "[ERROR] GNU parallel is required but not found." >&2
  exit 1
}

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# 解析输入列表，筛出需要继续运行的 FASTA。
PENDING_LIST="$(mktemp)"
total_count=0
skipped_count=0
while IFS= read -r line || [[ -n "$line" ]]; do
  line="${line%$'\r'}"   # 兼容 Windows 行尾
  [[ -z "$line" ]] && continue
  ((++total_count))

  sample="$(sample_name "$line")"
  gff="$OUTPUT_DIR/${sample}.gff"
  faa="$OUTPUT_DIR/${sample}.faa"
  ffn="$OUTPUT_DIR/${sample}.ffn"

  if [[ -s "$gff" && -s "$faa" && -s "$ffn" ]]; then
    ((++skipped_count))
    continue
  fi

  printf '%s\n' "$line" >>"$PENDING_LIST"
done <"$INPUT_LIST"

pending_count=$(wc -l <"$PENDING_LIST" | tr -d '[:space:]')
echo "[INFO] 共 ${total_count} 个序列，已完成 ${skipped_count} 个，待处理 ${pending_count} 个。"

if [[ "$pending_count" -eq 0 ]]; then
  echo "[INFO] 所有样本均已完成，退出。"
  exit 0
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SCRIPT_PATH="$SCRIPT_DIR/$(basename "$0")"
LOG_FILE="$LOG_DIR/prodigal_parallel.log"

parallel \
  --will-cite \
  --bar \
  --jobs "$THREADS" \
  --joblog "$LOG_FILE" \
  --halt soon,fail=1 \
  "$SCRIPT_PATH" --worker {} \
  :::: "$PENDING_LIST"
