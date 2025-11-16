#!/usr/bin/env bash

set -euo pipefail

#*=======注意激活conda========
# conda activate gtotree
#*=======注意激活conda========

# -------- 运行参数，可按需通过环境变量覆写 --------
PARALLEL_JOBS="16"
HMM_CPUS="16"
MUSCLE_THREADS="16"
EXPECTED_CONDA_ENV="gtotree"
BEST_HIT_MODE="true"

# -------- 目录设定 --------
PROJECT_ROOT="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/2-GToTree/example/2-目标HMM→命中序列"
CONF_DIR="${PROJECT_ROOT}/conf"
OUTPUT_ROOT="${PROJECT_ROOT}/output"
LOG_DIR="${OUTPUT_ROOT}/logs"

HMM_LIST="${CONF_DIR}/hmm_list.txt"
SEQUENCE_LIST="${CONF_DIR}/sequence_path.txt"
LOG_FILE="${LOG_DIR}/target_hmm_search.log"

mkdir -p "${OUTPUT_ROOT}" "${LOG_DIR}"
TMP_HMM_FILES=()

cleanup() {
  for tmp in "${TMP_HMM_FILES[@]}"; do
    [[ -f "${tmp}" ]] && rm -f "${tmp}"
  done
}
trap cleanup EXIT

timestamp() {
  date +"%Y-%m-%d %H:%M:%S"
}

log() {
  printf "[%s] %s\n" "$(timestamp)" "$*" | tee -a "${LOG_FILE}"
}

log "开始目标 HMM 搜索流程。"

# 检查 conda 环境与依赖
if [[ "${CONDA_DEFAULT_ENV:-}" != "${EXPECTED_CONDA_ENV}" ]]; then
  log "警告：当前 conda 环境为 \"${CONDA_DEFAULT_ENV:-未激活}\"，建议使用 \"${EXPECTED_CONDA_ENV}\"。"
fi

if ! command -v GToTree >/dev/null 2>&1; then
  log "错误：未找到 GToTree，可先激活 conda 环境 ${EXPECTED_CONDA_ENV}。"
  exit 1
fi

# 基础输入检查
if [[ ! -s "${HMM_LIST}" ]]; then
  log "错误：HMM 列表 \"${HMM_LIST}\" 不存在或为空。"
  exit 1
fi

if [[ ! -s "${SEQUENCE_LIST}" ]]; then
  log "错误：序列路径列表 \"${SEQUENCE_LIST}\" 不存在或为空。"
  exit 1
fi

# 读取 HMM 路径，过滤空行/回车
mapfile -t HMM_PATHS < <(sed 's/\r$//' "${HMM_LIST}" | grep -v '^[[:space:]]*$' || true)
if [[ ${#HMM_PATHS[@]} -eq 0 ]]; then
  log "错误：\"${HMM_LIST}\" 中没有有效的 HMM 路径。"
  exit 1
fi

# 验证氨基酸序列文件都存在
while IFS= read -r raw_seq || [[ -n "${raw_seq}" ]]; do
  seq_path="${raw_seq//$'\r'/}"
  [[ -z "${seq_path//[[:space:]]/}" ]] && continue
  if [[ ! -f "${seq_path}" ]]; then
    log "错误：序列文件不存在：${seq_path}"
    exit 1
  fi
done < "${SEQUENCE_LIST}"

# 针对每个 HMM 文件依次运行 GToTree
for hmm_file in "${HMM_PATHS[@]}"; do
  hmm_path="${hmm_file}"
  if [[ ! -f "${hmm_path}" ]]; then
    log "错误：HMM 文件不存在：${hmm_path}"
    exit 1
  fi

  hmm_basename="$(basename "${hmm_path}")"
  hmm_tag="${hmm_basename%.hmm}"
  [[ -z "${hmm_tag}" ]] && hmm_tag="custom_hmm"

  run_prefix="${OUTPUT_ROOT}/${hmm_tag}"
  mkdir -p "${run_prefix}"

  run_id="run_$(date +%Y%m%d_%H%M%S)"
  run_dir="${run_prefix}/${run_id}"
  while [[ -d "${run_dir}" ]]; do
    sleep 1
    run_id="run_$(date +%Y%m%d_%H%M%S)"
    run_dir="${run_prefix}/${run_id}"
  done

  # 若 HMM 缺少 GA 行（供 --cut_ga 使用），自动补上一条 GA 0 0;
  hmm_to_use="${hmm_path}"
  if ! grep -q '^GA' "${hmm_path}"; then
    tmp_hmm="$(mktemp -t "${hmm_tag}.XXXXXX.hmm")"
    awk '
      BEGIN { added = 0 }
      $1 == "GA" {
        added = 1
        print
        next
      }
      $1 == "HMM" && added == 0 {
        printf "GA      0 0;\n"
        added = 1
        print
        next
      }
      {
        print
      }
      END {
        if (added == 0) {
          printf "GA      0 0;\n"
        }
      }
    ' "${hmm_path}" > "${tmp_hmm}"
    hmm_to_use="${tmp_hmm}"
    TMP_HMM_FILES+=("${tmp_hmm}")
    log "HMM：${hmm_basename} 缺少 GA cutoff，已自动补上 \"GA 0 0;\"。"
  fi

  log "开始处理 HMM：${hmm_basename}，输出目录：${run_dir}"

  gtotree_args=(
    -A "${SEQUENCE_LIST}"
    -H "${hmm_to_use}"
    -o "${run_dir}"
    -j "${PARALLEL_JOBS}"
    -n "${HMM_CPUS}"
    -M "${MUSCLE_THREADS}"
    -N
    -k
  )

  if [[ "${BEST_HIT_MODE}" == "true" ]]; then
    gtotree_args+=(-B)
  fi

  GToTree "${gtotree_args[@]}"

  log "HMM：${hmm_basename} 处理完成，结果位于：${run_dir}"
done

log "所有 HMM 搜索任务完成。"
