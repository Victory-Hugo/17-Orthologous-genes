#!/usr/bin/env bash
################################################################################
# 作者 (Author): BigLin
# 脚本用途: 对指定目录下所有 *.aln.faa 多序列比对文件批量构建 HMM (hmmbuild)
#
# 依赖软件 (Required software/packages):
#   - bash >= 4
#   - coreutils (find, xargs, sort, mv, rm, etc.)
#   - HMMER v3+ (提供 hmmbuild)
#   - GNU parallel (提供 parallel，支持 --bar 进度条)
#
# 运行说明:
#   1) 直接运行，无需参数；变量均已写死
#   2) 自动断点续跑：已存在的 .hmm 视为完成，将被跳过
#   3) 并行构建：可调 JOBS
#   4) 中断(SIGINT/TERM)时，自动删除未完成的临时产物，仅保留已完成结果
################################################################################

set -Eeuo pipefail

############################
# 配置区（全部写死的变量）  #
############################
INPUT_DIR="/home/luolintao/7-Rv0194/data/1-HMM序列"
OUTPUT_DIR="/home/luolintao/7-Rv0194/data/2-HMM库"
FILE_GLOB_PATTERN="*.aln.faa"         # 匹配输入 MSA 的文件名模式
JOBS=8                                 # 并行任务数（根据机器调整）
RESULTS_DIR="${OUTPUT_DIR}/.parallel_results"
WORKDIR="${OUTPUT_DIR}/.work_hmm"      # 用于跟踪未完成任务
JOBLOG="${OUTPUT_DIR}/parallel_hmmbuild.joblog"
LISTFILE="${OUTPUT_DIR}/.hmmbuild_input.list"
HMMPROG="hmmbuild"                     # hmmbuild 程序名（在 PATH 中）

############################
# 彩色打印设置（仅装饰用途） #
############################
C_RESET=$'\033[0m'
C_RED=$'\033[31m'
C_GRN=$'\033[32m'
C_YLW=$'\033[33m'
C_BLU=$'\033[34m'
C_MAG=$'\033[35m'
C_CYN=$'\033[36m'

info()  { printf "%s[INFO]%s %s\n"  "$C_CYN" "$C_RESET" "$*"; }
good()  { printf "%s[DONE]%s %s\n"  "$C_GRN" "$C_RESET" "$*"; }
warn()  { printf "%s[WARN]%s %s\n"  "$C_YLW" "$C_RESET" "$*"; }
err()   { printf "%s[ERR ]%s %s\n"  "$C_RED" "$C_RESET" "$*"; }

################################
# 依赖检查                      #
################################
need_bin() {
  if ! command -v "$1" >/dev/null 2>&1; then
    err "缺少可执行文件：$1 ；请安装后重试。"
    exit 127
  fi
}
need_bin "$HMMPROG"
need_bin parallel
need_bin find
need_bin xargs

################################
# 准备环境                      #
################################
mkdir -p "$RESULTS_DIR" "$WORKDIR"

# 生成待处理文件清单：仅列出还没有对应 .hmm 的输入
# 对于 foo.aln.faa -> 输出 foo.hmm
info "扫描输入目录：$INPUT_DIR"
> "$LISTFILE"
while IFS= read -r -d '' f; do
  base="$(basename "$f")"
  stem="${base%.aln.faa}"
  out="${OUTPUT_DIR}/${stem}.hmm"
  if [[ -s "$out" ]]; then
    # 已完成则跳过（断点续跑）
    :
  else
    printf "%s\0" "$f" >> "$LISTFILE"
  fi
done < <(find "$INPUT_DIR" -maxdepth 1 -type f -name "$FILE_GLOB_PATTERN" -print0)

TOTAL_ALL=$(find "$INPUT_DIR" -maxdepth 1 -type f -name "$FILE_GLOB_PATTERN" | wc -l | tr -d ' ')
TOTAL_TODO=$(( $(wc -c < "$LISTFILE") > 0 ? 1 : 0 ))
# 上面只是测试文件是否空；真正的待处理数如下：
TODO_COUNT=$(tr -cd '\0' < "$LISTFILE" | wc -c | tr -d ' ')
DONE_COUNT=$(( TOTAL_ALL - TODO_COUNT ))

info "总文件数: ${TOTAL_ALL}；已完成: ${DONE_COUNT}；待处理: ${TODO_COUNT}"
if [[ "$TODO_COUNT" -eq 0 ]]; then
  good "没有待处理任务。所有 HMM 已存在。"
  exit 0
fi

################################
# 中断/退出时清理未完成任务     #
################################
PARALLEL_PID=""
cleanup_unfinished() {
  warn "收到中断信号，开始清理未完成的临时产物…"
  # 终止 parallel 及其子进程（如果仍在运行）
  if [[ -n "${PARALLEL_PID}" ]] && ps -p "$PARALLEL_PID" >/dev/null 2>&1; then
    warn "终止并行任务 (PID=$PARALLEL_PID)…"
    pkill -P "$PARALLEL_PID" 2>/dev/null || true
    kill  "$PARALLEL_PID" 2>/dev/null || true
  fi
  # 删除 .work 标记及与之对应的 .hmm.tmp
  find "$WORKDIR" -type f -name "*.work" -print0 2>/dev/null | while IFS= read -r -d '' w; do
    stem="$(basename "$w" .work)"
    tmp="${OUTPUT_DIR}/${stem}.hmm.tmp"
    if [[ -e "$tmp" ]]; then
      rm -f -- "$tmp"
      warn "已删除未完成: ${stem}.hmm.tmp"
    fi
    rm -f -- "$w"
  done
  warn "清理完成。已完成的 .hmm 文件保持不变。"
}
trap cleanup_unfinished INT TERM

################################
# 单个任务函数                  #
################################
build_one() {
  local in="$1"
  local base stem out tmp work
  base="$(basename "$in")"
  stem="${base%.aln.faa}"
  out="${OUTPUT_DIR}/${stem}.hmm"
  tmp="${OUTPUT_DIR}/${stem}.hmm.tmp"
  work="${WORKDIR}/${stem}.work"

  # 如果目标已存在且非空，直接跳过（稳妥）
  if [[ -s "$out" ]]; then
    printf "%s[SKIP]%s %s 已存在\n" "$C_BLU" "$C_RESET" "$out"
    return 0
  fi

  # 创建工作标记
  : > "$work"

  # 构建 HMM -> 写到临时文件
  # --amino：hmmbuild 会自动识别氨基酸，这里不强制；如需严格，可加 --amino
  # --cpu 在 parallel 环境下建议设为 1，避免过度竞争
  if "$HMMPROG" --cpu 1 "$tmp" "$in" >/dev/null 2>&1; then
    mv -f -- "$tmp" "$out"
    rm -f -- "$work"
    printf "%s[DONE]%s %s -> %s\n" "$C_GRN" "$C_RESET" "$in" "$out"
    return 0
  else
    rm -f -- "$tmp" "$work"
    printf "%s[FAIL]%s %s\n" "$C_RED" "$C_RESET" "$in"
    return 1
  fi
}

export -f build_one
export OUTPUT_DIR WORKDIR HMMPROG C_BLU C_RESET C_GRN C_RED

################################
# 并行执行（带进度条与日志）    #
################################
info "开始并行构建 HMM（JOBS=${JOBS}）…"
# 使用 --bar 作为进度条；-0 从 NUL 分隔的列表读取；--joblog 记录状态便于排错与续跑
# --results 仅用于调试保留 stdout/stderr（可删）
set +e  # 允许个别任务失败但整体继续；失败将在退出码中体现
parallel --bar -0 --jobs "$JOBS" --joblog "$JOBLOG" --results "$RESULTS_DIR" \
  build_one :::: "$LISTFILE" &
PARALLEL_PID=$!
wait "$PARALLEL_PID"
EXITCODE=$?
set -e

################################
# 结果汇总                      #
################################
# 重新统计
NEW_DONE=$(find "$INPUT_DIR" -maxdepth 1 -type f -name "$FILE_GLOB_PATTERN" -printf "%f\n" \
  | sed 's/\.aln\.faa$/.hmm/' \
  | sed "s|^|${OUTPUT_DIR}/|" \
  | xargs -I{} bash -c '[[ -s "$1" ]] && echo 1 || echo 0' _ {} \
  | awk '{s+=$1} END{print s}')

info "构建完成：成功 ${NEW_DONE} / ${TOTAL_ALL}"
if [[ "$EXITCODE" -eq 0 ]]; then
  good "所有待处理任务均已成功完成。"
else
  warn "部分任务失败；详情见日志：$JOBLOG 与 $RESULTS_DIR"
  exit "$EXITCODE"
fi
