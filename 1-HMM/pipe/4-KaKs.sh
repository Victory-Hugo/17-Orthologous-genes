#!/bin/bash
set -euo pipefail

################################################################################
# 从现有的密码子比对开始，分步执行：
# 1) 建立 *.codon.fna 软链
# 2) 构建 AXT
# 3) 运行 KaKs_Calculator
# 4) 合并 KaKs
# 5) 统计 & 绘图
# 依赖：KaKs_Calculator、python3（pandas、matplotlib）
################################################################################

SCRIPT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM"
PY_PREP="${SCRIPT_DIR}/script/kaks_prepare.py"
PY_AXT="${SCRIPT_DIR}/script/6-1-构建KaKs参考AXT.py"
PY_RUN="${SCRIPT_DIR}/script/kaks_run.py"
PY_MERGE="${SCRIPT_DIR}/script/kaks_merge.py"
PY_STATS="${SCRIPT_DIR}/script/7-1-统计KaKs结果.py"
PY_PLOT="${SCRIPT_DIR}/script/7-2-绘制KaKs分布.py"

#* ==================== 用户可修改区域 ====================
# 输入/输出目录（已包含 FimH 的密码子比对）
PIPEDIR="/mnt/d/5-NCBI-Reference/hmm分析示例/output/FimH.aln"
CODON_DIR="${PIPEDIR}/alignment"        #! 内含 *.codon.aln.fna（pal2nal 结果）
AXT_REF_DIR="${PIPEDIR}/axt"
KAKS_OUTPUT_DIR="${PIPEDIR}/kaks"
KAKS_STATS_DIR="${PIPEDIR}/kaks_stats"
KAKS_PLOT_DIR="${PIPEDIR}/kaks_plots"

# KaKs_Calculator
KAKS_BIN="${SCRIPT_DIR}/bin/KaKs"
KAKS_CPU=16                               # 并行运行 KaKs 任务的最大并发数（单任务依旧单线程）

# Query & 参考（参考名需与 FASTA header 的首个字段一致）
QUERIES=("FimH.aln")
REFERENCE="GCF_030487725.1_1|id=WP_301697224.1|score=30.20|evalue=5.70e-08|len=157|desc=fimbrial"

# KaKs 绘图
KAKS_METHOD=""  # 可选：指定 Method 名称

# KaKs 遗传密码表（对应 KaKs_Calculator 的 -c 选项，留空使用默认 1）
KAKS_CODON_TABLE=11

# 拆分 AXT（按 block 数）。0 表示不拆分；>0 则对每个 AXT 拆分后并发跑 KaKs
KAKS_SPLIT_BLOCKS=10
#* ==================== 用户可修改区域 ====================

log() {
    local level="$1"; shift
    echo "[$(date '+%F %T')][$level] $*"
}

need_file() {
    if [[ ! -e "$1" ]]; then
        log ERROR "缺少文件/目录: $1"
        exit 1
    fi
}

log INFO "使用脚本根目录: $SCRIPT_DIR"
need_file "$CODON_DIR"
need_file "$KAKS_BIN"

########################################
# 1) 准备 *.codon.fna 软链
########################################
log INFO "Step1: 创建密码子比对软链 (*.codon.fna)"
python3 "$PY_PREP" --codon-dir "$CODON_DIR"

########################################
# 2) 生成 AXT
########################################
log INFO "Step2: 构建 AXT -> $AXT_REF_DIR"
python3 "$PY_AXT" \
  --alignment-dir "$CODON_DIR" \
  --output-dir "$AXT_REF_DIR" \
  --reference "$REFERENCE" \
  --queries "${QUERIES[@]}"

########################################
# 3) 跑 KaKs（逐个 AXT，必要时拆分并并发）
########################################
log INFO "Step3: 运行 KaKs_Calculator -> $KAKS_OUTPUT_DIR"
for axt in "$AXT_REF_DIR"/*.axt; do
  [[ -e "$axt" ]] || continue
  run_cmd=(python3 "$PY_RUN" --axt "$axt" --kaks-dir "$KAKS_OUTPUT_DIR" --kaks-bin "$KAKS_BIN" --kaks-cpu "$KAKS_CPU")
  if [[ -n "${KAKS_CODON_TABLE:-}" ]]; then
    run_cmd+=(--kaks-codon-table "$KAKS_CODON_TABLE")
  fi
  if [[ "${KAKS_SPLIT_BLOCKS:-0}" -gt 0 ]]; then
    run_cmd+=(--split-blocks "$KAKS_SPLIT_BLOCKS")
  fi
  log INFO "  AXT: $axt"
  "${run_cmd[@]}"
done

########################################
# 4) 合并 KaKs 结果
########################################
log INFO "Step4: 合并 KaKs 结果"
python3 "$PY_MERGE" \
  --kaks-dir "$KAKS_OUTPUT_DIR" \
  --output "${KAKS_OUTPUT_DIR}/kaks_merged.tsv" \
  --keep-intermediate

########################################
# 5) 统计
########################################
log INFO "Step5: 统计 KaKs -> $KAKS_STATS_DIR"
python3 "$PY_STATS" \
  --kaks-dir "$KAKS_OUTPUT_DIR" \
  --output-dir "$KAKS_STATS_DIR"

########################################
# 6) 绘图
########################################
log INFO "Step6: 绘制 KaKs 分布 -> $KAKS_PLOT_DIR"
plot_cmd=(python3 "$PY_PLOT" --summary "${KAKS_STATS_DIR}/KaKs_summary.tsv" --output-dir "$KAKS_PLOT_DIR")
if [[ -n "$KAKS_METHOD" ]]; then
  plot_cmd+=("--method" "$KAKS_METHOD")
fi
"${plot_cmd[@]}"

log INFO "全部完成。关键输出："
log INFO "  AXT 目录: $AXT_REF_DIR"
log INFO "  KaKs: $KAKS_OUTPUT_DIR"
log INFO "  统计: $KAKS_STATS_DIR"
log INFO "  图像: $KAKS_PLOT_DIR"
