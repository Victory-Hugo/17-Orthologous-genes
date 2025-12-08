#!/bin/bash
set -euo pipefail
#*================================
#* 配置区
#*================================
FAA_FNA_DIR="/mnt/d/5-NCBI-Reference/hmm分析示例/output/FimH.aln/sequence" #? 输入目录（同时包含 .hits.faa 与 .hits.cds.fna）
OUTPUT_DIR="/mnt/d/5-NCBI-Reference/hmm分析示例/output/FimH.aln/alignment" #? 输出目录（会自动创建）
MAFFT_BIN="mafft"           #? MAFFT 可执行文件路径
PAL2NAL_BIN="pal2nal.pl"    #? PAL2NAL 可执行文件路径
CODON_TABLE=11              #? 密码子表：11 = Bacterial, archaeal and plant plastid code
#*================================

#######################################
# 主逻辑
#######################################
log() {
  local level="$1"; shift
  echo "[$(date '+%F %T')][$level] $*"
}

log INFO "输入目录: $FAA_FNA_DIR"
log INFO "输出目录: $OUTPUT_DIR"

if [[ ! -d "$FAA_FNA_DIR" ]]; then
  log ERROR "输入目录不存在: $FAA_FNA_DIR"
  exit 1
fi

for bin in "$MAFFT_BIN" "$PAL2NAL_BIN"; do
  if ! command -v "$bin" >/dev/null 2>&1; then
    log ERROR "未找到可执行文件: $bin"
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

# 以父目录名作为 query 名（示例: FimH.aln）
QUERY_NAME="$(basename "$(dirname "$OUTPUT_DIR")")"
PROTEIN_FASTA="${OUTPUT_DIR}/${QUERY_NAME}.combined.faa"
CDS_FASTA="${OUTPUT_DIR}/${QUERY_NAME}.combined.cds.fna"
PROTEIN_ALN="${OUTPUT_DIR}/${QUERY_NAME}.protein.aln.faa"
CODON_ALN="${OUTPUT_DIR}/${QUERY_NAME}.codon.aln.fna"

log INFO "收集蛋白与 CDS 序列..."
> "$PROTEIN_FASTA"
> "$CDS_FASTA"
mapfile -t FAA_FILES < <(find "$FAA_FNA_DIR" -maxdepth 1 -type f -name "*.hits.faa" | sort)

if [[ ${#FAA_FILES[@]} -eq 0 ]]; then
  log ERROR "未在 $FAA_FNA_DIR 中找到任何 *.hits.faa 文件"
  exit 1
fi

count=0
for faa in "${FAA_FILES[@]}"; do
  sample_base="$(basename "$faa" .hits.faa)"
  cds_file="${FAA_FNA_DIR}/${sample_base}.hits.cds.fna"
  if [[ ! -s "$cds_file" ]]; then
    log ERROR "缺少匹配的 CDS 文件: $cds_file"
    exit 1
  fi
  cat "$faa" >> "$PROTEIN_FASTA"
  cat "$cds_file" >> "$CDS_FASTA"
  count=$((count + 1))
done
log INFO "已合并 $count 个样本的蛋白与 CDS 序列。"

log INFO "Step1: MAFFT 蛋白多序列比对 -> $PROTEIN_ALN"
"$MAFFT_BIN" --localpair --maxiterate 1000 "$PROTEIN_FASTA" > "$PROTEIN_ALN"

log INFO "Step2: PAL2NAL 转换为密码子比对 -> $CODON_ALN"
"$PAL2NAL_BIN" "$PROTEIN_ALN" "$CDS_FASTA" -codontable "$CODON_TABLE" -output fasta > "$CODON_ALN"

log INFO "完成。关键输出："
log INFO "  蛋白 MSA: $PROTEIN_ALN"
log INFO "  密码子 MSA: $CODON_ALN"
