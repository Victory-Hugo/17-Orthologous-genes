#!/usr/bin/env bash
set -euo pipefail

# ===============================================
# Reciprocal Best Hit (RBH) 流程 - BLASTP 版
# 用法: ./1-rbh.sh [threads]
# ===============================================

# ---------- 配置加载 ----------
CONF_FILE="/mnt/f/onedrive/文档（科研）/脚本/Download/17-Orthologous-genes/7-RBH/conf/Rv1819c-RBH.conf"
[[ -f "$CONF_FILE" ]] || { echo "错误: 配置文件不存在: $CONF_FILE" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONF_FILE"

# ---------- 参数覆盖 ----------
THREADS=${1:-$THREADS}
FORCE=${FORCE:-0}

# ---------- 前置检查 ----------
command -v makeblastdb >/dev/null 2>&1 || { echo "错误: 未找到 makeblastdb，请先安装 BLAST+。" >&2; exit 1; }
command -v blastp >/dev/null 2>&1 || { echo "错误: 未找到 blastp，请先安装 BLAST+。" >&2; exit 1; }

[[ -f "$QUERY_FASTA" ]] || { echo "错误: 查询序列不存在: $QUERY_FASTA" >&2; exit 1; }
[[ -f "$REF_PROTEOME" ]] || { echo "错误: 参考蛋白组不存在: $REF_PROTEOME" >&2; exit 1; }
[[ -d "$CANDIDATE_DIR" ]] || { echo "错误: 候选蛋白目录不存在: $CANDIDATE_DIR" >&2; exit 1; }

mkdir -p "$OUTPUT_DIR" "$DB_DIR"

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

# ---------- 通用建库函数 ----------
build_blast_db() {
    local in_faa="$1"
    local out_prefix="$2"
    makeblastdb -in "$in_faa" -dbtype prot -parse_seqids -out "$out_prefix" >/dev/null
}

# ---------- 从原始 FASTA 提取序列 ----------
extract_fasta_by_id() {
    local seq_id="$1"
    local fasta_file="$2"
    local out_file="$3"
    awk -v id="$seq_id" '
        BEGIN{found=0}
        /^>/{
            if(found){exit}
            header=$0
            sub(/^>/,"",header)
            split(header,a,/[^A-Za-z0-9_.|:-]/)
            cur=a[1]
            found=(cur==id)
        }
        found{print}
        END{exit(found?0:1)}
    ' "$fasta_file" > "$out_file"
}

# ---------- 标准化 BLAST 序列 ID ----------
normalize_blast_id() {
    local raw_id="$1"
    if [[ "$raw_id" == *"|"* ]]; then
        IFS='|' read -r tag acc rest <<<"$raw_id"
        case "$tag" in
            ref|gb|sp|tr|emb|dbj|tpg|tpe|tpd|pir|prf)
                if [[ -n "${acc:-}" ]]; then
                    echo "$acc"
                    return
                fi
                ;;
        esac
    fi
    echo "$raw_id"
}

# ---------- 获取 Rv1819c 查询 ID ----------
if [[ -n "${RV_QUERY_HEADER:-}" ]]; then
    RV_QUERY_ID=$(echo "$RV_QUERY_HEADER" | awk '{print $1}')
else
    RV_QUERY_ID=$(awk 'NR==1{gsub(/^>/,"");split($0,a,/[^A-Za-z0-9_.|:-]/);print a[1];exit}' "$QUERY_FASTA")
fi
if [[ -z "$RV_QUERY_ID" ]]; then
    echo "错误: 无法解析 Rv1819c FASTA ID。" >&2
    exit 1
fi

# ---------- 建立参考蛋白组数据库 ----------
REF_DB="${DB_DIR}/ref_proteome"
if [[ ! -f "${REF_DB}.pin" ]]; then
    echo "▶ 构建参考蛋白组数据库..."
    build_blast_db "$REF_PROTEOME" "$REF_DB"
fi

# ---------- 解析参考蛋白组中 Rv1819c 的 ID ----------
if [[ -z "${RV_REF_ID:-}" ]]; then
    echo "▶ 解析参考蛋白组中的 Rv1819c ID..."
    rv_map_out="${TMP_DIR}/rv1819c_in_ref.tsv"
    blastp \
        -query "$QUERY_FASTA" \
        -db "$REF_DB" \
        -evalue "$EVALUE" \
        -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
        -num_threads "$THREADS" \
        > "$rv_map_out" || true
    RV_REF_ID=$(sort -k12,12gr -k11,11g "$rv_map_out" | \
        awk -v qcov="$QCOV_THRESHOLD" -v scov="$SCOV_THRESHOLD" 'BEGIN{OFS="\t"} {qcovp=$4/$5*100; scovp=$4/$6*100; if($11<=1e-10 && qcovp>=qcov && scovp>=scov){print $2; exit}}')
    if [[ -z "$RV_REF_ID" ]]; then
        echo "⚠ 未能在参考蛋白组中解析 Rv1819c ID，回退到 QUERY ID: $RV_QUERY_ID" >&2
        RV_REF_ID="$RV_QUERY_ID"
    else
        echo "  ✓ Rv1819c 参考 ID: $RV_REF_ID"
    fi
fi

# ---------- 断点续跑准备 ----------
declare -A DONE
if [[ "$FORCE" -eq 1 ]]; then
    echo "▶ FORCE=1，重新计算并覆盖结果..."
    echo -e "genome_faa\tforward_best_hit\tforward_evalue\treverse_best_hit\treverse_evalue\tforward_qcov\treverse_qcov\tRBH_status" > "$OUT_TSV"
elif [[ -s "$OUT_TSV" ]]; then
    echo "▶ 检测到已有结果，启用断点续跑..."
    while IFS=$'\t' read -r genome_faa _; do
        [[ -z "$genome_faa" || "$genome_faa" == "genome_faa" ]] && continue
        DONE["$genome_faa"]=1
    done < "$OUT_TSV"
else
    echo -e "genome_faa\tforward_best_hit\tforward_evalue\treverse_best_hit\treverse_evalue\tforward_qcov\treverse_qcov\tRBH_status" > "$OUT_TSV"
fi

# ---------- 批量 RBH ----------
shopt -s nullglob
faa_files=("$CANDIDATE_DIR"/*.faa)
if [[ ${#faa_files[@]} -eq 0 ]]; then
    echo "错误: 目录中未找到 .faa 文件: $CANDIDATE_DIR" >&2
    exit 1
fi

for faa_file in "${faa_files[@]}"; do
    filename="$(basename "$faa_file")"
    if [[ -n "${DONE[$filename]:-}" ]]; then
        echo "  ↩︎ 跳过（已完成）: $filename"
        continue
    fi

    echo "▶ 处理: $filename"
    db_prefix="${DB_DIR}/${filename%.faa}"
    if [[ ! -f "${db_prefix}.pin" ]]; then
        build_blast_db "$faa_file" "$db_prefix"
    fi

    forward_out="${TMP_DIR}/${filename}.forward.tsv"
    blastp \
        -query "$QUERY_FASTA" \
        -db "$db_prefix" \
        -evalue "$EVALUE" \
        -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
        -num_threads "$THREADS" \
        > "$forward_out" || true

    forward_best=$(sort -k12,12gr -k11,11g "$forward_out" | \
        awk -v qcov="$QCOV_THRESHOLD" -v scov="$SCOV_THRESHOLD" 'BEGIN{OFS="\t"} {qcovp=$4/$5*100; scovp=$4/$6*100; if($11<=1e-10 && qcovp>=qcov && scovp>=scov){print $0, qcovp; exit}}')

    if [[ -z "$forward_best" ]]; then
        echo -e "$filename\tNA\tNA\tNA\tNA\tNA\tNA\tFALSE" >> "$OUT_TSV"
        continue
    fi

    forward_sseqid=$(echo "$forward_best" | awk '{print $2}')
    forward_evalue=$(echo "$forward_best" | awk '{print $11}')
    forward_qcov=$(echo "$forward_best" | awk '{print $13}')

    # 反向 BLAST：从原始 FASTA 提取 best hit 序列
    hit_faa="${TMP_DIR}/${filename}.hit.faa"
    if ! extract_fasta_by_id "$forward_sseqid" "$faa_file" "$hit_faa"; then
        echo "  ⚠ 无法从原始 FASTA 提取序列: $forward_sseqid" >&2
        echo -e "$filename\t$forward_sseqid\t$forward_evalue\tNA\tNA\t$forward_qcov\tNA\tFALSE" >> "$OUT_TSV"
        continue
    fi

    reverse_out="${TMP_DIR}/${filename}.reverse.tsv"
    blastp \
        -query "$hit_faa" \
        -db "$REF_DB" \
        -evalue "$EVALUE" \
        -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
        -num_threads "$THREADS" \
        > "$reverse_out" || true

    reverse_best=$(sort -k12,12gr -k11,11g "$reverse_out" | \
        awk -v qcov="$QCOV_THRESHOLD" -v scov="$SCOV_THRESHOLD" 'BEGIN{OFS="\t"} {qcovp=$4/$5*100; scovp=$4/$6*100; if($11<=1e-10 && qcovp>=qcov && scovp>=scov){print $0, qcovp; exit}}')

    if [[ -z "$reverse_best" ]]; then
        echo -e "$filename\t$forward_sseqid\t$forward_evalue\tNA\tNA\t$forward_qcov\tNA\tFALSE" >> "$OUT_TSV"
        continue
    fi

    reverse_sseqid=$(echo "$reverse_best" | awk '{print $2}')
    reverse_sseqid_norm=$(normalize_blast_id "$reverse_sseqid")
    reverse_evalue=$(echo "$reverse_best" | awk '{print $11}')
    reverse_qcov=$(echo "$reverse_best" | awk '{print $13}')

    if [[ "$reverse_sseqid_norm" == "$RV_REF_ID" ]]; then
        rbh_status="TRUE"
    else
        rbh_status="FALSE"
    fi

    echo -e "$filename\t$forward_sseqid\t$forward_evalue\t$reverse_sseqid\t$reverse_evalue\t$forward_qcov\t$reverse_qcov\t$rbh_status" >> "$OUT_TSV"

done

echo ""
echo "✓ RBH 流程完成。"
echo "输出: $OUT_TSV"
