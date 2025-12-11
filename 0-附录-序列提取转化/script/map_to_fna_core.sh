#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "用法: $0 <map_file> <out_dir> [threads] [feature_type]" >&2
    echo "示例: $0 map.csv /path/to/out 16 CDS" >&2
}

map_file="${1:-}"
out_dir="${2:-}"
threads="${3:-16}"
# 可选：限制第三列feature类型（如CDS），为空则不限制
feature_type="${4:-}"

if [[ -z "$map_file" || -z "$out_dir" ]]; then
    usage
    exit 1
fi

if [[ ! -f "$map_file" ]]; then
    echo "未找到映射文件: $map_file" >&2
    exit 1
fi
if ! command -v gffread >/dev/null 2>&1; then
    echo "未找到gffread，请先安装" >&2
    exit 1
fi

mkdir -p "$out_dir"

process_line() {
    local line="$1"
    IFS=, read -r asm_id protein_id fasta_path gff_path <<< "$line"

    [[ -z "${asm_id:-}" || -z "${protein_id:-}" || -z "${fasta_path:-}" || -z "${gff_path:-}" ]] && return 0

    if [[ ! -f "$fasta_path" ]]; then
        echo "缺少fasta文件: $fasta_path" >&2
        return 0
    fi
    if [[ ! -f "$gff_path" ]]; then
        echo "缺少gff文件: $gff_path" >&2
        return 0
    fi

    local tmp_dir
    tmp_dir="$(mktemp -d)"
    cleanup() { rm -rf "$tmp_dir"; }
    trap cleanup RETURN

    local target_gff="${tmp_dir}/target.gff"
    # 精确匹配属性字段，避免子串误伤，保留注释头；可选限制feature类型
    awk -F'\t' -v pid="$protein_id" -v ft="$feature_type" '
        /^#/ { print; next }
        (ft == "" || $3 == ft) && $9 ~ ("(^|;)(ID|Name|Parent)=" pid "([;]|$)") { print }
    ' "$gff_path" > "$target_gff"

    out_fna="${out_dir}/${asm_id}.fna"
    header=">${asm_id}|${protein_id}"

    # 断点续跑：已有非空fna则跳过
    if [[ -s "$out_fna" ]]; then
        echo "已存在，跳过: $out_fna"
        return 0
    fi

    if [[ ! -s "$target_gff" ]]; then
        echo "未找到蛋白${protein_id}的注释: ${asm_id}" >&2
        return 0
    fi

    local tmp_fna="${tmp_dir}/out.fna"
    if gffread "$target_gff" -g "$fasta_path" -x - | {
        read -r _first_line || exit 1
        printf '%s\n' "$header"
        cat
    } > "$tmp_fna"; then
        if [[ -s "$tmp_fna" ]]; then
            mv "$tmp_fna" "$out_fna"
        else
            echo "gffread未生成序列: ${asm_id}" >&2
            rm -f "$tmp_fna"
        fi
    else
        echo "gffread执行失败: ${asm_id}" >&2
        rm -f "$tmp_fna"
    fi
}
export -f process_line
export out_dir

if command -v parallel >/dev/null 2>&1; then
    parallel --halt soon,fail=1 --line-buffer --bar -j "$threads" process_line :::: "$map_file"
else
    while IFS= read -r line || [[ -n "$line" ]]; do
        printf '%s\0' "$line"
    done < "$map_file" | xargs -0 -P "$threads" -I{} bash -c 'process_line "$@"' _ "{}"
fi
