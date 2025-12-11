: '
脚本功能说明：
本脚本用于根据映射文件（map.txt）批量提取目标蛋白的核酸序列（fna），并输出到指定目录。支持多线程并发处理。

主要流程：
1. 读取映射文件，每行包含组装ID、蛋白ID、fasta路径、gff路径。
2. 生成目标蛋白ID列表，用于后续gff文件过滤。
3. 对每一行映射信息，执行以下操作：
    - 检查fasta和gff文件是否存在。
    - 过滤gff文件，提取目标蛋白注释行。
    - 使用gffread根据gff注释和fasta文件提取对应的核酸序列。
    - 输出为fna文件，文件名为组装ID.fna，序列头格式为 >组装ID|蛋白ID。
    - 支持断点续跑，已存在且非空的fna文件会跳过。
4. 支持GNU parallel或xargs多线程加速。

参数说明：
- map_file：映射文件路径，格式为csv，包含组装ID、蛋白ID、fasta路径、gff路径。
- out_dir：输出目录。
- threads：并发线程数。

'
#!/usr/bin/env bash
set -euo pipefail

map_file="/mnt/l/18-Rv0194-Gene/3-SingleGene-tree/data/map.txt"
out_dir="/mnt/l/18-Rv0194-Gene/3-SingleGene-tree/data/Rv0194"
threads="16"

mkdir -p "$out_dir"

protein_list="$(mktemp)"
trap 'rm -f "$protein_list"' EXIT
cut -d, -f2 "$map_file" | sed '/^$/d' > "$protein_list"

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

    tmp_dir="$(mktemp -d)"
    filtered_gff="${tmp_dir}/filtered.gff"
    target_gff="${tmp_dir}/target.gff"

    # 一次过滤gff获取所有目标蛋白行
    grep -F -f "$protein_list" "$gff_path" > "$filtered_gff" || true
    # 针对当前蛋白再切分
    grep -F "$protein_id" "$filtered_gff" > "$target_gff" || true

    out_fna="${out_dir}/${asm_id}.fna"
    header=">${asm_id}|${protein_id}"

    # 断点续跑：已有非空fna则跳过
    if [[ -s "$out_fna" ]]; then
        echo "已存在，跳过: $out_fna"
        rm -rf "$tmp_dir"
        return 0
    fi

    if [[ -s "$target_gff" ]]; then
        gffread "$target_gff" -g "$fasta_path" -x "$out_fna"
        if [[ -s "$out_fna" ]]; then
            { echo "$header"; tail -n +2 "$out_fna"; } > "${out_fna}.tmp" && mv "${out_fna}.tmp" "$out_fna"
        else
            echo "gffread未生成序列: ${asm_id}" >&2
            rm -f "$out_fna"
        fi
    else
        echo "未找到蛋白${protein_id}的注释: ${asm_id}" >&2
    fi

    rm -rf "$tmp_dir"
}
export -f process_line
export protein_list out_dir

readarray -t map_lines < "$map_file"

if command -v parallel >/dev/null 2>&1; then
    printf '%s\n' "${map_lines[@]}" | parallel --halt soon,fail=1 -j "$threads" process_line {}
else
    printf '%s\n' "${map_lines[@]}" | xargs -P "$threads" -I{} bash -c 'process_line "$@"' _ {}
fi
