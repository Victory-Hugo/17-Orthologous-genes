#!/bin/bash
set -euo pipefail

# 按样本循环运行 Rv0194 筛选
# 路径固定配置（如需调整请直接编辑本脚本）
pfam_dir="/mnt/l/18-Rv0194-Gene/1-BLAST/output/pfamA"
tmhmm_dir="/mnt/l/18-Rv0194-Gene/1-BLAST/output/tmhmm"
output_root="/mnt/l/18-Rv0194-Gene/1-BLAST/output/pfamA-filter"
# 如无需邻域过滤，可将 neighbor_tsv/neighbor_keywords 置为 "none"
fasta_dir="none"
neighbor_tsv="none"
neighbor_keywords="none"
python_bin="python"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

rv0194_filter="${script_dir}/rv0194_filter.py"

echo "[INFO] 批量运行 Rv0194 筛选..."
echo "[INFO] Pfam 目录: ${pfam_dir}"
echo "[INFO] TMHMM 目录: ${tmhmm_dir}"
echo "[INFO] FASTA 目录: ${fasta_dir}"
echo "[INFO] 输出根目录: ${output_root}"
if [[ -n "${neighbor_tsv}" ]]; then
    echo "[INFO] 邻域 TSV: ${neighbor_tsv}"
    [[ -n "${neighbor_keywords}" ]] && echo "[INFO] 邻域关键词: ${neighbor_keywords}"
fi

shopt -s nullglob
pfam_files=("${pfam_dir}"/*_pfam.csv)
if [[ ${#pfam_files[@]} -eq 0 ]]; then
    echo "[WARN] 未在 ${pfam_dir} 发现 *_pfam.csv，退出。" >&2
    exit 1
fi

combined_abc="${output_root}/abc_tran_summary.tsv"
combined_tm="${output_root}/tm_summary.tsv"
combined_report="${output_root}/detailed_report.tsv"
combined_candidates="${output_root}/rv0194_like_candidates.txt"

rm -f "${combined_abc}" "${combined_tm}" "${combined_report}" "${combined_candidates}"
mkdir -p "${output_root}"
headers_written=false

for pfam_csv in "${pfam_files[@]}"; do
    sample="$(basename "${pfam_csv}" "_pfam.csv")"
    tmhmm_file="${tmhmm_dir}/${sample}.txt"

    # 找到样本 FASTA（常见扩展名）
    fasta_path=""
    for ext in .faa .fa .fasta; do
        cand="${fasta_dir}/${sample}${ext}"
        if [[ -f "${cand}" ]]; then
            fasta_path="${cand}"
            break
        fi
    done

    echo "[INFO] 样本 ${sample}"
    echo "       Pfam:   ${pfam_csv}"
    echo "       TMHMM:  ${tmhmm_file}"
    echo "       FASTA:  ${fasta_path:-未找到，将跳过基序检查}"
    echo "       输出:   汇总到 ${output_root}"

    if [[ ! -f "${tmhmm_file}" ]]; then
        echo "[WARN] 缺少 TMHMM 文件，跳过样本 ${sample}" >&2
        continue
    fi

    tmp_dir="$(mktemp -d)"
    cmd=(
        "${python_bin}" "${rv0194_filter}"
        --hmmscan-csv "${pfam_csv}"
        --tmhmm "${tmhmm_file}"
        --output-dir "${tmp_dir}"
    )

    if [[ -n "${fasta_path}" ]]; then
        cmd+=(--fasta "${fasta_path}" --require-motif)
    fi

    if [[ -n "${neighbor_tsv}" ]]; then
        cmd+=(--neighbor-tsv "${neighbor_tsv}")
        if [[ -n "${neighbor_keywords}" ]] && [[ "${neighbor_keywords}" != "none" ]]; then
            cmd+=(--neighbor-keywords "${neighbor_keywords}")
        fi
    fi

    echo "[INFO] command: ${cmd[*]}"
    if ! "${cmd[@]}"; then
        echo "[WARN] 样本 ${sample} 运行失败，跳过汇总" >&2
        rm -rf "${tmp_dir}"
        continue
    fi

    if [[ "${headers_written}" = false ]]; then
        cat "${tmp_dir}/abc_tran_summary.tsv" > "${combined_abc}"
        cat "${tmp_dir}/tm_summary.tsv" > "${combined_tm}"
        cat "${tmp_dir}/detailed_report.tsv" > "${combined_report}"
        headers_written=true
    else
        tail -n +2 "${tmp_dir}/abc_tran_summary.tsv" >> "${combined_abc}"
        tail -n +2 "${tmp_dir}/tm_summary.tsv" >> "${combined_tm}"
        tail -n +2 "${tmp_dir}/detailed_report.tsv" >> "${combined_report}"
    fi
    cat "${tmp_dir}/rv0194_like_candidates.txt" >> "${combined_candidates}"

    rm -rf "${tmp_dir}"
done

echo "[INFO] 全部样本处理完成，汇总文件位于 ${output_root}。"
