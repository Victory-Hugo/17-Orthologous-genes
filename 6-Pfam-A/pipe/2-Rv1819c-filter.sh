#!/bin/bash
set -euo pipefail

pfam_csv="/mnt/l/19-Rv1819c-Gene/1-BLAST/output/pfamA-merge/pfam_combined.csv"
tmhmm_merge="/mnt/l/19-Rv1819c-Gene/1-BLAST/output/tmhmm-merge/tmhmm_merge.tsv"
signal_merge="/mnt/l/19-Rv1819c-Gene/1-BLAST/output/signal-merge/signal_merge.tsv"
output_dir="/mnt/l/19-Rv1819c-Gene/1-BLAST/output/Final"
python_bin="python"

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
rv1819c_filter="${script_dir}/rv1819c_filter.py"

echo "[INFO] Rv1819c 过滤"
echo "[INFO] Pfam CSV: ${pfam_csv}"
echo "[INFO] TMHMM 合并: ${tmhmm_merge}"
echo "[INFO] SignalP 合并: ${signal_merge}"
echo "[INFO] 输出目录: ${output_dir}"

mkdir -p "${output_dir}"

exec "${python_bin}" "${rv1819c_filter}" \
    --pfam-csv "${pfam_csv}" \
    --tmhmm-merge "${tmhmm_merge}" \
    --signal-merge "${signal_merge}" \
    --output-dir "${output_dir}"
