#!/bin/bash
set -euo pipefail

# 独立运行 Rv0194 筛选
# 参数（位置）：
#   1=hmmscan_csv 2=tmhmm_parsed 3=output_dir 4=python_bin

python_bin="${4:-python}"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

hmmscan_csv="${1:-/mnt/l/18-Rv0194-Gene/1-BLAST/output/pfamA-filter/pfam_combined.csv}"
tmhmm_parsed="${2:-/mnt/l/18-Rv0194-Gene/1-BLAST/output/tmhmm/tmhmm_parsed.tsv}"
output_dir="${3:-/mnt/l/18-Rv0194-Gene/1-BLAST/output/pfamA-filter}"

rv0194_filter="${script_dir}/rv0194_filter.py"

echo "[INFO] 运行 Rv0194 筛选..."
echo "[INFO] hmmscan CSV: ${hmmscan_csv}"
echo "[INFO] TMHMM parsed: ${tmhmm_parsed}"
echo "[INFO] 输出目录: ${output_dir}"

"${python_bin}" "${rv0194_filter}" \
    --hmmscan-csv "${hmmscan_csv}" \
    --tmhmm "${tmhmm_parsed}" \
    --output-dir "${output_dir}"

echo "[INFO] 完成。结果目录：${output_dir}"
