#!/bin/bash


# 功能简介:
# 本脚本用于批量处理蛋白质序列（.faa 文件），依次进行 Pfam-A 域注释、跨膜区段预测（TMHMM）、信号肽预测（SignalP）适用于高通量蛋白功能注释与筛选分析。

# 主要流程:
# 1. 参数与路径初始化，支持命令行参数或默认路径。
# 2. 检查 TMHMM 和 SignalP 可执行文件是否可用。
# 3. 查找输入目录下所有 .faa 文件。
# 4. Pfam-A 扫描：调用 hmmscan 对每个 .faa 文件进行 Pfam-A 域注释，并格式化输出为 CSV。
# 5. 合并所有 Pfam-A CSV 文件为一个总表。
# 6. TMHMM 批处理：并行预测所有蛋白的跨膜区段，合并并解析结果为标准格式。
# 7. SignalP 批处理：并行预测所有蛋白的信号肽，输出结果。

script_dir="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/6-Pfam-A"
format_clean="${script_dir}/pipe/domtblout_to_csv.py"
pfam_db_default="${script_dir}/data/Pfam-A.hmm"
python="python"
tmhmm_bin_default="tmhmm"
signalp_bin_default="signalp4"

input_dir_default="/mnt/l/18-Rv0194-Gene/1-BLAST/output/sequence"
pfam_out_default="/mnt/l/18-Rv0194-Gene/1-BLAST/output/pfamA"
tmhmm_out_default="/mnt/l/18-Rv0194-Gene/1-BLAST/output/tmhmm"
signal_out_default="/mnt/l/18-Rv0194-Gene/1-BLAST/output/signal"
filter_out_default="/mnt/l/18-Rv0194-Gene/1-BLAST/output/pfamA-filter"
threads_default=8

input_dir="${1:-${input_dir_default}}"
pfam_db="${2:-${pfam_db_default}}"
threads="${3:-${threads_default}}"
pfam_out="${4:-${pfam_out_default}}"
tmhmm_out="${5:-${tmhmm_out_default}}"
signal_out="${6:-${signal_out_default}}"
filter_out="${7:-${filter_out_default}}"
tmhmm_bin="${8:-${tmhmm_bin_default}}"
signalp_bin="${9:-${signalp_bin_default}}"

echo "[INFO] 输入目录: ${input_dir}"
echo "[INFO] Pfam DB: ${pfam_db}"
echo "[INFO] 线程数: ${threads}"
echo "[INFO] 输出目录: pfam=${pfam_out}, tmhmm=${tmhmm_out}, signal=${signal_out}, combined=${filter_out}"
echo "[INFO] 可执行: tmhmm=${tmhmm_bin}, signalp=${signalp_bin}"

if ! command -v "${tmhmm_bin}" >/dev/null 2>&1; then
    echo "未找到 TMHMM 可执行文件：${tmhmm_bin}" >&2
    exit 1
fi
if ! command -v "${signalp_bin}" >/dev/null 2>&1; then
    echo "未找到 SignalP 可执行文件：${signalp_bin}" >&2
    exit 1
fi

mapfile -t faa_files < <(find "${input_dir}" -maxdepth 1 -type f -name "*.faa" | sort)
if [[ ${#faa_files[@]} -eq 0 ]]; then
    echo "未在 ${input_dir} 下找到 .faa 文件" >&2
    exit 1
fi

# 为了在输出中保留来源文件信息（如 GCA/GCF），给每条序列 header 前加上文件名前缀，使用临时目录避免改动原始文件。
prefixed_dir="$(mktemp -d)"
prefixed_faa_files=()
for faa in "${faa_files[@]}"; do
    base="$(basename "${faa}" .faa)"
    prefixed_faa="${prefixed_dir}/${base}.faa"
    awk -v prefix="${base}|" 'BEGIN{OFS=""} /^>/{sub(/^>/,">"prefix)} {print}' "${faa}" > "${prefixed_faa}"
    prefixed_faa_files+=("${prefixed_faa}")
done
trap 'rm -rf "${prefixed_dir}"' EXIT

mkdir -p "${pfam_out}" "${tmhmm_out}" "${signal_out}" "${filter_out}"

# Pfam 扫描
echo "[INFO] 开始 Pfam 扫描..."
for faa in "${prefixed_faa_files[@]}"; do
    base="$(basename "${faa}" .faa)"
    hmmscan --cpu "${threads}" --domtblout "${pfam_out}/${base}_pfam.domtblout" "${pfam_db}" "${faa}"
done
${python} "${format_clean}" -i "${pfam_out}" -o "${pfam_out}"


# 合并 Pfam CSV
combined_csv="${filter_out}/pfam_combined.csv"
first_csv=true
for csv_file in "${pfam_out}"/*.csv; do
    [ -e "${csv_file}" ] || continue
    if ${first_csv}; then
        cp "${csv_file}" "${combined_csv}"
        first_csv=false
    else
        tail -n +2 "${csv_file}" >> "${combined_csv}"
    fi
done

# TMHMM 批处理
echo "[INFO] 运行 TMHMM..."
tmhmm_tmp="${tmhmm_out}/tmp"
mkdir -p "${tmhmm_tmp}"
process_tmhmm() {
    local faa_path="$1"
    local stem="${faa_path##*/}"
    stem="${stem%.*}"
    local temp_output
    temp_output="$(mktemp "${tmhmm_tmp}/${stem}.XXXXXX")"
    local final_output="${tmhmm_out}/${stem}.txt"
    "${tmhmm_bin}" --short < "${faa_path}" > "${temp_output}"
    mv -f "${temp_output}" "${final_output}"
}
export tmhmm_bin tmhmm_tmp tmhmm_out
export -f process_tmhmm
parallel --jobs "${threads}" process_tmhmm ::: "${prefixed_faa_files[@]}"
# 合并 tmhmm 输出并解析出 TM 段
tmhmm_combined="${tmhmm_out}/tmhmm_short.out"
tmhmm_parsed="${tmhmm_out}/tmhmm_parsed.tsv"
cat "${tmhmm_out}"/*.txt > "${tmhmm_combined}"
awk '
{
    id=$1
    topo=""
    for(i=1;i<=NF;i++){
        if($i ~ /^Topology=/){sub(/^Topology=/,"",$i); topo=$i}
    }
    if(topo=="") next
    # 逐个匹配 (M|i|o)xx-yy 模式，Topology 可能无空格
    while (match(topo, /(M|i|o)[0-9]+-[0-9]+/)) {
        seg = substr(topo, RSTART, RLENGTH)
        topo = substr(topo, RSTART + RLENGTH)  # 继续匹配后续片段
        sub(/^[Mio]/, "", seg)
        split(seg, coords, "-")
        if(coords[1]!="" && coords[2]!=""){
            print id "\tTMhelix\t" coords[1] "\t" coords[2]
        }
    }
}
' "${tmhmm_combined}" > "${tmhmm_parsed}"
rm -rf "${tmhmm_tmp}" "${tmhmm_out}"/TMHMM_*

# SignalP 批处理
echo "[INFO] 运行 SignalP..."
signal_tmp="${signal_out}/tmp"
mkdir -p "${signal_tmp}"
process_signalp() {
    local faa_path="$1"
    local stem="${faa_path##*/}"
    stem="${stem%.*}"
    local temp_output
    temp_output="$(mktemp "${signal_tmp}/${stem}.XXXXXX")"
    local final_output="${signal_out}/${stem}.txt"
    "${signalp_bin}" -f short -t gram- -T "${signal_tmp}" "${faa_path}" > "${temp_output}"
    mv -f "${temp_output}" "${final_output}"
}
export signalp_bin signal_tmp signal_out
export -f process_signalp
parallel --jobs "${threads}" process_signalp ::: "${prefixed_faa_files[@]}"
rm -rf "${signal_tmp}"
