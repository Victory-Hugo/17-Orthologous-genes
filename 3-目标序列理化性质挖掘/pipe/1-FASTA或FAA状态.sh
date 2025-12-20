#!/bin/bash
#? 统计多个 FASTA 文件的多种理化性质指标，并将结果分别输出到指定目录中，用于分析基因在不同物种中的长度和理化特征变化。

set -euo pipefail

PYTHON_SCRIPT="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/3-目标序列理化性质挖掘/src/FASTA或faa状态.py"

FASTA_LIST_TXT="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/3-目标序列理化性质挖掘/conf/搜索完.list.txt"
OUTPUT_BASE="/mnt/l/18-Rv0194-Gene/3-SingleGene-tree/output/理化性质/"
LENGTH_DIR="${OUTPUT_BASE}/1-序列长度"
GC_DIR="${OUTPUT_BASE}/2-GC含量"
PI_DIR="${OUTPUT_BASE}/3-等电点"
HYDRO_DIR="${OUTPUT_BASE}/4-疏水性指数"
CAI_DIR="${OUTPUT_BASE}/5-CAI"

# 配置 CAI 计算所需的宿主高表达基因 CDS 参考序列（DNA）
CAI_REFERENCE_FASTA="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/3-目标序列理化性质挖掘/data/鲍曼参考_GCF_014672775.1_CDS.fna"   # 直接指定一个 FASTA 文件路径
CAI_REFERENCE_LIST=""    # 或者指定一个列出多个 FASTA 文件路径的文本文件

PROGRESS_WIDTH=40

if [[ ! -f "${PYTHON_SCRIPT}" ]]; then
    echo "缺少 Python 脚本: ${PYTHON_SCRIPT}" >&2
    exit 1
fi

if [[ ! -f "${FASTA_LIST_TXT}" ]]; then
    echo "FASTA 列表文件不存在: ${FASTA_LIST_TXT}" >&2
    exit 1
fi

TOTAL=$(grep -cve '^\s*$' "${FASTA_LIST_TXT}" || true)
if [[ "${TOTAL}" -eq 0 ]]; then
    echo "FASTA 列表文件为空: ${FASTA_LIST_TXT}" >&2
    exit 0
fi

if ! command -v python3 >/dev/null 2>&1; then
    echo "缺少软件: python3。请先安装 python3。" >&2
    exit 1
fi

if ! python3 -c "import Bio" >/dev/null 2>&1; then
    echo "缺少软件: biopython。请运行 'pip install biopython' 或 'conda install biopython' 进行安装。" >&2
    exit 1
fi

mkdir -p "${LENGTH_DIR}" "${GC_DIR}" "${PI_DIR}" "${HYDRO_DIR}" "${CAI_DIR}"

LENGTH_SUMMARY="${LENGTH_DIR}/length.tsv"
GC_SUMMARY="${GC_DIR}/gc_content.tsv"
PI_SUMMARY="${PI_DIR}/isoelectric_point.tsv"
HYDRO_SUMMARY="${HYDRO_DIR}/hydrophobicity.tsv"
CAI_SUMMARY="${CAI_DIR}/cai.tsv"

echo -e "Source_Fasta\tSequence_ID\tLength_all\tLength_no_gap" > "${LENGTH_SUMMARY}"
echo -e "Source_Fasta\tSequence_ID\tGC_percent" > "${GC_SUMMARY}"
echo -e "Source_Fasta\tSequence_ID\tIsoelectric_point" > "${PI_SUMMARY}"
echo -e "Source_Fasta\tSequence_ID\tHydrophobicity" > "${HYDRO_SUMMARY}"
echo -e "Source_Fasta\tSequence_ID\tCAI" > "${CAI_SUMMARY}"

python3 "${PYTHON_SCRIPT}" \
    --fasta-list "${FASTA_LIST_TXT}" \
    --length-out "${LENGTH_SUMMARY}" \
    --gc-out "${GC_SUMMARY}" \
    --pi-out "${PI_SUMMARY}" \
    --hydro-out "${HYDRO_SUMMARY}" \
    --cai-out "${CAI_SUMMARY}" \
    --progress-width "${PROGRESS_WIDTH}" \
    --cai-reference-fasta "${CAI_REFERENCE_FASTA}" \
    --cai-reference-list "${CAI_REFERENCE_LIST}"

echo "结果输出到以下文件："
echo "- ${LENGTH_SUMMARY}"
echo "- ${GC_SUMMARY}"
echo "- ${PI_SUMMARY}"
echo "- ${HYDRO_SUMMARY}"
echo "- ${CAI_SUMMARY}"

echo "Done!"
