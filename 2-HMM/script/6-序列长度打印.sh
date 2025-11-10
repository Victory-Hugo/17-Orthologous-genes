#!/bin/bash
#? 统计多个 FASTA 文件的多种理化性质指标，并将结果分别输出到指定目录中，用于分析基因在不同物种中的长度和理化特征变化。

set -euo pipefail

FASTA_LIST_TXT="/home/luolintao/0_Github/17-Orthologous-genes/2-HMM/conf/4-搜索完_list.txt"
OUTPUT_BASE="/home/luolintao/5-AB-Baoman/2-Lol家族基因搜索/output"
LENGTH_DIR="${OUTPUT_BASE}/1-序列长度"
GC_DIR="${OUTPUT_BASE}/2-GC含量"
PI_DIR="${OUTPUT_BASE}/3-等电点"
HYDRO_DIR="${OUTPUT_BASE}/4-疏水性指数"
CAI_DIR="${OUTPUT_BASE}/5-CAI"

# 配置 CAI 计算所需的宿主高表达基因 CDS 参考序列（DNA）
CAI_REFERENCE_FASTA=""   # 直接指定一个 FASTA 文件路径
CAI_REFERENCE_LIST=""    # 或者指定一个列出多个 FASTA 文件路径的文本文件

PROGRESS_WIDTH=40

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

export FASTA_LIST_TXT
export LENGTH_SUMMARY
export GC_SUMMARY
export PI_SUMMARY
export HYDRO_SUMMARY
export CAI_SUMMARY
export PROGRESS_WIDTH
export CAI_REFERENCE_FASTA
export CAI_REFERENCE_LIST

python3 <<'PYTHON'
import os
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex

fasta_list_path = Path(os.environ['FASTA_LIST_TXT'])
length_path = Path(os.environ['LENGTH_SUMMARY'])
gc_path = Path(os.environ['GC_SUMMARY'])
pi_path = Path(os.environ['PI_SUMMARY'])
hydro_path = Path(os.environ['HYDRO_SUMMARY'])
cai_path = Path(os.environ['CAI_SUMMARY'])
progress_width = int(os.environ.get('PROGRESS_WIDTH', '40'))

DNA_ALPHABET = set("ACGTRYKMSWBDHVN")
DNA_STRICT = set("ACGT")
PROTEIN_ALLOWED = set("ACDEFGHIKLMNPQRSTVWY")


def detect_sequence_type(seq: str) -> str:
    seq_upper = seq.upper()
    letters = [c for c in seq_upper if c.isalpha()]
    if not letters:
        return "unknown"
    dna_letters = sum(1 for c in letters if c in DNA_ALPHABET)
    dna_ratio = dna_letters / len(letters)
    return "dna" if dna_ratio >= 0.9 else "protein"


def gc_percent(seq: str):
    seq_upper = seq.upper().replace('U', 'T')
    filtered = [c for c in seq_upper if c in DNA_STRICT]
    if not filtered:
        return None
    g = filtered.count('G')
    c = filtered.count('C')
    return (g + c) * 100.0 / len(filtered)


def clean_protein(seq: str) -> str:
    seq_upper = seq.upper().replace('*', '')
    return ''.join(c for c in seq_upper if c in PROTEIN_ALLOWED)


def sanitise_cds(seq: str):
    seq_upper = seq.upper().replace('U', 'T')
    if not seq_upper:
        return None
    if any(c not in DNA_STRICT for c in seq_upper):
        return None
    remainder = len(seq_upper) % 3
    trimmed = seq_upper if remainder == 0 else seq_upper[:len(seq_upper) - remainder]
    if trimmed.endswith(('TAA', 'TAG', 'TGA')):
        trimmed = trimmed[:-3]
    return trimmed if trimmed else None


def load_cai_reference():
    ref_fasta = os.environ.get('CAI_REFERENCE_FASTA', '').strip()
    ref_list = os.environ.get('CAI_REFERENCE_LIST', '').strip()
    sequences = []
    missing = []

    def collect_from_file(path_str: str):
        path = Path(path_str)
        if not path.is_file():
            missing.append(path_str)
            return
        try:
            for record in SeqIO.parse(str(path), 'fasta'):
                cds = sanitise_cds(str(record.seq))
                if cds:
                    sequences.append(cds)
        except Exception as exc:  # pragma: no cover - defensive
            missing.append(f"{path_str} ({exc})")

    if ref_fasta:
        collect_from_file(ref_fasta)

    if ref_list:
        try:
            with open(ref_list, encoding='utf-8') as handle:
                for line in handle:
                    entry = line.strip()
                    if entry:
                        collect_from_file(entry)
        except FileNotFoundError:
            missing.append(f"列表文件不存在: {ref_list}")

    return sequences, missing


def render_progress(current: int, total: int):
    if total <= 0:
        return
    filled = current * progress_width // total
    bar = '#' * filled + ' ' * (progress_width - filled)
    percent = current * 100 // total
    print(f"\r[{bar}] {percent:3d}% ({current}/{total})", end='', file=sys.stderr, flush=True)


with open(fasta_list_path, encoding='utf-8') as handle:
    fasta_paths = [line.strip() for line in handle if line.strip()]

total_files = len(fasta_paths)
if total_files == 0:
    print("FASTA 列表文件为空。", file=sys.stderr)
    sys.exit(0)

reference_sequences, missing_reference = load_cai_reference()
cai_calculator = None

if reference_sequences:
    cai_calculator = CodonAdaptationIndex()
    try:
        cai_calculator.generate_index(reference_sequences)
    except Exception as exc:  # pragma: no cover - defensive
        print(f"\n[warning] CAI 参考序列初始化失败: {exc}", file=sys.stderr)
        cai_calculator = None
elif os.environ.get('CAI_REFERENCE_FASTA', '').strip() or os.environ.get('CAI_REFERENCE_LIST', '').strip():
    print("\n[warning] 未从提供的 CAI 参考路径中获得有效 CDS 序列，CAI 将输出 NA。", file=sys.stderr)

if missing_reference:
    for item in missing_reference:
        print(f"\n[warning] CAI 参考缺失或读取失败: {item}", file=sys.stderr)


with (
    open(length_path, 'a', encoding='utf-8') as length_fh,
    open(gc_path, 'a', encoding='utf-8') as gc_fh,
    open(pi_path, 'a', encoding='utf-8') as pi_fh,
    open(hydro_path, 'a', encoding='utf-8') as hydro_fh,
    open(cai_path, 'a', encoding='utf-8') as cai_fh,
):
    for idx, fasta_str in enumerate(fasta_paths, start=1):
        fasta_path = Path(fasta_str)
        if not fasta_path.is_file():
            print(f"\n[warning] 找不到 FASTA 文件: {fasta_path}", file=sys.stderr)
            render_progress(idx, total_files)
            continue

        try:
            records = list(SeqIO.parse(str(fasta_path), 'fasta'))
        except Exception as exc:
            print(f"\n[warning] 解析失败 {fasta_path}: {exc}", file=sys.stderr)
            render_progress(idx, total_files)
            continue

        if not records:
            print(f"\n[warning] FASTA 文件中没有序列: {fasta_path}", file=sys.stderr)

        for record in records:
            raw_seq = str(record.seq)
            compact_seq = raw_seq.replace('\n', '').replace('\r', '').replace('-', '')

            length_fh.write(f"{fasta_path.name}\t{record.id}\t{len(raw_seq)}\t{len(compact_seq)}\n")

            seq_type = detect_sequence_type(compact_seq)

            gc_value = 'NA'
            if seq_type == 'dna':
                gc_calc = gc_percent(compact_seq)
                if gc_calc is not None:
                    gc_value = f"{gc_calc:.4f}"
            gc_fh.write(f"{fasta_path.name}\t{record.id}\t{gc_value}\n")

            pi_value = 'NA'
            hydro_value = 'NA'
            if seq_type == 'protein':
                protein_seq = clean_protein(compact_seq)
                if protein_seq:
                    try:
                        analysis = ProteinAnalysis(protein_seq)
                        pi_value = f"{analysis.isoelectric_point():.4f}"
                        hydro_value = f"{analysis.gravy():.4f}"
                    except Exception:
                        pi_value = 'NA'
                        hydro_value = 'NA'
            pi_fh.write(f"{fasta_path.name}\t{record.id}\t{pi_value}\n")
            hydro_fh.write(f"{fasta_path.name}\t{record.id}\t{hydro_value}\n")

            cai_value = 'NA'
            if cai_calculator is not None and seq_type == 'dna':
                cds_seq = sanitise_cds(compact_seq)
                if cds_seq:
                    try:
                        cai_score = cai_calculator.cai_for_gene(cds_seq)
                        cai_value = f"{cai_score:.4f}"
                    except Exception:
                        cai_value = 'NA'
            cai_fh.write(f"{fasta_path.name}\t{record.id}\t{cai_value}\n")

        render_progress(idx, total_files)

print('\n', end='', file=sys.stderr)
PYTHON

echo "结果输出到以下文件："
echo "- ${LENGTH_SUMMARY}"
echo "- ${GC_SUMMARY}"
echo "- ${PI_SUMMARY}"
echo "- ${HYDRO_SUMMARY}"
echo "- ${CAI_SUMMARY}"

echo "Done!"
