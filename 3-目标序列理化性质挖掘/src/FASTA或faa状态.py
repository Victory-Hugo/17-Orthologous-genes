#!/usr/bin/env python3
"""
Compute sequence statistics (length, GC, pI, hydrophobicity, CAI) for each record in
multiple FASTA files. This module is invoked by the companion shell script.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

DNA_ALPHABET = set("ACGTRYKMSWBDHVN")
DNA_STRICT = set("ACGT")
PROTEIN_ALLOWED = set("ACDEFGHIKLMNPQRSTVWY")

CODON_TO_AA = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

AA_TO_CODONS = {}
for codon, aa in CODON_TO_AA.items():
    AA_TO_CODONS.setdefault(aa, []).append(codon)


class SimpleCodonAdaptationIndex:
    """Minimal CAI implementation to remove dependency on Bio.SeqUtils.CodonUsage."""

    def __init__(self):
        self.weights: dict[str, float] = {}

    def _iter_codons(self, sequence: str):
        for idx in range(0, len(sequence) - 2, 3):
            codon = sequence[idx : idx + 3]
            if len(codon) == 3 and all(base in DNA_STRICT for base in codon):
                yield codon

    def generate_index(self, sequences: Sequence[str]) -> None:
        from collections import Counter

        codon_counts = Counter()
        aa_totals = Counter()

        for seq in sequences:
            for codon in self._iter_codons(seq):
                aa = CODON_TO_AA.get(codon)
                if aa and aa != "*":
                    codon_counts[codon] += 1
                    aa_totals[aa] += 1

        weights: dict[str, float] = {}
        for aa, total in aa_totals.items():
            if total == 0:
                continue
            max_freq = 0.0
            codon_freqs = {}
            for codon in AA_TO_CODONS.get(aa, []):
                freq = codon_counts.get(codon, 0) / total
                codon_freqs[codon] = freq
                if freq > max_freq:
                    max_freq = freq
            if max_freq == 0.0:
                continue
            for codon, freq in codon_freqs.items():
                if freq > 0.0:
                    weights[codon] = freq / max_freq

        if not weights:
            raise ValueError("无法根据参考序列生成 CAI 权重。")
        self.weights = weights

    def cai_for_gene(self, sequence: str):
        import math

        logs = []
        for codon in self._iter_codons(sequence):
            weight = self.weights.get(codon)
            if not weight or weight <= 0.0:
                return None
            logs.append(math.log(weight))
        if not logs:
            return None
        return math.exp(sum(logs) / len(logs))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyse FASTA files listed in a text file and write TSV summaries."
    )
    parser.add_argument("--fasta-list", required=True, help="Text file listing FASTA paths.")
    parser.add_argument("--length-out", required=True, help="TSV file for length metrics.")
    parser.add_argument("--gc-out", required=True, help="TSV file for GC percentages.")
    parser.add_argument("--pi-out", required=True, help="TSV file for isoelectric points.")
    parser.add_argument("--hydro-out", required=True, help="TSV file for hydrophobicity.")
    parser.add_argument("--cai-out", required=True, help="TSV file for CAI values.")
    parser.add_argument("--progress-width", type=int, default=40, help="Progress bar width.")
    parser.add_argument("--cai-reference-fasta", default="", help="Reference CDS FASTA for CAI.")
    parser.add_argument(
        "--cai-reference-list",
        default="",
        help="Text file listing multiple reference FASTA files for CAI.",
    )
    return parser.parse_args()


def read_fasta_list(path: Path) -> List[Path]:
    with path.open(encoding="utf-8") as handle:
        return [Path(line.strip()) for line in handle if line.strip()]


def detect_sequence_type(seq: str) -> str:
    seq_upper = seq.upper()
    letters = [c for c in seq_upper if c.isalpha()]
    if not letters:
        return "unknown"
    dna_letters = sum(1 for c in letters if c in DNA_ALPHABET)
    return "dna" if dna_letters / len(letters) >= 0.9 else "protein"


def gc_percent(seq: str):
    seq_upper = seq.upper().replace("U", "T")
    filtered = [c for c in seq_upper if c in DNA_STRICT]
    if not filtered:
        return None
    g = filtered.count("G")
    c = filtered.count("C")
    return (g + c) * 100.0 / len(filtered)


def clean_protein(seq: str) -> str:
    seq_upper = seq.upper().replace("*", "")
    return "".join(c for c in seq_upper if c in PROTEIN_ALLOWED)


def sanitise_cds(seq: str):
    seq_upper = seq.upper().replace("U", "T")
    if not seq_upper:
        return None
    if any(c not in DNA_STRICT for c in seq_upper):
        return None
    remainder = len(seq_upper) % 3
    trimmed = seq_upper if remainder == 0 else seq_upper[: len(seq_upper) - remainder]
    if trimmed.endswith(("TAA", "TAG", "TGA")):
        trimmed = trimmed[:-3]
    return trimmed if trimmed else None


def render_progress(current: int, total: int, width: int):
    if total <= 0:
        return
    filled = current * width // total
    bar = "#" * filled + " " * (width - filled)
    percent = current * 100 // total
    print(f"\r[{bar}] {percent:3d}% ({current}/{total})", end="", file=sys.stderr, flush=True)


def collect_reference_sequences(path_str: str, missing: List[str]) -> List[str]:
    sequences: List[str] = []
    path = Path(path_str)
    if not path.is_file():
        missing.append(path_str)
        return sequences
    try:
        for record in SeqIO.parse(str(path), "fasta"):
            cds = sanitise_cds(str(record.seq))
            if cds:
                sequences.append(cds)
    except Exception as exc:  # pragma: no cover - defensive
        missing.append(f"{path_str} ({exc})")
    return sequences


def load_cai_reference(ref_fasta: str, ref_list: str) -> Tuple[List[str], List[str]]:
    sequences: List[str] = []
    missing: List[str] = []

    if ref_fasta.strip():
        sequences.extend(collect_reference_sequences(ref_fasta.strip(), missing))

    if ref_list.strip():
        list_path = Path(ref_list.strip())
        try:
            with list_path.open(encoding="utf-8") as handle:
                for line in handle:
                    entry = line.strip()
                    if entry:
                        sequences.extend(collect_reference_sequences(entry, missing))
        except FileNotFoundError:
            missing.append(f"列表文件不存在: {ref_list}")

    return sequences, missing


def init_cai_calculator(reference_sequences: Sequence[str]):
    if not reference_sequences:
        return None
    calculator = SimpleCodonAdaptationIndex()
    try:
        calculator.generate_index(reference_sequences)
    except Exception as exc:  # pragma: no cover - defensive
        print(f"\n[warning] CAI 参考序列初始化失败: {exc}", file=sys.stderr)
        return None
    return calculator


def iter_records(fasta_path: Path):
    try:
        yield from SeqIO.parse(str(fasta_path), "fasta")
    except Exception as exc:
        print(f"\n[warning] 解析失败 {fasta_path}: {exc}", file=sys.stderr)


def process_fasta_files(
    fasta_paths: Iterable[Path],
    outputs: Tuple[Path, Path, Path, Path, Path],
    progress_width: int,
    cai_reference_fasta: str,
    cai_reference_list: str,
):
    fasta_paths = list(fasta_paths)
    total_files = len(fasta_paths)
    if total_files == 0:
        print("FASTA 列表文件为空。", file=sys.stderr)
        return

    reference_sequences, missing_reference = load_cai_reference(
        cai_reference_fasta, cai_reference_list
    )
    cai_calculator = init_cai_calculator(reference_sequences)

    if not reference_sequences and (cai_reference_fasta.strip() or cai_reference_list.strip()):
        print("\n[warning] 未从提供的 CAI 参考路径中获得有效 CDS 序列，CAI 将输出 NA。", file=sys.stderr)
    if missing_reference:
        for item in missing_reference:
            print(f"\n[warning] CAI 参考缺失或读取失败: {item}", file=sys.stderr)

    length_path, gc_path, pi_path, hydro_path, cai_path = outputs

    with (
        length_path.open("a", encoding="utf-8") as length_fh,
        gc_path.open("a", encoding="utf-8") as gc_fh,
        pi_path.open("a", encoding="utf-8") as pi_fh,
        hydro_path.open("a", encoding="utf-8") as hydro_fh,
        cai_path.open("a", encoding="utf-8") as cai_fh,
    ):
        for idx, fasta_path in enumerate(fasta_paths, start=1):
            if not fasta_path.is_file():
                print(f"\n[warning] 找不到 FASTA 文件: {fasta_path}", file=sys.stderr)
                render_progress(idx, total_files, progress_width)
                continue

            records = list(iter_records(fasta_path))
            if not records:
                print(f"\n[warning] FASTA 文件中没有序列: {fasta_path}", file=sys.stderr)

            for record in records:
                raw_seq = str(record.seq)
                compact_seq = raw_seq.replace("\n", "").replace("\r", "").replace("-", "")

                length_fh.write(f"{fasta_path.name}\t{record.id}\t{len(raw_seq)}\t{len(compact_seq)}\n")

                seq_type = detect_sequence_type(compact_seq)

                gc_value = "NA"
                if seq_type == "dna":
                    gc_calc = gc_percent(compact_seq)
                    if gc_calc is not None:
                        gc_value = f"{gc_calc:.4f}"
                gc_fh.write(f"{fasta_path.name}\t{record.id}\t{gc_value}\n")

                pi_value = "NA"
                hydro_value = "NA"
                if seq_type == "protein":
                    protein_seq = clean_protein(compact_seq)
                    if protein_seq:
                        try:
                            analysis = ProteinAnalysis(protein_seq)
                            pi_value = f"{analysis.isoelectric_point():.4f}"
                            hydro_value = f"{analysis.gravy():.4f}"
                        except Exception:
                            pi_value = "NA"
                            hydro_value = "NA"
                pi_fh.write(f"{fasta_path.name}\t{record.id}\t{pi_value}\n")
                hydro_fh.write(f"{fasta_path.name}\t{record.id}\t{hydro_value}\n")

                cai_value = "NA"
                if cai_calculator is not None and seq_type == "dna":
                    cds_seq = sanitise_cds(compact_seq)
                    if cds_seq:
                        try:
                            cai_score = cai_calculator.cai_for_gene(cds_seq)
                            cai_value = f"{cai_score:.4f}"
                        except Exception:
                            cai_value = "NA"
                cai_fh.write(f"{fasta_path.name}\t{record.id}\t{cai_value}\n")

            render_progress(idx, total_files, progress_width)

    print("\n", end="", file=sys.stderr)


def main():
    args = parse_args()
    fasta_list_path = Path(args.fasta_list)
    if not fasta_list_path.is_file():
        print(f"FASTA 列表文件不存在: {fasta_list_path}", file=sys.stderr)
        sys.exit(1)

    fasta_paths = read_fasta_list(fasta_list_path)
    outputs = (
        Path(args.length_out),
        Path(args.gc_out),
        Path(args.pi_out),
        Path(args.hydro_out),
        Path(args.cai_out),
    )

    process_fasta_files(
        fasta_paths=fasta_paths,
        outputs=outputs,
        progress_width=args.progress_width,
        cai_reference_fasta=args.cai_reference_fasta,
        cai_reference_list=args.cai_reference_list,
    )


if __name__ == "__main__":
    main()
