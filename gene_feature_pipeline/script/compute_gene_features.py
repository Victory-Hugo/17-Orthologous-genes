"""Compute gene and protein feature metrics from nucleotide FASTA files."""
from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


CODON_TABLE: Dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


HYDROPATHY_INDEX = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5,
    "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8,
    "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5,
    "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}

PKA_VALUES = {
    "N_TERM": 9.69,
    "C_TERM": 2.34,
    "C": 8.33,
    "D": 3.86,
    "E": 4.25,
    "H": 6.0,
    "K": 10.53,
    "R": 12.48,
    "Y": 10.07,
}

CHARGE_GROUPS = {
    "positive": {
        "K": 1, "R": 1, "H": 1
    },
    "negative": {
        "D": -1, "E": -1, "C": -1, "Y": -1
    },
}

STOP_CODONS = {"TAA", "TAG", "TGA"}


@dataclass
class GeneRecord:
    identifier: str
    description: str
    sequence: str

    @property
    def cleaned_sequence(self) -> str:
        return "".join(base for base in self.sequence.upper() if base in {"A", "T", "G", "C"})


@dataclass
class GeneFeature:
    gene_id: str
    description: str
    length_nt: int
    gc_content: float
    protein_length: int
    is_valid_orf: bool
    pi: float | None
    gravy: float | None
    cai: float | None


def parse_fasta(path: Path) -> Iterable[GeneRecord]:
    with path.open() as handle:
        identifier = None
        description = ""
        seq_chunks: List[str] = []
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    yield GeneRecord(identifier, description, "".join(seq_chunks))
                header = line[1:]
                parts = header.split(None, 1)
                identifier = parts[0]
                description = parts[1] if len(parts) > 1 else ""
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if identifier is not None:
            yield GeneRecord(identifier, description, "".join(seq_chunks))


def gc_content(sequence: str) -> float:
    cleaned = [base for base in sequence.upper() if base in {"A", "T", "G", "C"}]
    if not cleaned:
        return 0.0
    gc_count = sum(1 for base in cleaned if base in {"G", "C"})
    return gc_count / len(cleaned)


def translate(sequence: str) -> Tuple[str, bool]:
    seq = sequence.upper()
    protein: List[str] = []
    is_valid = True
    if len(seq) % 3 != 0:
        is_valid = False
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if len(codon) < 3:
            is_valid = False
            break
        aa = CODON_TABLE.get(codon, "X")
        if aa == "*":
            if i < len(seq) - 3:
                # internal stop
                is_valid = False
            break
        if aa == "X":
            is_valid = False
        protein.append(aa)
    return "".join(protein), is_valid


def gravy(protein: str) -> float | None:
    if not protein:
        return None
    values = [HYDROPATHY_INDEX.get(aa) for aa in protein]
    if any(value is None for value in values):
        return None
    return sum(values) / len(values)


def net_charge_at_ph(protein: str, ph: float) -> float:
    n_term = 1 / (1 + 10 ** (ph - PKA_VALUES["N_TERM"]))
    c_term = -1 / (1 + 10 ** (PKA_VALUES["C_TERM"] - ph))
    charge = n_term + c_term
    for aa in protein:
        if aa in CHARGE_GROUPS["positive"]:
            pka = PKA_VALUES[aa]
            charge += 1 / (1 + 10 ** (ph - pka))
        elif aa in CHARGE_GROUPS["negative"]:
            pka = PKA_VALUES[aa]
            charge += -1 / (1 + 10 ** (pka - ph))
    return charge


def isoelectric_point(protein: str, precision: float = 0.01) -> float | None:
    if not protein:
        return None
    low, high = 0.0, 14.0
    while high - low > precision:
        mid = (low + high) / 2
        charge = net_charge_at_ph(protein, mid)
        if charge > 0:
            low = mid
        else:
            high = mid
    return (low + high) / 2


def load_codon_usage(path: Path) -> Dict[str, float]:
    with path.open() as handle:
        usage = json.load(handle)
    normalized = {codon.upper(): float(freq) for codon, freq in usage.items()}
    return normalized


def calculate_cai(sequence: str, usage: Dict[str, float]) -> float | None:
    seq = sequence.upper()
    weights: List[float] = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if codon in STOP_CODONS:
            break
        aa = CODON_TABLE.get(codon)
        if aa is None:
            return None
        synonymous = [c for c, translated in CODON_TABLE.items() if translated == aa and c in usage]
        if not synonymous:
            return None
        max_freq = max(usage[c] for c in synonymous)
        if max_freq == 0:
            return None
        weights.append(usage.get(codon, 0.0) / max_freq)
    if not weights:
        return None
    log_sum = sum(math.log(w if w > 0 else 1e-9) for w in weights)
    return math.exp(log_sum / len(weights))


def summarize_features(records: Iterable[GeneRecord], codon_usage: Dict[str, float] | None) -> List[GeneFeature]:
    features: List[GeneFeature] = []
    for record in records:
        sequence = record.cleaned_sequence
        length_nt = len(sequence)
        gc = gc_content(sequence)
        protein, valid_orf = translate(sequence)
        pi = isoelectric_point(protein)
        hydropathy = gravy(protein)
        cai = calculate_cai(sequence, codon_usage) if codon_usage else None
        features.append(
            GeneFeature(
                gene_id=record.identifier,
                description=record.description,
                length_nt=length_nt,
                gc_content=gc,
                protein_length=len(protein),
                is_valid_orf=valid_orf,
                pi=pi,
                gravy=hydropathy,
                cai=cai,
            )
        )
    return features


def write_csv(features: Iterable[GeneFeature], output_path: Path) -> None:
    fieldnames = [
        "gene_id",
        "description",
        "length_nt",
        "gc_content",
        "protein_length",
        "is_valid_orf",
        "pi",
        "gravy",
        "cai",
    ]
    with output_path.open("w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for feature in features:
            writer.writerow({
                "gene_id": feature.gene_id,
                "description": feature.description,
                "length_nt": feature.length_nt,
                "gc_content": f"{feature.gc_content:.4f}",
                "protein_length": feature.protein_length,
                "is_valid_orf": feature.is_valid_orf,
                "pi": f"{feature.pi:.2f}" if feature.pi is not None else "",
                "gravy": f"{feature.gravy:.3f}" if feature.gravy is not None else "",
                "cai": f"{feature.cai:.3f}" if feature.cai is not None else "",
            })


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compute nucleotide and protein features for genes.")
    parser.add_argument("fasta", type=Path, help="Input FASTA file containing coding sequences.")
    parser.add_argument("--codon-usage", type=Path, help="JSON file with codon usage frequencies.")
    parser.add_argument("--output", type=Path, default=Path("gene_features.csv"), help="Output CSV path.")
    return parser


def main(args: List[str] | None = None) -> None:
    parser = build_argument_parser()
    namespace = parser.parse_args(args=args)
    records = list(parse_fasta(namespace.fasta))
    codon_usage = load_codon_usage(namespace.codon_usage) if namespace.codon_usage else None
    features = summarize_features(records, codon_usage)
    write_csv(features, namespace.output)


if __name__ == "__main__":
    main()
