#!/usr/bin/env python3
"""
Sequence converter using pre-built index (no translation).

Given a single-sequence FASTA (*.faa or *.fna) named with the assembly accession
(e.g., GCF_000006765.1.faa), locate the paired CDS/protein via the index and
write the counterpart sequence. Supports single query or batch CSV/TSV.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert FAA<->FNA using pre-built index (no translation).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", help="Input FASTA (filename must be assembly accession).")
    parser.add_argument("--output", help="Output path for converted sequence.")
    parser.add_argument("--batch", help="CSV/TSV with columns input,output.")
    parser.add_argument(
        "--config",
        default=str(Path(__file__).with_name("config.json")),
        help="JSON config with index_dir (and assemblies_dir, unused).",
    )
    return parser.parse_args()


def load_config(path: Path) -> Dict[str, str]:
    if path.exists():
        with path.open() as fh:
            return json.load(fh)
    return {"index_dir": "index"}


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    header = None
    seq_parts: List[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts).replace(" ", "").upper()))
                header = line[1:].strip().split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq_parts).replace(" ", "").upper()))
    return records


def md5_12(seq: str) -> str:
    return hashlib.md5(seq.encode()).hexdigest()[:12]


def bucket_for_hash(h: str) -> str:
    return h.split("_", 1)[-1][:2]


def bucket_for_assembly(assembly: str) -> str:
    return assembly.split(".")[0][:6]


def load_indices(index_dir: Path, assembly: str) -> Tuple[Dict, Dict]:
    bucket = bucket_for_assembly(assembly)
    gene_path = index_dir / "gene_index" / bucket / f"{assembly}.json"
    alias_path = index_dir / "alias_index" / bucket / f"{assembly}.json"
    if not gene_path.exists() or not alias_path.exists():
        raise SystemExit(f"No index files for assembly {assembly} at {gene_path}")
    with gene_path.open() as fh:
        gene_index = json.load(fh)
    return gene_index


def build_maps(gene_index: Dict[str, Dict]) -> Tuple[Dict[str, Tuple[str, Dict]], Dict[str, Tuple[str, Dict]]]:
    pmap: Dict[str, Tuple[str, Dict]] = {}
    nmap: Dict[str, Tuple[str, Dict]] = {}
    for gene_key, entry in gene_index.items():
        pmap[entry["protein_id"]] = (gene_key, entry)
        nmap[entry["nucleotide_id"]] = (gene_key, entry)
    return pmap, nmap


def load_seq_by_id(index_dir: Path, seq_id: str, is_protein: bool) -> str:
    store = "proteins" if is_protein else "nucleotides"
    bucket = bucket_for_hash(seq_id)
    seq_path = index_dir / "sequence_store" / store / bucket / f"{seq_id}.json"
    if not seq_path.exists():
        raise SystemExit(f"Sequence file missing: {seq_path}")
    with seq_path.open() as fh:
        return json.load(fh)["sequence"]


def write_fasta(path: Path, header: str, seq: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write(header + "\n")
        fh.write(seq + "\n")


def assembly_from_filename(path: Path) -> str:
    stem = path.stem
    # handle double extension like .fasta.gz? Not required, assume .faa/.fna
    if stem.endswith(".fasta"):
        stem = stem[:-6]
    if stem.endswith(".fa"):
        stem = stem[:-3]
    if not (stem.startswith("GCA_") or stem.startswith("GCF_")):
        raise SystemExit(f"Input filename must start with assembly accession (GCA_/GCF_): {path.name}")
    return stem


def convert_one(index_dir: Path, input_path: Path, output_path: Path) -> None:
    if not input_path.exists():
        raise SystemExit(f"Input file not found: {input_path}")
    assembly = assembly_from_filename(input_path)
    records = read_fasta(input_path)
    if not records:
        raise SystemExit(f"Input file has no sequence: {input_path}")
    is_protein = input_path.suffix.lower() in {".faa", ".aa"}
    is_nuc = input_path.suffix.lower() in {".fna", ".fa", ".fasta", ".fsa"}
    if not (is_protein or is_nuc):
        raise SystemExit(f"Unsupported extension (expect .faa/.fna): {input_path}")

    gene_index = load_indices(index_dir, assembly)
    pmap, nmap = build_maps(gene_index)

    lines: List[str] = []
    for hdr, seq in records:
        if not seq:
            continue
        if is_protein:
            pid = f"p_{md5_12(seq)}"
            if pid not in pmap:
                raise SystemExit(f"No protein hash match in index: {pid}")
            gene_key, entry = pmap[pid]
            target_id = entry["nucleotide_id"]
            target_seq = load_seq_by_id(index_dir, target_id, is_protein=False)
            seqid, start, end, strand = entry["coords"]
            header = f">{assembly}|{gene_key}|fna|{seqid}:{start}-{end}({strand})"
        else:
            nid = f"n_{md5_12(seq)}"
            if nid not in nmap:
                raise SystemExit(f"No nucleotide hash match in index: {nid}")
            gene_key, entry = nmap[nid]
            target_id = entry["protein_id"]
            target_seq = load_seq_by_id(index_dir, target_id, is_protein=True)
            seqid, start, end, strand = entry["coords"]
            header = f">{assembly}|{gene_key}|faa|{seqid}:{start}-{end}({strand})"
        lines.extend([header, target_seq])

    if not lines:
        raise SystemExit(f"No valid sequences converted from {input_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as fh:
        for i in range(0, len(lines), 2):
            fh.write(lines[i] + "\n")
            fh.write(lines[i + 1] + "\n")


def read_batch(path: Path) -> Iterable[Tuple[Path, Path]]:
    text = path.read_text().splitlines()
    if not text:
        return []
    delimiter = "\t" if ("\t" in text[0]) else ","
    reader = csv.reader(text, delimiter=delimiter)
    for row in reader:
        if not row or len(row) < 2:
            continue
        yield Path(row[0].strip()), Path(row[1].strip())


def main() -> None:
    args = parse_args()
    cfg = load_config(Path(args.config))
    index_dir = Path(cfg.get("index_dir", "index"))

    if args.batch:
        batch_path = Path(args.batch)
        if not batch_path.exists():
            raise SystemExit(f"Batch file not found: {batch_path}")
        for inp, outp in read_batch(batch_path):
            convert_one(index_dir, inp, outp)
    else:
        if not args.input or not args.output:
            raise SystemExit("Single mode requires --input and --output")
        convert_one(index_dir, Path(args.input), Path(args.output))


if __name__ == "__main__":
    main()
