#!/usr/bin/env python3
"""
Query tool for the built index with single and batch modes.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from build_index import (
    bucket_for_assembly,
    bucket_for_nucleotide,
    bucket_for_protein,
    clean_alias,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract sequences from pre-built index.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--assembly", help="Assembly accession (e.g., GCA_XXXX).")
    parser.add_argument("--query", help="Alias or identifier to search for.")
    parser.add_argument(
        "--format",
        choices=["faa", "fna", "both"],
        default="both",
        help="Output format; default prints both.",
    )
    parser.add_argument("--out", help="Output file for single-query mode; default stdout.")
    parser.add_argument(
        "--batch",
        help="CSV/TSV file with columns Assembly,Query,Format[,OutFile]. Header optional.",
    )
    parser.add_argument("--index", default="index", help="Index directory.")
    return parser.parse_args()


def load_json(path: Path) -> Dict:
    with path.open() as fh:
        return json.load(fh)


def read_sequence(path: Path) -> str:
    with path.open() as fh:
        data = json.load(fh)
        return data["sequence"]


def load_indices(index_dir: Path, assembly: str) -> Tuple[Dict, Dict]:
    assembly_bucket = bucket_for_assembly(assembly)
    alias_path = index_dir / "alias_index" / assembly_bucket / f"{assembly}.json"
    gene_path = index_dir / "gene_index" / assembly_bucket / f"{assembly}.json"
    if not alias_path.exists() or not gene_path.exists():
        raise SystemExit(f"No index files for assembly {assembly}")
    return load_json(alias_path), load_json(gene_path)


def resolve_gene(alias_index: Dict, gene_index: Dict, assembly: str, query: str) -> Tuple[str, Dict]:
    cleaned = clean_alias(query)
    alias_key = f"{assembly}::{cleaned}"
    if alias_key not in alias_index:
        raise SystemExit(f"No match for alias {query!r} in assembly {assembly}.")
    mapped = alias_index[alias_key]
    if isinstance(mapped, list):
        if len(mapped) == 1:
            gene_key = mapped[0]
        else:
            options = ", ".join(mapped[:5])
            raise SystemExit(
                f"Alias {query!r} in {assembly} maps to multiple genes: {options}"
            )
    else:
        gene_key = mapped
    return gene_key, gene_index[gene_key]


def fetch_sequences(
    index_dir: Path, assembly: str, query: str, fmt: str, entry: Dict
) -> List[str]:
    prot_id = entry["protein_id"]
    nuc_id = entry["nucleotide_id"]
    seqid, start, end, strand = entry["coords"]

    lines: List[str] = []
    if fmt in ("faa", "both"):
        prot_path = (
            index_dir
            / "sequence_store"
            / "proteins"
            / bucket_for_protein(prot_id)
            / f"{prot_id}.json"
        )
        protein_seq = read_sequence(prot_path)
        header = f">{assembly}|{query}|faa|{seqid}:{start}-{end}({strand})"
        lines.extend([header, protein_seq])
    if fmt in ("fna", "both"):
        nuc_path = (
            index_dir
            / "sequence_store"
            / "nucleotides"
            / bucket_for_nucleotide(nuc_id)
            / f"{nuc_id}.json"
        )
        nucleotide_seq = read_sequence(nuc_path)
        header = f">{assembly}|{query}|fna|{seqid}:{start}-{end}({strand})"
        lines.extend([header, nucleotide_seq])
    return lines


def main() -> None:
    args = parse_args()
    index_dir = Path(args.index)

    if args.batch:
        run_batch(args, index_dir)
        return

    if not args.assembly or not args.query:
        raise SystemExit("Single mode requires --assembly and --query")

    alias_index, gene_index = load_indices(index_dir, args.assembly)
    gene_key, entry = resolve_gene(alias_index, gene_index, args.assembly, args.query)
    lines = fetch_sequences(index_dir, args.assembly, args.query, args.format, entry)

    if args.out:
        Path(args.out).write_text("\n".join(lines) + "\n")
    else:
        sys.stdout.write("\n".join(lines) + "\n")


def run_batch(args: argparse.Namespace, index_dir: Path) -> None:
    batch_path = Path(args.batch)
    if not batch_path.exists():
        raise SystemExit(f"Batch file not found: {batch_path}")

    rows = list(read_table(batch_path))
    if not rows:
        print("[warn] batch file is empty", file=sys.stderr)
        return

    index_cache: Dict[str, Tuple[Dict, Dict]] = {}
    # 如果使用 OutFile，先清空同名文件，避免多次运行累积
    outfiles_to_clear = {r[3] for r in rows if r[3]}
    for of in outfiles_to_clear:
        try:
            Path(of).unlink()
        except FileNotFoundError:
            pass

    for asm, query, fmt, outfile in rows:
        fmt = fmt.lower() if fmt else args.format
        if fmt not in ("faa", "fna", "both"):
            print(f"[warn] skip {asm},{query}: invalid format {fmt}", file=sys.stderr)
            continue

        if asm not in index_cache:
            try:
                index_cache[asm] = load_indices(index_dir, asm)
            except SystemExit as e:
                print(f"[warn] {e}", file=sys.stderr)
                continue
        alias_index, gene_index = index_cache[asm]
        try:
            gene_key, entry = resolve_gene(alias_index, gene_index, asm, query)
            lines = fetch_sequences(index_dir, asm, query, fmt, entry)
        except SystemExit as e:
            print(f"[warn] {e}", file=sys.stderr)
            continue

        if outfile:
            out_path = Path(outfile)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            with out_path.open("a") as fh:
                fh.write("\n".join(lines) + "\n")
        else:
            # per-row file naming: assembly_query.<ext>
            targets = []
            if fmt == "both":
                targets = [("faa", lines_for_type(lines, "faa")), ("fna", lines_for_type(lines, "fna"))]
            else:
                targets = [(fmt, lines)]
            for ftype, content in targets:
                name = f"{asm}_{query}.{ftype}"
                Path(name).write_text("\n".join(content) + "\n")


def lines_for_type(lines: List[str], ftype: str) -> List[str]:
    result: List[str] = []
    for i in range(0, len(lines), 2):
        header = lines[i]
        if f"|{ftype}|" in header:
            result.extend(lines[i : i + 2])
    return result


def read_table(path: Path) -> Iterable[Tuple[str, str, str, str]]:
    text = path.read_text().splitlines()
    if not text:
        return []
    delimiter = "\t" if ("\t" in text[0]) else ","
    reader = csv.reader(text, delimiter=delimiter)
    for row in reader:
        if not row:
            continue
        if row[0].strip().lower() == "assembly":
            continue
        assembly = row[0].strip()
        query = row[1].strip() if len(row) > 1 else ""
        fmt = row[2].strip() if len(row) > 2 else ""
        outfile = row[3].strip() if len(row) > 3 else ""
        if not assembly or not query:
            continue
        yield assembly, query, fmt, outfile


if __name__ == "__main__":
    main()
