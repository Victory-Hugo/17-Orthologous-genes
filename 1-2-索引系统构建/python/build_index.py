#!/usr/bin/env python3
"""
Index builder for genome assemblies.

Reads assemblies under the configured database directory and produces:
- alias_index.json
- gene_index.json
- sequence_store/{proteins,nucleotides}/<assembly>/

No external dependencies are required.
"""
from __future__ import annotations

import argparse
import os
import hashlib
import json
import re
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build index for assemblies (FAA/FNA/GFF).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--assemblies",
        default="database",
        help="Directory containing assemblies (each in its own folder).",
    )
    parser.add_argument(
        "--index",
        default="index",
        help="Output directory for generated index.",
    )
    parser.add_argument(
        "--config",
        default=None,
        help="Optional JSON config with assemblies_dir and index_dir.",
    )
    parser.add_argument(
        "--log-skipped",
        action="store_true",
        help="Print reasons for skipped CDS entries (missing protein/coords).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Delete and rebuild the entire index directory.",
    )
    parser.add_argument(
        "--overwrite-existing",
        action="store_true",
        help="Rebuild assemblies that already have alias/gene index files.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of worker processes to use (<=0 to use all available cores).",
    )
    return parser.parse_args()


def read_fasta(path: Path) -> Iterable[Tuple[str, str]]:
    header = None
    chunks: List[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield header, "".join(chunks)


def parse_gff_attributes(attr: str) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {}
    for part in attr.split(";"):
        if not part.strip():
            continue
        if "=" in part:
            key, value = part.split("=", 1)
        else:
            key, value = part, ""
        values = [v for v in value.split(",") if v]
        out[key] = values or [value]
    return out


def clean_alias(raw: str) -> str:
    text = raw.strip().lower()
    text = text.replace('"', "").replace("'", "")
    if "=" in text:
        text = text.split("=", 1)[1]
    text = text.replace("-", "_")
    text = re.sub(r"\s+", "", text)
    text = re.sub(r"^[_:]+", "", text)
    text = re.sub(r"[^a-z0-9_]+", "_", text)
    text = re.sub(r"_+", "_", text)
    return text


def md5_12(data: str) -> str:
    return hashlib.md5(data.encode()).hexdigest()[:12]


def bucket_for_assembly(assembly: str) -> str:
    root = assembly.split(".")[0]
    return root[:6]


def bucket_for_protein(protein_id: str) -> str:
    core = protein_id.split("_", 1)[-1]
    return core[:2]


def bucket_for_nucleotide(nuc_id: str) -> str:
    core = nuc_id.split("_", 1)[-1]
    return core[:2]


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(table)[::-1]


def collect_aliases_from_header(header: str) -> List[str]:
    aliases: List[str] = []
    main = header.split()[0]
    aliases.append(main)
    for token in re.split(r"[\s|]+", header):
        token = token.strip("[]")
        if not token or "=" in token:
            continue
        aliases.append(token)
    match = re.findall(r"\[([a-zA-Z0-9_]+)=([^\]]+)\]", header)
    for key, value in match:
        aliases.append(value)
    return aliases


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def load_genome_fasta(path: Path) -> Dict[str, str]:
    genome = {}
    for header, seq in read_fasta(path):
        seqid = header.split()[0]
        genome[seqid] = seq.upper()
    return genome


def load_faa(path: Path) -> Dict[str, Tuple[str, List[str]]]:
    records: Dict[str, Tuple[str, List[str]]] = {}
    for header, seq in read_fasta(path):
        aliases = collect_aliases_from_header(header)
        key = aliases[0]
        records[key] = (seq, aliases)
    return records


def load_cds_fna(path: Path) -> Dict[str, Tuple[str, List[str]]]:
    records: Dict[str, Tuple[str, List[str]]] = {}
    pattern = re.compile(r"protein_id=([^\s\]]+)")
    for header, seq in read_fasta(path):
        matches = pattern.findall(header)
        key = matches[0] if matches else header.split()[0]
        aliases = collect_aliases_from_header(header)
        records[key] = (seq.upper(), aliases)
    return records


def build_for_assembly(
    assembly_dir: Path,
    index_dir: Path,
    log_skipped: bool,
) -> None:
    assembly = assembly_dir.name
    faa_path = assembly_dir / f"{assembly}.faa"
    fasta_path = assembly_dir / f"{assembly}.fasta"
    gff_path = assembly_dir / f"{assembly}.gff"
    cds_path = assembly_dir / f"{assembly}_CDS.fna"

    genome = load_genome_fasta(fasta_path)
    faa_records = load_faa(faa_path) if faa_path.exists() else {}
    cds_records = load_cds_fna(cds_path) if cds_path.exists() else {}

    prot_store = index_dir / "sequence_store" / "proteins"
    nuc_store = index_dir / "sequence_store" / "nucleotides"
    ensure_dir(prot_store)
    ensure_dir(nuc_store)

    alias_index: Dict[str, List[str]] = {}
    gene_index: Dict[str, Dict[str, object]] = {}

    gene_counter = 1
    protein_lookup = {k: v for k, v in faa_records.items()}
    cds_lookup = {k: v for k, v in cds_records.items()}

    skipped_no_protein = 0
    skipped_bad_coords = 0
    mismatched_cds = 0

    # 先按 protein_id 聚合多段 CDS，避免有 ribosomal slippage/多段 CDS 被拆成多个基因
    cds_features: Dict[str, Dict[str, object]] = {}
    with gff_path.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) != 9:
                continue
            seqid, _, feature_type, start, end, _, strand, phase, attrs = parts
            if feature_type != "CDS":
                continue
            attr_map = parse_gff_attributes(attrs)
            protein_id = None
            for key in ("protein_id", "Name", "ID"):
                if key in attr_map:
                    protein_id = attr_map[key][0]
                    break
            if not protein_id:
                continue
            entry = cds_features.setdefault(
                protein_id,
                {
                    "seqid": seqid,
                    "strand": strand,
                    "segments": [],
                    "aliases": [],
                },
            )
            # 如果同一个 protein_id 在不同 contig/链上，直接跳过，避免拼错
            if entry["seqid"] != seqid or entry["strand"] != strand:
                skipped_bad_coords += 1
                if log_skipped:
                    print(
                        f"[skip] {assembly} protein_id={protein_id} spans multiple contigs/strands"
                    )
                cds_features.pop(protein_id, None)
                continue
            entry["segments"].append((int(start), int(end), int(phase) if phase.isdigit() else 0))
            for vals in attr_map.values():
                entry["aliases"].extend(vals)

    def extract_segment(seqid: str, start_i: int, end_i: int, strand: str) -> str | None:
        seq = genome.get(seqid, "")
        if not seq or start_i < 1:
            return None
        seq_len = len(seq)
        if end_i > seq_len:
            wrap_len = end_i - seq_len
            if wrap_len >= seq_len:
                return None
            subseq = seq[start_i - 1 :] + seq[:wrap_len]
        else:
            subseq = seq[start_i - 1 : end_i]
        if strand == "-":
            subseq = reverse_complement(subseq)
        return subseq

    for protein_id, feature in cds_features.items():
        seqid = feature["seqid"]
        strand = feature["strand"]
        segments = sorted(feature["segments"], key=lambda x: x[0])

        subseqs: List[str] = []
        bad_segment = False
        for start_i, end_i, phase in segments:
            start_adj = start_i + phase
            subseq = extract_segment(seqid, start_adj, end_i, strand)
            if subseq is None:
                skipped_bad_coords += 1
                bad_segment = True
                if log_skipped:
                    print(
                        f"[skip] {assembly} {seqid}:{start_i}-{end_i} missing genome sequence or bad coords"
                    )
                break
            subseqs.append(subseq)
        if bad_segment or not subseqs:
            continue

        if strand == "-":
            subseqs = list(reversed(subseqs))
        joined_seq = "".join(subseqs)

        if protein_id in cds_lookup:
            cds_seq, _ = cds_lookup[protein_id]
            if cds_seq != joined_seq:
                mismatched_cds += 1

        protein_seq = None
        protein_aliases: List[str] = []
        if protein_id in protein_lookup:
            protein_seq, protein_aliases = protein_lookup[protein_id]
        else:
            skipped_no_protein += 1
            if log_skipped:
                print(
                    f"[skip] {assembly} protein_id={protein_id} "
                    "not found in FAA; manual translation disabled"
                )
            continue

        nuc_hash = f"n_{md5_12(joined_seq)}"
        prot_hash = f"p_{md5_12(protein_seq)}"

        gene_key = f"{assembly}::cds{gene_counter:05d}"
        gene_counter += 1

        for alias in feature["aliases"] + protein_aliases:
            cleaned = clean_alias(alias)
            if not cleaned:
                continue
            alias_key = f"{assembly}::{cleaned}"
            alias_index.setdefault(alias_key, [])
            if gene_key not in alias_index[alias_key]:
                alias_index[alias_key].append(gene_key)

        gene_index[gene_key] = {
            "protein_id": prot_hash,
            "nucleotide_id": nuc_hash,
            "coords": [seqid, min(s[0] for s in segments), max(s[1] for s in segments), strand],
        }

        prot_bucket = bucket_for_protein(prot_hash)
        nuc_bucket = bucket_for_nucleotide(nuc_hash)

        prot_path = prot_store / prot_bucket / f"{prot_hash}.json"
        if not prot_path.exists():
            ensure_dir(prot_path.parent)
            with prot_path.open("w") as out:
                json.dump({"id": prot_hash, "sequence": protein_seq}, out)

        nuc_path = nuc_store / nuc_bucket / f"{nuc_hash}.json"
        if not nuc_path.exists():
            ensure_dir(nuc_path.parent)
            with nuc_path.open("w") as out:
                json.dump({"id": nuc_hash, "sequence": joined_seq}, out)

    conflict_count = sum(1 for v in alias_index.values() if len(v) > 1)

    alias_bucket = bucket_for_assembly(assembly)
    alias_dir = index_dir / "alias_index" / alias_bucket
    gene_dir = index_dir / "gene_index" / alias_bucket
    ensure_dir(alias_dir)
    ensure_dir(gene_dir)

    with (alias_dir / f"{assembly}.json").open("w") as out:
        json.dump(alias_index, out, indent=2, sort_keys=True)

    with (gene_dir / f"{assembly}.json").open("w") as out:
        json.dump(gene_index, out, indent=2, sort_keys=True)

    missing_seq = 0
    for entry in gene_index.values():
        prot_path = (
            prot_store / bucket_for_protein(entry["protein_id"]) / f"{entry['protein_id']}.json"
        )
        nuc_path = (
            nuc_store / bucket_for_nucleotide(entry["nucleotide_id"]) / f"{entry['nucleotide_id']}.json"
        )
        if not prot_path.exists() or not nuc_path.exists():
            missing_seq += 1

    print(
        f"[{assembly}] genes={gene_counter-1} skipped_no_protein={skipped_no_protein} "
        f"skipped_bad_coords={skipped_bad_coords} mismatched_cds={mismatched_cds} "
        f"alias_conflicts={conflict_count} missing_seq={missing_seq}"
    )


def build_wrapper(assembly_dir: str, index_dir: str, log_skipped: bool) -> str:
    build_for_assembly(Path(assembly_dir), Path(index_dir), log_skipped)
    return Path(assembly_dir).name


def main() -> None:
    args = parse_args()
    config_data = {}
    if args.config:
        config_path = Path(args.config)
        if config_path.exists():
            with config_path.open() as fh:
                config_data = json.load(fh)
    assemblies_dir = Path(config_data.get("assemblies_dir", args.assemblies)).resolve()
    index_dir = Path(config_data.get("index_dir", args.index)).resolve()

    if index_dir.exists() and args.force:
        shutil.rmtree(index_dir)

    ensure_dir(index_dir)
    ensure_dir(index_dir / "sequence_store" / "proteins")
    ensure_dir(index_dir / "sequence_store" / "nucleotides")
    ensure_dir(index_dir / "alias_index")
    ensure_dir(index_dir / "gene_index")

    todo: List[Path] = []
    for assembly_dir in sorted(assemblies_dir.iterdir()):
        if not assembly_dir.is_dir():
            continue
        assembly = assembly_dir.name
        bucket = bucket_for_assembly(assembly)
        alias_path = index_dir / "alias_index" / bucket / f"{assembly}.json"
        gene_path = index_dir / "gene_index" / bucket / f"{assembly}.json"
        if not args.overwrite_existing and alias_path.exists() and gene_path.exists():
            print(f"[skip] {assembly} already indexed (use --overwrite-existing to rebuild)")
            continue
        todo.append(assembly_dir.resolve())

    if not todo:
        return

    workers = args.workers if args.workers != 0 else (os.cpu_count() or 1)
    workers = max(1, workers)

    if workers == 1 or len(todo) == 1:
        for assembly_dir in todo:
            build_wrapper(str(assembly_dir), str(index_dir), args.log_skipped)
        return

    with ProcessPoolExecutor(max_workers=workers) as executor:
        future_map = {
            executor.submit(build_wrapper, str(assembly_dir), str(index_dir), args.log_skipped): assembly_dir.name
            for assembly_dir in todo
        }
        for future in as_completed(future_map):
            assembly = future_map[future]
            try:
                future.result()
            except Exception as exc:
                print(f"[error] {assembly} failed: {exc}")


if __name__ == "__main__":
    main()
