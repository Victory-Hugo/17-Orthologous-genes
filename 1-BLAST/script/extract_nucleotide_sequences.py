#!/usr/bin/env python3
"""
Extract nucleotide sequences from FNA files based on BLASTN results (streaming version).

Usage:
    python3 extract_nucleotide_sequences.py <blastn_results> <target_path_list> <output_dir>
"""

import os
import sys
from collections import defaultdict


def resolve_fna_path(path):
    path = path.strip()
    if not path:
        return None
    if os.path.isfile(path) and path.endswith(".fna"):
        return path
    if path.endswith(".faa"):
        candidate = path[:-4] + ".fna"
    else:
        root, _ = os.path.splitext(path)
        candidate = root + ".fna"
    return candidate if os.path.isfile(candidate) else None


def build_fna_path_map(target_path_list):
    """
    Build a mapping of filename to full FNA path.
    """
    fna_path_map = {}

    with open(target_path_list, "r") as handle:
        for line in handle:
            resolved = resolve_fna_path(line)
            if resolved:
                filename = os.path.basename(resolved)
                fna_path_map[filename] = resolved

    return fna_path_map


def collect_required_sequences(blastn_results):
    """
    Collect which sequences we need from which files.
    """
    required = defaultdict(set)

    with open(blastn_results, "r") as handle:
        next(handle)
        for line in handle:
            fields = line.strip().split("\t")
            if len(fields) < 4:
                continue
            source_file = fields[0]
            sseqid = fields[2]
            required[source_file].add(sseqid)

    return required


def extract_sequence_from_file(fna_file, target_ids, output_dir):
    """
    Stream through a single FNA file and extract required sequences.
    """
    source_file = os.path.basename(fna_file)
    source_id = os.path.splitext(source_file)[0]

    extracted_count = 0
    current_id = None
    current_seq = []
    writing = False

    with open(fna_file, "r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if writing and current_id and current_seq:
                    output_file = os.path.join(output_dir, f"{source_id}.fna")
                    with open(output_file, "w") as out_handle:
                        out_handle.write(f">{source_id}|{current_id}\n")
                        seq = "".join(current_seq)
                        for i in range(0, len(seq), 60):
                            out_handle.write(seq[i:i + 60] + "\n")
                    extracted_count += 1

                current_id = line[1:].split()[0]
                current_seq = []
                writing = current_id in target_ids
            else:
                if writing:
                    current_seq.append(line)

        if writing and current_id and current_seq:
            output_file = os.path.join(output_dir, f"{source_id}.fna")
            with open(output_file, "w") as out_handle:
                out_handle.write(f">{source_id}|{current_id}\n")
                seq = "".join(current_seq)
                for i in range(0, len(seq), 60):
                    out_handle.write(seq[i:i + 60] + "\n")
            extracted_count += 1

    return extracted_count


def extract_sequences(blastn_results, target_path_list, output_dir):
    """
    Extract matching nucleotide sequences from FNA files based on BLASTN results.
    """
    os.makedirs(output_dir, exist_ok=True)

    print("Building file path mapping...")
    fna_path_map = build_fna_path_map(target_path_list)

    print("Analyzing required sequences...")
    required_seqs = collect_required_sequences(blastn_results)

    print("Extracting sequences (streaming)...")
    total_extracted = 0

    for source_file, target_ids in required_seqs.items():
        if source_file not in fna_path_map:
            print(f"Warning: File {source_file} not found in path list", file=sys.stderr)
            continue

        fna_file = fna_path_map[source_file]
        count = extract_sequence_from_file(fna_file, target_ids, output_dir)
        total_extracted += count
        print(f"  Processed {source_file}: {count} sequences extracted")

    print(f"Total: Extracted {total_extracted} sequences to {output_dir}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: extract_nucleotide_sequences.py <blastn_results> <target_path_list> <output_dir>")
        sys.exit(1)

    blastn_results = sys.argv[1]
    target_path_list = sys.argv[2]
    output_dir = sys.argv[3]

    extract_sequences(blastn_results, target_path_list, output_dir)
