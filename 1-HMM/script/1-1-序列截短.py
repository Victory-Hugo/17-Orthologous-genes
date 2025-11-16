#!/usr/bin/env python3
"""
Truncate long protein sequences into overlapping windows to satisfy hmmsearch's 100k aa limit.
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Iterable, Iterator, List, Tuple

LIST_FILE = Path("/home/luolintao/0-tmp/3-Bam_Tam/conf/remained.txt")
OUTPUT_ROOT = Path("/data_ssd3/7-luolintao-ssd/0-GTDB-Database/remained_faa_short")
SOURCE_ROOT = Path("/data_ssd3/7-luolintao-ssd/0-GTDB-Database/remained_faa")

LIMIT = 100_000
OVERLAP = 5_000
WINDOW = 95_000  # Each fragment length (<= LIMIT)
STEP = WINDOW - OVERLAP
LINE_WIDTH = 60

if STEP <= 0:
    raise SystemExit("STEP must be positive")


def read_list(file_path: Path) -> List[Path]:
    items: List[Path] = []
    with file_path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            entry = raw.strip()
            if not entry or entry.startswith("#"):
                continue
            items.append(Path(entry))
    return items


def parse_fasta(file_path: Path) -> Iterator[Tuple[str, str]]:
    header: str | None = None
    seq_parts: List[str] = []
    with file_path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)


def chunk_sequence(header: str, sequence: str) -> Iterator[Tuple[str, str]]:
    length = len(sequence)
    if length <= LIMIT:
        yield header, sequence
        return

    start = 0
    index = 1
    while start < length:
        end = min(start + WINDOW, length)
        fragment = sequence[start:end]
        chunk_header = f"{header}|chunk{index}:{start + 1}-{end}"
        if len(fragment) > LIMIT:
            raise ValueError(f"Fragment still exceeds limit: {len(fragment)}")
        yield chunk_header, fragment
        if end == length:
            break
        start += STEP
        index += 1


def wrap_sequence(sequence: str, width: int = LINE_WIDTH) -> Iterable[str]:
    for idx in range(0, len(sequence), width):
        yield sequence[idx : idx + width]


def ensure_parent(target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)


def derive_output_path(source: Path) -> Path:
    try:
        relative = source.relative_to(SOURCE_ROOT)
    except ValueError:
        relative = Path(source.name)
    return OUTPUT_ROOT / relative


def process_file(source: Path) -> Tuple[int, int, int]:
    output_path = derive_output_path(source)
    ensure_parent(output_path)

    total = 0
    truncated = 0
    fragments = 0

    with output_path.open("w", encoding="utf-8") as out_handle:
        for header, sequence in parse_fasta(source):
            total += 1
            long_flag = len(sequence) > LIMIT
            fragment_count = 0
            for fragment_header, fragment_seq in chunk_sequence(header, sequence):
                out_handle.write(f">{fragment_header}\n")
                line_written = False
                for line in wrap_sequence(fragment_seq):
                    out_handle.write(f"{line}\n")
                    line_written = True
                if not line_written:
                    out_handle.write("\n")
                fragment_count += 1
            if long_flag:
                truncated += 1
            fragments += fragment_count

    return total, truncated, fragments


def main() -> None:
    sources = read_list(LIST_FILE)
    if not sources:
        print(f"No input files found in {LIST_FILE}", file=sys.stderr)
        sys.exit(1)

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    summary: List[str] = []
    processed = 0
    missing = 0
    for idx, source in enumerate(sources, start=1):
        if not source.exists():
            missing += 1
            note = f"[{idx}/{len(sources)}] missing: {source}"
            summary.append(note)
            print(note, file=sys.stderr)
            continue

        total, truncated, fragments = process_file(source)
        processed += 1
        delta = fragments - total
        delta_repr = f"+{delta}" if delta > 0 else "0"
        note = (
            f"[{idx}/{len(sources)}] {source.name}: total={total}, truncated={truncated}, fragments={fragments} ({delta_repr})"
        )
        summary.append(note)
        print(note)

    report_path = OUTPUT_ROOT / "truncation_report.txt"
    report_path.write_text("\n".join(summary) + "\n", encoding="utf-8")

    print(
        f"Completed: processed={processed}, missing={missing}. Report saved to {report_path}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
