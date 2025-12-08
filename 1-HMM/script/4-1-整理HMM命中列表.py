#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
根据 4-提取匹配序列.py 的输出，整理任意 HMM 命中列表（适用于 LolD、LolF 及其他基因）。

主要功能：
  1. 解析各 query 目录下的 *.hits.faa 文件，读取 header 中的元信息
  2. 导出标准化的 TSV 清单，方便后续提取 CDS、蛋白序列或进行统计
  3. 可选择特定 query，默认遍历全部 query

使用示例：
  python 4-1-整理HMM命中列表.py \\
      --hits-dir /path/to/output \\
      --output-tsv /path/to/pipeline/hits_manifest.tsv \\
      --queries lolD.trim
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict, OrderedDict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="整理 HMM 命中列表，输出标准化 TSV",
    )
    parser.add_argument(
        "--hits-dir",
        required=True,
        help="4-提取匹配序列.py 的输出目录（包含 query 子目录）",
    )
    parser.add_argument(
        "--output-tsv",
        required=True,
        help="输出 TSV 文件路径",
    )
    parser.add_argument(
        "--queries",
        nargs="*",
        default=None,
        help="仅处理指定 query 名（默认处理全部）",
    )
    return parser.parse_args()


def parse_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    header = None
    seq_lines: List[str] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines).replace(" ", "")
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header is not None:
            yield header, "".join(seq_lines).replace(" ", "")


def parse_metadata(header: str) -> Tuple[str, Dict[str, str]]:
    parts = header.split("|")
    record_id = parts[0]
    metadata: Dict[str, str] = OrderedDict()
    for part in parts[1:]:
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        metadata[key.strip()] = value.strip()
    return record_id, metadata


def float_or_none(value: Optional[str]) -> Optional[float]:
    if value is None or value == "" or value.lower() == "na":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def collect_hits(hits_dir: Path, queries: Optional[Iterable[str]]) -> List[Dict[str, object]]:
    result: List[Dict[str, object]] = []
    available_queries = sorted([p.name for p in hits_dir.iterdir() if p.is_dir()])
    query_filter = set(queries) if queries else None

    for query_name in available_queries:
        if query_filter and query_name not in query_filter:
            continue
        query_dir = hits_dir / query_name
        sequence_dir = query_dir / "sequence"
        candidates = []
        if sequence_dir.is_dir():
            candidates.extend(sorted(sequence_dir.glob("*.hits.faa")))
        candidates.extend(sorted(query_dir.glob("*.hits.faa")))
        seen = set()
        fasta_files = []
        for path in candidates:
            if path not in seen:
                fasta_files.append(path)
                seen.add(path)
        for fasta_file in fasta_files:
            sample_name = fasta_file.name.replace(".hits.faa", "")
            for record_idx, (header, sequence) in enumerate(parse_fasta(fasta_file), start=1):
                clean_sequence = sequence.replace("*", "")
                record_id, metadata = parse_metadata(header)
                score = float_or_none(metadata.get("score"))
                evalue = float_or_none(metadata.get("evalue"))
                coverage = float_or_none(metadata.get("cov"))
                aa_length = int(metadata.get("len") or len(clean_sequence))
                result.append(
                    {
                        "query": query_name,
                        "sample": sample_name,
                        "record": record_id,
                        "target_id": metadata.get("id", ""),
                        "score": score,
                        "evalue": evalue,
                        "coverage": coverage,
                        "aa_length": aa_length,
                        "aa_sequence": clean_sequence,
                        "source_fasta": str(fasta_file),
                        "sequence_index": record_idx,
                    }
                )
    return result


def write_tsv(records: List[Dict[str, object]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "query",
        "sample",
        "record",
        "target_id",
        "score",
        "evalue",
        "coverage",
        "aa_length",
        "aa_sequence",
        "sequence_index",
        "source_fasta",
    ]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in records:
            writer.writerow(row)


def print_summary(records: List[Dict[str, object]]) -> None:
    if not records:
        print("⚠️  未找到任何命中记录", file=sys.stderr)
        return
    print(f"✓ 收集命中 {len(records)} 条")
    per_query: Dict[str, int] = defaultdict(int)
    per_sample: Dict[str, int] = defaultdict(int)
    for rec in records:
        per_query[rec["query"]] += 1
        per_sample[rec["sample"]] += 1
    print("└─ 每个 query 的命中数：")
    for query, count in sorted(per_query.items()):
        print(f"   · {query}: {count}")
    print(f"└─ 覆盖样本数: {len(per_sample)}")


def main() -> None:
    args = parse_args()
    hits_dir = Path(args.hits_dir).resolve()
    if not hits_dir.exists():
        raise SystemExit(f"命中目录不存在: {hits_dir}")
    records = collect_hits(hits_dir, args.queries)
    write_tsv(records, Path(args.output_tsv))
    print_summary(records)


if __name__ == "__main__":
    main()
