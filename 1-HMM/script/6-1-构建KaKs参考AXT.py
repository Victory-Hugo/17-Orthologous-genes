#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
根据密码子比对结果生成参考序列专用的 AXT 文件（KaKs_Calculator 输入）。
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="密码子对齐 -> 参考序列 AXT")
    parser.add_argument("--alignment-dir", required=True, help="5-2 脚本输出目录（*.codon.fna）")
    parser.add_argument("--output-dir", required=True, help="AXT 输出目录")
    parser.add_argument("--reference", default=None, help="参考序列名称（与 FASTA header 完全一致）")
    parser.add_argument("--queries", nargs="*", default=None, help="限定处理的 query 名称")
    return parser.parse_args()


def read_fasta(path: Path) -> "OrderedDict[str, str]":
    seqs: "OrderedDict[str, str]" = OrderedDict()
    header = None
    seq_lines: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(seq_lines)
                header = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header is not None:
            seqs[header] = "".join(seq_lines)
    return seqs


def choose_reference(names: List[str], preferred: Optional[str]) -> str:
    if preferred and preferred in names:
        return preferred
    if preferred and preferred not in names:
        print(f"⚠️  指定的参考 {preferred} 不存在，自动改用 {names[0]}")
    return names[0]


def main() -> None:
    args = parse_args()
    alignment_dir = Path(args.alignment_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    query_filter = set(args.queries) if args.queries else None

    summary_rows: List[Dict[str, str]] = []

    for codon_file in sorted(alignment_dir.glob("*.codon.fna")):
        query_name = codon_file.name.replace(".codon.fna", "")
        if query_filter and query_name not in query_filter:
            continue
        sequences = read_fasta(codon_file)
        if len(sequences) < 2:
            print(f"⚠️  {query_name} 仅有 {len(sequences)} 条序列，跳过")
            continue
        names = list(sequences.keys())
        reference = choose_reference(names, args.reference)
        output_path = output_dir / f"{query_name}.ref_{reference}.axt"
        block_idx = 0
        with output_path.open("w", encoding="utf-8") as handle:
            for name, seq in sequences.items():
                if name == reference:
                    continue
                block_idx += 1
                handle.write(f"{block_idx} {reference} {name}\n")
                handle.write(sequences[reference] + "\n")
                handle.write(seq + "\n\n")
                summary_rows.append(
                    {
                        "query": query_name,
                        "reference": reference,
                        "target": name,
                        "blocks": block_idx,
                        "alignment_path": str(output_path),
                    }
                )
        print(f"✓ {query_name} -> 生成 {block_idx} 个 AXT 区块 ({output_path.name})")

    if not summary_rows:
        raise SystemExit("未生成任何 AXT。")

    summary_path = output_dir / "axt_summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["query", "reference", "target", "blocks", "alignment_path"],
            delimiter="\t",
        )
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)

    print(f"AXT 输出目录: {output_dir}")
    print(f"参考序列: {args.reference or summary_rows[0]['reference']}")


if __name__ == "__main__":
    main()
