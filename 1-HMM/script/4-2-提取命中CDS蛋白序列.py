#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
根据 4-1-整理HMM命中列表.py 生成的清单，提取任意基因的蛋白与 CDS 序列。

输出：
  output_dir/
    protein/<query>.faa
    cds/<query>.ffn
    summary.tsv
"""

from __future__ import annotations

import argparse
import csv
import textwrap
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="根据命中清单提取 CDS + 蛋白序列（适用于所有 HMM 基因）",
    )
    parser.add_argument("--manifest", required=True, help="4-1 脚本输出的 TSV")
    parser.add_argument("--sample-dir", required=True, help="包含 .faa/.ffn 的目录")
    parser.add_argument("--output-dir", required=True, help="输出目录")
    parser.add_argument("--queries", nargs="*", default=None, help="限定处理的 query 名称")
    return parser.parse_args()


def read_manifest(path: Path, queries: Optional[Iterable[str]]) -> List[Dict[str, str]]:
    records: List[Dict[str, str]] = []
    query_filter = set(queries) if queries else None
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if query_filter and row["query"] not in query_filter:
                continue
            records.append(row)
    return records


def load_fasta(path: Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    if not path.exists():
        raise FileNotFoundError(f"FASTA 不存在: {path}")
    header = None
    seq_lines: List[str] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(seq_lines).replace(" ", "")
                header = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header is not None:
            seqs[header] = "".join(seq_lines).replace(" ", "")
    return seqs


def format_fasta(records: List[Tuple[str, str]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        for header, seq in records:
            handle.write(f">{header}\n")
            for chunk in textwrap.wrap(seq, width=60):
                handle.write(chunk + "\n")


def main() -> None:
    args = parse_args()
    manifest_path = Path(args.manifest)
    sample_dir = Path(args.sample_dir)
    output_dir = Path(args.output_dir)
    if not manifest_path.exists():
        raise SystemExit(f"未找到清单: {manifest_path}")
    records = read_manifest(manifest_path, args.queries)
    if not records:
        raise SystemExit("清单为空，请确认输入。")

    protein_records: Dict[str, List[Tuple[str, str]]] = defaultdict(list)
    cds_records: Dict[str, List[Tuple[str, str]]] = defaultdict(list)
    summary_rows: List[Dict[str, object]] = []
    ffn_cache: Dict[str, Dict[str, str]] = {}

    for row in records:
        query = row["query"]
        sample = row["sample"]
        target_id = row["target_id"]
        if not target_id:
            print(f"⚠️  {sample} 缺少 target_id, 跳过")
            continue
        protein_header = f"{sample}|{target_id}|score={row['score']}|cov={row['coverage']}"
        protein_seq = (row.get("aa_sequence") or "").replace("*", "")
        if not protein_seq:
            print(f"⚠️  {sample}/{query} 缺少蛋白序列，跳过")
            continue
        protein_records[query].append((protein_header, protein_seq))

        if sample not in ffn_cache:
            # 兼容多种文件命名：sample.ffn、sample_CDS.fna、sample.fna
            candidates = [
                sample_dir / f"{sample}.ffn",
                sample_dir / f"{sample}_CDS.fna",
                sample_dir / f"{sample}.fna",
            ]
            ffn_path = next((p for p in candidates if p.exists()), None)
            if not ffn_path:
                print(f"⚠️  未找到 {sample} 的 CDS 文件（尝试: .ffn/_CDS.fna/.fna），跳过 {sample}")
                continue
            try:
                ffn_cache[sample] = load_fasta(ffn_path)
            except FileNotFoundError:
                print(f"⚠️  未找到 {ffn_path}，跳过 {sample}")
                continue
        cds_seq = ffn_cache[sample].get(target_id)
        if not cds_seq:
            print(f"⚠️  {sample} 的 CDS 中未找到 {target_id}")
            continue
        cds_records[query].append((protein_header, cds_seq))
        summary_rows.append(
            {
                "query": query,
                "sample": sample,
                "target_id": target_id,
                "protein_len": len(protein_seq),
                "cds_len": len(cds_seq),
                "score": row.get("score"),
                "evalue": row.get("evalue"),
                "coverage": row.get("coverage"),
            }
        )

    for query, recs in protein_records.items():
        format_fasta(recs, output_dir / "protein" / f"{query}.faa")
    for query, recs in cds_records.items():
        format_fasta(recs, output_dir / "cds" / f"{query}.ffn")

    summary_path = output_dir / "summary.tsv"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "query",
            "sample",
            "target_id",
            "protein_len",
            "cds_len",
            "score",
            "evalue",
            "coverage",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)

    print(f"✓ 共导出蛋白序列 {sum(len(v) for v in protein_records.values())} 条")
    print(f"✓ 共导出 CDS 序列 {sum(len(v) for v in cds_records.values())} 条")
    print(f"输出目录: {output_dir}")


if __name__ == "__main__":
    main()
