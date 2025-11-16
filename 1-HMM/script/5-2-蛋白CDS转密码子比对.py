#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
调用专业工具 PAL2NAL：利用蛋白 MSA + CDS FASTA 生成密码子对齐，可选输出 all-vs-all AXT。
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from collections import OrderedDict
from pathlib import Path
from typing import List, Sequence


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="使用 PAL2NAL 生成密码子对齐")
    parser.add_argument("--alignment-dir", required=True, help="蛋白比对目录（*.aln.faa）")
    parser.add_argument("--cds-dir", required=True, help="CDS FASTA 目录（*.ffn）")
    parser.add_argument("--output-dir", required=True, help="密码子比对输出目录")
    parser.add_argument("--queries", nargs="*", default=None, help="仅处理指定 query")
    parser.add_argument("--axt-dir", default=None, help="可选：输出 all-vs-all AXT 文件目录")
    parser.add_argument("--pal2nal-bin", default="pal2nal.pl", help="PAL2NAL 可执行文件路径或命令名")
    parser.add_argument("--codon-table", type=int, default=11, help="遗传密码表 (PAL2NAL -codontable 参数)")
    parser.add_argument(
        "--pal2nal-extra",
        nargs="*",
        default=None,
        help="PAL2NAL 额外参数（按原样追加，例如 -nogap）",
    )
    return parser.parse_args()


def ensure_pal2nal(bin_name: str) -> str:
    resolved = shutil.which(bin_name)
    if resolved:

        return resolved
    candidate = Path(bin_name)
    if candidate.exists():
        return str(candidate)
    raise SystemExit(f"未找到 PAL2NAL 可执行文件: {bin_name}，请确认已安装并在 PATH 或指定绝对路径。")


def read_fasta(path: Path) -> "OrderedDict[str, str]":
    sequences: "OrderedDict[str, str]" = OrderedDict()
    header = None
    seq_lines: List[str] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    sequences[header] = "".join(seq_lines).replace(" ", "")
                header = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header is not None:
            sequences[header] = "".join(seq_lines).replace(" ", "")
    return sequences


def write_all_pairs_axt(codon_fasta: Path, output_path: Path) -> None:
    seqs = read_fasta(codon_fasta)
    names = list(seqs.keys())
    block_idx = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                block_idx += 1
                s1 = names[i]
                s2 = names[j]
                handle.write(f"{block_idx} {s1} {s2}\n")
                handle.write(seqs[s1] + "\n")
                handle.write(seqs[s2] + "\n\n")


def run_pal2nal(
    pal2nal_bin: str,
    protein_alignment: Path,
    cds_fasta: Path,
    output_file: Path,
    codon_table: int,
    extra_args: Sequence[str] | None = None,
) -> None:
    output_file.parent.mkdir(parents=True, exist_ok=True)
    cmd: List[str] = [
        pal2nal_bin,
        str(protein_alignment),
        str(cds_fasta),
        "-output",
        "fasta",
        "-codontable",
        str(codon_table),
    ]
    if extra_args:
        cmd.extend(extra_args)
    with output_file.open("w", encoding="utf-8") as out_handle:
        subprocess.run(cmd, stdout=out_handle, stderr=subprocess.PIPE, check=True, text=True)


def main() -> None:
    args = parse_args()
    pal2nal_bin = ensure_pal2nal(args.pal2nal_bin)
    alignment_dir = Path(args.alignment_dir)
    cds_dir = Path(args.cds_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    axt_dir = Path(args.axt_dir) if args.axt_dir else None
    query_filter = set(args.queries) if args.queries else None

    processed = 0
    for aln_file in sorted(alignment_dir.glob("*.aln.faa")):
        query_name = aln_file.name.replace(".aln.faa", "")
        if query_filter and query_name not in query_filter:
            continue
        cds_path = cds_dir / f"{query_name}.ffn"
        if not cds_path.exists():
            print(f"⚠️  未找到 CDS: {cds_path}，跳过 {query_name}")
            continue

        codon_output = output_dir / f"{query_name}.codon.fna"
        try:
            run_pal2nal(
                pal2nal_bin=pal2nal_bin,
                protein_alignment=aln_file,
                cds_fasta=cds_path,
                output_file=codon_output,
                codon_table=args.codon_table,
                extra_args=args.pal2nal_extra,
            )
        except subprocess.CalledProcessError as exc:
            print(f"✗ PAL2NAL 运行失败 ({query_name}): {exc.stderr.strip()}")
            continue

        if axt_dir:
            write_all_pairs_axt(codon_output, axt_dir / f"{query_name}.all_pairs.axt")
        processed += 1
        print(f"✓ {query_name} -> PAL2NAL 密码子对齐完成: {codon_output}")

    if processed == 0:
        raise SystemExit("未处理任何文件，请检查输入或 PAL2NAL 设置。")

    print(f"输出目录: {output_dir}")
    if axt_dir:
        print(f"AXT 输出目录: {axt_dir}")


if __name__ == "__main__":
    main()
