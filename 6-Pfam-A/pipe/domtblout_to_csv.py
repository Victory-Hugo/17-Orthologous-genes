#!/usr/bin/env python3
"""
将 hmmscan --domtblout 的结果转换为 CSV，并按需删除原始 domtblout 文件。
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable, List


DOMTBL_FIELDS: List[str] = [
    "target_name",
    "target_accession",
    "tlen",
    "query_name",
    "query_accession",
    "qlen",
    "full_sequence_E_value",
    "full_sequence_score",
    "full_sequence_bias",
    "domain_number",
    "domain_total",
    "c_E_value",
    "i_E_value",
    "domain_score",
    "domain_bias",
    "hmm_from",
    "hmm_to",
    "ali_from",
    "ali_to",
    "env_from",
    "env_to",
    "acc",
    "description_of_target",
]


def parse_domtblout_lines(lines: Iterable[str]) -> Iterable[List[str]]:
    """迭代 domtblout 行，返回按列拆分的字段列表。"""
    for line in lines:
        if not line or line.startswith("#"):
            continue
        # 22 次分割得到 23 列，末尾列保留含空格的描述。
        parts = line.strip().split(maxsplit=22)
        if len(parts) != len(DOMTBL_FIELDS):
            raise ValueError(f"解析列数不符（期望 {len(DOMTBL_FIELDS)} 列，得到 {len(parts)} 列）：{line}")
        yield parts


def convert_domtblout_to_csv(domtbl_path: Path, output_dir: Path) -> Path:
    """将单个 domtblout 文件转换为 CSV，返回生成的 CSV 路径。"""
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / (domtbl_path.stem + ".csv")
    with domtbl_path.open("r", encoding="utf-8") as infile, csv_path.open("w", newline="", encoding="utf-8") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(DOMTBL_FIELDS)
        for row in parse_domtblout_lines(infile):
            writer.writerow(row)
    return csv_path


def run(input_dir: Path, output_dir: Path, delete_original: bool = True) -> None:
    """
    批量转换目录内的 domtblout 文件为 CSV。

    :param input_dir: 存放 domtblout 文件的目录
    :param output_dir: 输出 CSV 的目录
    :param delete_original: 是否删除原始 domtblout 文件
    """
    domtbl_files = sorted(input_dir.glob("*.domtblout"))
    if not domtbl_files:
        raise FileNotFoundError(f"未在目录 {input_dir} 中找到 domtblout 文件")

    for domtbl_path in domtbl_files:
        csv_path = convert_domtblout_to_csv(domtbl_path, output_dir)
        if delete_original:
            domtbl_path.unlink()
        print(f"{domtbl_path.name} -> {csv_path.name}")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="将 hmmscan domtblout 结果转换为 CSV，可选择删除原始文件。"
    )
    parser.add_argument(
        "-i",
        "--input-dir",
        type=Path,
        required=True,
        help="存放 domtblout 的目录",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        required=True,
        help="CSV 输出目录",
    )
    parser.add_argument(
        "--keep-original",
        action="store_true",
        help="保留原始 domtblout 文件（默认删除）",
    )
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()
    run(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        delete_original=not args.keep_original,
    )


if __name__ == "__main__":
    main()
