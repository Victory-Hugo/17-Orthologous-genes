#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
合并多个 *.kaks.tsv，默认保留中间文件，可选删除。
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable, List


def log(level: str, msg: str) -> None:
    print(f"[{level}] {msg}")


def merge(files: Iterable[Path], merged_path: Path, keep_intermediate: bool) -> None:
    files = list(files)
    if not files:
        raise SystemExit("未找到 KaKs 结果可合并")

    merged_path.parent.mkdir(parents=True, exist_ok=True)
    with merged_path.open("w", encoding="utf-8") as fout:
        for idx, f in enumerate(files):
            with f.open("r", encoding="utf-8") as fin:
                for line_no, line in enumerate(fin):
                    if idx > 0 and line_no == 0:
                        continue
                    fout.write(line)

    log("INFO", f"合并完成 -> {merged_path}")
    if not keep_intermediate:
        for f in files:
            if f != merged_path:
                f.unlink(missing_ok=True)
        log("INFO", "已删除中间 KaKs 结果文件")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="合并 KaKs_Calculator 输出")
    parser.add_argument("--kaks-dir", required=True, help="包含 *.kaks.tsv 的目录")
    parser.add_argument(
        "--output",
        default=None,
        help="合并输出路径（默认: kaks_dir/kaks_merged.tsv）",
    )
    parser.add_argument(
        "--keep-intermediate",
        action="store_true",
        help="保留单个 kaks 结果文件",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    kaks_dir = Path(args.kaks_dir)
    if not kaks_dir.is_dir():
        raise SystemExit(f"缺少 KaKs 目录: {kaks_dir}")

    merged_path = Path(args.output) if args.output else kaks_dir / "kaks_merged.tsv"
    files: List[Path] = [
        f for f in sorted(kaks_dir.glob("*.kaks.tsv")) if f != merged_path
    ]
    merge(files, merged_path, keep_intermediate=args.keep_intermediate)


if __name__ == "__main__":
    try:
        main()
    except SystemExit as exc:
        raise
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise
