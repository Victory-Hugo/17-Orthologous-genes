#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
准备 KaKs 输入：在密码子对齐目录内为 *.codon.aln.fna 创建兼容的 *.codon.fna 软链。
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def log(level: str, msg: str) -> None:
    print(f"[{level}] {msg}")


def ensure_symlink(src: Path, dst: Path, force: bool) -> None:
    if dst.exists():
        if dst.is_symlink():
            target = dst.resolve()
            if target == src.resolve():
                log("INFO", f"已存在软链: {dst.name} -> {src.name}")
                return
            if force:
                dst.unlink()
                log("WARN", f"移除指向错误的软链: {dst} (原指向 {target})")
            else:
                log("WARN", f"软链已存在但指向 {target}，跳过 (use --force 可修正)")
                return
        else:
            log("WARN", f"目标已存在且不是软链: {dst}，跳过")
            return
    dst.symlink_to(src)
    log("INFO", f"创建软链: {dst.name} -> {src.name}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="为 pal2nal 输出创建 *.codon.fna 软链")
    parser.add_argument("--codon-dir", required=True, help="包含 *.codon.aln.fna 的目录")
    parser.add_argument(
        "--force",
        action="store_true",
        help="如软链存在但指向错误则覆盖",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    codon_dir = Path(args.codon_dir)
    if not codon_dir.is_dir():
        raise SystemExit(f"缺少密码子对齐目录: {codon_dir}")

    codon_files = sorted(codon_dir.glob("*.codon.aln.fna"))
    if not codon_files:
        raise SystemExit(f"在 {codon_dir} 未找到 *.codon.aln.fna")

    for src in codon_files:
        dst = src.with_name(src.name.replace(".codon.aln.fna", ".codon.fna"))
        ensure_symlink(src, dst, force=args.force)

    log("INFO", f"处理完成，共 {len(codon_files)} 个密码子比对。")


if __name__ == "__main__":
    try:
        main()
    except SystemExit as exc:
        raise
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise
