#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
并行执行 KaKs_Calculator：读取 AXT，输出对应的 *.kaks.tsv，并记录日志。
"""

from __future__ import annotations

import argparse
import concurrent.futures
import subprocess
import sys
from pathlib import Path
from typing import List, Tuple


def log(level: str, msg: str) -> None:
    print(f"[{level}] {msg}")


def run_one(
    kaks_bin: Path,
    axt: Path,
    kaks_dir: Path,
    codon_table: int | None,
    log_dir: Path,
) -> Tuple[Path, str]:
    output = kaks_dir / (axt.stem + ".kaks.tsv")
    log_file = log_dir / (axt.stem + ".log")

    cmd = [str(kaks_bin), "-i", str(axt), "-o", str(output)]
    if codon_table:
        cmd += ["-c", str(codon_table)]

    result = subprocess.run(cmd, text=True, capture_output=True)
    log_lines: List[str] = []
    log_lines.append("CMD: " + " ".join(cmd))
    if result.stdout:
        log_lines.append("STDOUT:\n" + result.stdout)
    if result.stderr:
        log_lines.append("STDERR:\n" + result.stderr)
    log_file.parent.mkdir(parents=True, exist_ok=True)
    log_file.write_text("\n".join(log_lines), encoding="utf-8")

    if result.returncode != 0:
        raise RuntimeError(f"KaKs 失败: {axt.name}，详见 {log_file}")

    if not output.exists():
        raise RuntimeError(f"KaKs 未生成输出: {output}")

    return output, str(log_file)


def split_axt_file(axt_path: Path, split_blocks: int, work_root: Path) -> List[Path]:
    """
    将单个 AXT 按 block 数拆分为若干小 AXT，返回拆分后的文件列表。
    若 split_blocks <= 0 或 block 数不超过 split_blocks，则返回原文件。
    """
    if split_blocks <= 0:
        return [axt_path]

    blocks: List[Tuple[str, str, str]] = []
    with axt_path.open("r", encoding="utf-8") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            header = header.strip()
            if not header:
                continue
            seq1 = handle.readline().strip()
            seq2 = handle.readline().strip()
            # 跳过可能的空行
            handle.readline()
            blocks.append((header, seq1, seq2))

    if len(blocks) <= split_blocks:
        return [axt_path]

    out_dir = work_root / axt_path.stem
    out_dir.mkdir(parents=True, exist_ok=True)
    parts: List[Path] = []

    def write_part(part_idx: int, part_blocks: List[Tuple[str, str, str]]) -> Path:
        out_path = out_dir / f"{axt_path.stem}.part{part_idx}.axt"
        with out_path.open("w", encoding="utf-8") as out:
            for idx, (header, seq1, seq2) in enumerate(part_blocks, start=1):
                tokens = header.split()
                suffix = " ".join(tokens[1:]) if len(tokens) > 1 else ""
                out.write(f"{idx} {suffix}\n")
                out.write(seq1 + "\n")
                out.write(seq2 + "\n\n")
        return out_path

    current: List[Tuple[str, str, str]] = []
    part_idx = 0
    for block in blocks:
        current.append(block)
        if len(current) >= split_blocks:
            part_idx += 1
            parts.append(write_part(part_idx, current))
            current = []
    if current:
        part_idx += 1
        parts.append(write_part(part_idx, current))

    return parts


def merge_kaks(files: List[Path], merged_path: Path) -> None:
    merged_path.parent.mkdir(parents=True, exist_ok=True)
    with merged_path.open("w", encoding="utf-8") as fout:
        for idx, f in enumerate(files):
            with f.open("r", encoding="utf-8") as fin:
                for line_no, line in enumerate(fin):
                    if idx > 0 and line_no == 0:
                        continue
                    fout.write(line)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="批量运行 KaKs_Calculator（支持单个 AXT 拆分并发）")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--axt", help="单个 AXT 文件（绝对路径推荐）")
    group.add_argument("--axt-dir", help="AXT 输入目录（处理全部 *.axt）")
    parser.add_argument("--kaks-dir", required=True, help="KaKs 输出目录")
    parser.add_argument("--kaks-bin", required=True, help="KaKs_Calculator 可执行文件")
    parser.add_argument("--kaks-codon-table", type=int, default=None, help="KaKs_Calculator -c 参数")
    parser.add_argument("--kaks-cpu", type=int, default=1, help="并行任务数")
    parser.add_argument("--log-dir", default=None, help="日志目录（默认: kaks_dir/logs）")
    parser.add_argument("--split-blocks", type=int, default=0, help="将单个 AXT 按 block 数拆分并并发处理，0 表示不拆分")
    parser.add_argument("--keep-chunks", action="store_true", help="保留拆分后的 AXT 与对应的 KaKs 结果")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    kaks_dir = Path(args.kaks_dir)
    kaks_bin = Path(args.kaks_bin)
    log_dir = Path(args.log_dir) if args.log_dir else kaks_dir / "logs"

    if not kaks_bin.exists():
        raise SystemExit(f"未找到 KaKs 可执行文件: {kaks_bin}")

    kaks_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    if args.axt:
        axt_files = [Path(args.axt)]
    else:
        axt_dir = Path(args.axt_dir)
        if not axt_dir.is_dir():
            raise SystemExit(f"缺少 AXT 目录: {axt_dir}")
        axt_files = sorted(axt_dir.glob("*.axt"))

    if not axt_files:
        raise SystemExit("未找到需要处理的 AXT 文件")

    max_workers = max(1, args.kaks_cpu)

    for axt in axt_files:
        if not axt.exists():
            raise SystemExit(f"AXT 文件不存在: {axt}")
        log("INFO", f"处理 AXT: {axt} (拆分阈值 {args.split_blocks} blocks, 并发 {max_workers})")

        parts = split_axt_file(axt, args.split_blocks, work_root=kaks_dir / "axt_split")
        part_outputs: List[Path] = []
        errors: List[str] = []

        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
            futs = {
                ex.submit(
                    run_one,
                    kaks_bin=kaks_bin,
                    axt=part,
                    kaks_dir=kaks_dir,
                    codon_table=args.kaks_codon_table,
                    log_dir=log_dir,
                ): part
                for part in parts
            }
            for fut in concurrent.futures.as_completed(futs):
                part = futs[fut]
                try:
                    out_file, log_file = fut.result()
                    part_outputs.append(out_file)
                    log("INFO", f"完成: {part.name} -> {out_file.name} (日志: {log_file})")
                except Exception as exc:  # noqa: BLE001
                    errors.append(f"{part.name}: {exc}")

        if errors:
            log("ERROR", f"{axt.name} 部分拆分任务失败：")
            for err in errors:
                log("ERROR", f"  - {err}")
            raise SystemExit(1)

        if len(parts) > 1:
            merged_path = kaks_dir / (axt.stem + ".kaks.tsv")
            merge_kaks(part_outputs, merged_path)
            if not args.keep_chunks:
                for f in part_outputs:
                    f.unlink(missing_ok=True)
                for part in parts:
                    part.unlink(missing_ok=True)
            log("INFO", f"{axt.name} 拆分并合并完成 -> {merged_path}")
        else:
            log("INFO", f"{axt.name} 无需拆分，结果: {part_outputs[0]}")

    log("INFO", f"全部完成，共处理 {len(axt_files)} 个 AXT。")


if __name__ == "__main__":
    try:
        main()
    except SystemExit as exc:
        raise
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise
