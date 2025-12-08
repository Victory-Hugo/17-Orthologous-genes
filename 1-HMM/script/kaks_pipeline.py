#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
轻量 KaKs 流水线：从已有的密码子比对开始，构建 AXT -> KaKs -> 合并 -> 统计 -> 绘图。
Shell 仅负责传参，核心逻辑集中在此。
"""

from __future__ import annotations

import argparse
import concurrent.futures
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List


def log(level: str, msg: str) -> None:
    print(f"[{level}] {msg}")


def run_cmd(cmd: List[str]) -> None:
    log("INFO", "执行命令: " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def ensure_symlinks(codon_dir: Path) -> None:
    """pal2nal 输出为 *.codon.aln.fna，旧脚本读取 *.codon.fna，这里建立兼容软链。"""
    codon_files = sorted(codon_dir.glob("*.codon.aln.fna"))
    if not codon_files:
        raise SystemExit(f"在 {codon_dir} 未找到 *.codon.aln.fna")
    for f in codon_files:
        target = f.with_name(f.stem.replace(".codon.aln", ".codon") + ".fna")
        if target.exists():
            continue
        os.symlink(f, target)
        log("INFO", f"软链: {target.name} -> {f.name}")


def build_axt(py_axt: Path, codon_dir: Path, axt_dir: Path, reference: str, queries: List[str]) -> None:
    cmd = [
        sys.executable,
        str(py_axt),
        "--alignment-dir",
        str(codon_dir),
        "--output-dir",
        str(axt_dir),
    ]
    if reference:
        cmd += ["--reference", reference]
    if queries:
        cmd += ["--queries", *queries]
    run_cmd(cmd)


def run_kaks_one(kaks_bin: Path, axt: Path, output: Path, codon_table: int | None) -> None:
    cmd = [str(kaks_bin), "-i", str(axt), "-o", str(output)]
    if codon_table:
        cmd += ["-c", str(codon_table)]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def run_kaks_parallel(
    kaks_bin: Path, axt_dir: Path, kaks_dir: Path, codon_table: int | None, max_workers: int
) -> List[Path]:
    kaks_dir.mkdir(parents=True, exist_ok=True)
    axt_files = sorted(axt_dir.glob("*.axt"))
    if not axt_files:
        raise SystemExit(f"在 {axt_dir} 未找到 AXT 文件")
    outputs: List[Path] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
        futs = []
        for axt in axt_files:
            out = kaks_dir / (axt.stem + ".kaks.tsv")
            outputs.append(out)
            futs.append(ex.submit(run_kaks_one, kaks_bin, axt, out, codon_table))
        for fut in concurrent.futures.as_completed(futs):
            fut.result()
    log("INFO", f"KaKs 完成：生成 {len(outputs)} 个结果文件")
    return outputs


def merge_kaks(files: Iterable[Path], merged_path: Path) -> None:
    merged_path.parent.mkdir(parents=True, exist_ok=True)
    files = list(files)
    if not files:
        raise SystemExit("未找到 KaKs 结果可合并")
    with merged_path.open("w", encoding="utf-8") as fout:
        for idx, f in enumerate(files):
            with f.open("r", encoding="utf-8") as fin:
                for line_no, line in enumerate(fin):
                    if idx > 0 and line_no == 0:
                        continue  # 跳过后续文件表头
                    fout.write(line)
    log("INFO", f"合并完成 -> {merged_path}")
    for f in files:
        if f != merged_path:
            f.unlink(missing_ok=True)
    log("INFO", "已删除中间 KaKs 文件")


def run_stats_and_plot(py_stats: Path, py_plot: Path, kaks_dir: Path, stats_dir: Path, plot_dir: Path, method: str | None) -> None:
    run_cmd(
        [
            sys.executable,
            str(py_stats),
            "--kaks-dir",
            str(kaks_dir),
            "--output-dir",
            str(stats_dir),
        ]
    )
    summary = stats_dir / "KaKs_summary.tsv"
    if not summary.exists():
        raise SystemExit("统计结果缺失: KaKs_summary.tsv")
    plot_cmd = [
        sys.executable,
        str(py_plot),
        "--summary",
        str(summary),
        "--output-dir",
        str(plot_dir),
    ]
    if method:
        plot_cmd += ["--method", method]
    run_cmd(plot_cmd)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="KaKs 轻量流水线（从密码子比对开始）")
    p.add_argument("--codon-dir", required=True)
    p.add_argument("--axt-dir", required=True)
    p.add_argument("--kaks-dir", required=True)
    p.add_argument("--kaks-stats-dir", required=True)
    p.add_argument("--kaks-plot-dir", required=True)
    p.add_argument("--kaks-bin", required=True)
    p.add_argument("--kaks-codon-table", type=int, default=None)
    p.add_argument("--kaks-cpu", type=int, default=1)
    p.add_argument("--reference", default="")
    p.add_argument("--queries", nargs="*", default=[])
    p.add_argument("--py-axt", required=True)
    p.add_argument("--py-stats", required=True)
    p.add_argument("--py-plot", required=True)
    p.add_argument("--kaks-method", default=None)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    codon_dir = Path(args.codon_dir)
    axt_dir = Path(args.axt_dir)
    kaks_dir = Path(args.kaks_dir)
    stats_dir = Path(args.kaks_stats_dir)
    plot_dir = Path(args.kaks_plot_dir)
    kaks_bin = Path(args.kaks_bin)

    if not codon_dir.is_dir():
        raise SystemExit(f"缺少密码子对齐目录: {codon_dir}")
    if not kaks_bin.exists():
        raise SystemExit(f"未找到 KaKs 可执行文件: {kaks_bin}")

    ensure_symlinks(codon_dir)
    axt_dir.mkdir(parents=True, exist_ok=True)
    build_axt(Path(args.py_axt), codon_dir, axt_dir, args.reference, args.queries)

    outputs = run_kaks_parallel(
        kaks_bin=kaks_bin,
        axt_dir=axt_dir,
        kaks_dir=kaks_dir,
        codon_table=args.kaks_codon_table,
        max_workers=max(1, args.kaks_cpu),
    )

    merged_path = kaks_dir / "kaks_merged.tsv"
    merge_kaks(outputs, merged_path)

    run_stats_and_plot(
        py_stats=Path(args.py_stats),
        py_plot=Path(args.py_plot),
        kaks_dir=kaks_dir,
        stats_dir=stats_dir,
        plot_dir=plot_dir,
        method=args.kaks_method,
    )
    log("INFO", "✓ KaKs 流水线完成")


if __name__ == "__main__":
    main()
