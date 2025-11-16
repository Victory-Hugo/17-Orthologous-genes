#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
读取 KaKs_Calculator 结果并输出统计摘要
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import List

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Ka/Ks 结果统计")
    parser.add_argument("--kaks-dir", required=True, help="KaKs 输出所在目录")
    parser.add_argument("--output-dir", required=True, help="统计结果输出目录")
    return parser.parse_args()


def read_kaks_files(kaks_dir: Path) -> pd.DataFrame:
    frames: List[pd.DataFrame] = []
    pattern = re.compile(r"(?P<query>.+)\.ref_(?P<ref>.+)\.axt\.kaks\.tsv$")
    for path in sorted(kaks_dir.glob("*.kaks.tsv")):
        match = pattern.match(path.name)
        query = match.group("query") if match else "NA"
        reference = match.group("ref") if match else "NA"
        df = pd.read_csv(path, sep="\t")
        if "Sequence" not in df.columns:
            continue
        df["Query"] = query
        df["Reference"] = reference
        frames.append(df)
    if not frames:
        raise RuntimeError("未读取到任何 KaKs 结果")
    return pd.concat(frames, ignore_index=True)


def summarize(df: pd.DataFrame, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "KaKs_summary.tsv"
    df.to_csv(summary_path, sep="\t", index=False)

    stats = {
        "总比较数": len(df),
        "Method 数量": df["Method"].nunique(),
        "Ka 中位数": df["Ka"].median(),
        "Ks 中位数": df["Ks"].median(),
        "Ka/Ks 中位数": df["Ka/Ks"].median(),
        "Ka/Ks > 1": (df["Ka/Ks"] > 1).sum(),
        "Ka/Ks 0.5~1": ((df["Ka/Ks"] > 0.5) & (df["Ka/Ks"] <= 1)).sum(),
    }
    stats_path = output_dir / "KaKs_stats.tsv"
    pd.Series(stats).to_csv(stats_path, sep="\t", header=False)

    per_method = (
        df.groupby("Method")[["Ka", "Ks", "Ka/Ks"]]
        .agg(["mean", "median"])
        .reset_index()
    )
    per_method.to_csv(output_dir / "KaKs_per_method.tsv", sep="\t", index=False)

    print("===== Ka/Ks 摘要 =====")
    for key, value in stats.items():
        print(f"{key}: {value}")
    print(f"详细结果: {summary_path}")


def main() -> None:
    args = parse_args()
    kaks_dir = Path(args.kaks_dir)
    df = read_kaks_files(kaks_dir)
    summarize(df, Path(args.output_dir))


if __name__ == "__main__":
    main()
