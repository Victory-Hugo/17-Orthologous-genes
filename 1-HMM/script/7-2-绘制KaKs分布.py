#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
根据 KaKs_summary.tsv 绘制常用分布图
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="绘制 Ka/Ks 分布图")
    parser.add_argument("--summary", required=True, help="7-1 生成的 KaKs_summary.tsv")
    parser.add_argument("--output-dir", required=True, help="图像输出目录")
    parser.add_argument("--method", default=None, help="仅使用指定 Method 的结果")
    return parser.parse_args()


def load_data(summary_path: Path, method: str | None) -> pd.DataFrame:
    df = pd.read_csv(summary_path, sep="\t")
    if method:
        df = df[df["Method"] == method]
    return df


def plot_histogram(df: pd.DataFrame, output_path: Path) -> None:
    plt.figure(figsize=(6, 4))
    plt.hist(df["Ka/Ks"], bins=30, color="#1f77b4", alpha=0.8)
    plt.xlabel("Ka/Ks")
    plt.ylabel("Count")
    plt.title("Ka/Ks 分布")
    plt.axvline(1, color="red", linestyle="--", label="Ka/Ks = 1")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_scatter(df: pd.DataFrame, output_path: Path) -> None:
    plt.figure(figsize=(5, 5))
    plt.scatter(df["Ks"], df["Ka"], s=30, alpha=0.7, c=df["Ka/Ks"], cmap="viridis")
    plt.xlabel("Ks")
    plt.ylabel("Ka")
    plt.title("Ka vs Ks")
    plt.colorbar(label="Ka/Ks")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_category_bar(df: pd.DataFrame, output_path: Path) -> None:
    max_val = max(df["Ka/Ks"].max(), 1.01)
    bins = pd.cut(
        df["Ka/Ks"],
        bins=[0, 0.5, 1, max_val],
        labels=["Ka/Ks < 0.5", "0.5-1", "> 1"],
        include_lowest=True,
    )
    counts = bins.value_counts().reindex(["Ka/Ks < 0.5", "0.5-1", "> 1"]).fillna(0)
    plt.figure(figsize=(6, 4))
    counts.plot(kind="bar", color="#ff7f0e")
    plt.ylabel("Count")
    plt.title("Ka/Ks 分类统计")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    df = load_data(Path(args.summary), args.method)
    if df.empty:
        raise SystemExit("输入数据为空，无法绘图")
    plot_histogram(df, output_dir / "kaks_hist.png")
    plot_scatter(df, output_dir / "ka_vs_ks.png")
    plot_category_bar(df, output_dir / "kaks_category.png")
    print(f"图像输出目录: {output_dir}")


if __name__ == "__main__":
    main()
