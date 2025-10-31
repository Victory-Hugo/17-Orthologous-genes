#!/usr/bin/env python3
# Author: BigLin
# Dependencies: python3 (standard library: os, pathlib, sys)

import os
import sys
from pathlib import Path

# Hard-coded configuration paths
BASE_DIR = Path("/home/luolintao/0-tmp/3-Bam_Tam/output")
OUTPUT_TSV = Path("/home/luolintao/0-tmp/3-Bam_Tam/output/presence_absence_matrix.tsv")

# Console styling (ANSI color codes)
COLOR_INFO = "\033[94m"
COLOR_SUCCESS = "\033[92m"
COLOR_WARNING = "\033[93m"
COLOR_RESET = "\033[0m"


def collect_tbl_files(base_dir: Path):
	"""Walk the base directory and collect (gene_name, genome_name) pairs for each .tbl file found."""
	tbl_paths = [path for path in base_dir.rglob("*.tbl") if path.is_file()]
	total = len(tbl_paths)
	if total == 0:
		print(f"{COLOR_WARNING}[WARN]{COLOR_RESET} 未在 {base_dir} 找到 .tbl 文件")
		return []
	tbl_pairs = []
	for index, file_path in enumerate(sorted(tbl_paths), start=1):
		gene_name = file_path.parent.name
		genome_name = file_path.stem
		tbl_pairs.append((gene_name, genome_name))
		render_progress(index, total, prefix="[COLLECT]")
	if total:
		sys.stdout.write("\n")
		sys.stdout.flush()
	return tbl_pairs


def build_presence_absence_matrix(tbl_pairs):
	"""Build mapping structures required for the presence/absence matrix."""
	gene_names = sorted({gene for gene, _ in tbl_pairs})
	genome_names = sorted({genome for _, genome in tbl_pairs})

	presence = {genome: {gene: 0 for gene in gene_names} for genome in genome_names}
	for gene, genome in tbl_pairs:
		presence[genome][gene] = 1
	return gene_names, genome_names, presence


def render_progress(current, total, prefix=""):
	"""Render a simple ASCII progress bar with colored prefix."""
	if total == 0:
		sys.stdout.write(f"{COLOR_WARNING}{prefix} no items to process{COLOR_RESET}\n")
		sys.stdout.flush()
		return

	bar_length = 40
	ratio = current / total
	filled = int(bar_length * ratio)
	bar = "#" * filled + "." * (bar_length - filled)
	percent = ratio * 100
	sys.stdout.write(
		f"\r{COLOR_INFO}{prefix}{COLOR_RESET} |{bar}| {percent:6.2f}% ({current}/{total})"
	)
	sys.stdout.flush()


def write_matrix(output_path: Path, gene_names, genome_names, presence):
	"""Write the presence/absence matrix to the TSV file."""
	output_path.parent.mkdir(parents=True, exist_ok=True)
	with output_path.open("w", encoding="utf-8") as handle:
		header = ["文件名", *gene_names]
		handle.write("\t".join(header) + "\n")
		for genome in genome_names:
			row = [genome] + [str(presence[genome][gene]) for gene in gene_names]
			handle.write("\t".join(row) + "\n")


def main():
	print(f"{COLOR_INFO}[INFO]{COLOR_RESET} 开始扫描 {BASE_DIR}")
	if not BASE_DIR.exists():
		print(f"{COLOR_WARNING}[WARN]{COLOR_RESET} 目录不存在: {BASE_DIR}")
		sys.exit(1)

	tbl_pairs = collect_tbl_files(BASE_DIR)
	total_files = len(tbl_pairs)
	print(f"{COLOR_INFO}[INFO]{COLOR_RESET} 共发现 {total_files} 个 .tbl 文件")

	gene_names, genome_names, presence = build_presence_absence_matrix(tbl_pairs)
	print(f"{COLOR_INFO}[INFO]{COLOR_RESET} 基因数量: {len(gene_names)} | 样本数量: {len(genome_names)}")

	write_matrix(OUTPUT_TSV, gene_names, genome_names, presence)
	print(f"{COLOR_SUCCESS}[DONE]{COLOR_RESET} 结果已写入: {OUTPUT_TSV}")


if __name__ == "__main__":
	main()
