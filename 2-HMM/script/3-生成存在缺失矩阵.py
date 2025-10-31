#!/usr/bin/env python3
# 作者: BigLin
# 依赖: python3 (标准库: os, pathlib, sys)

import os
import sys
from pathlib import Path

# 硬编码的配置路径
#* 输入文件结构如下
#* ├── geneA/
#* │   ├── genome_001.tbl
#* │   ├── genome_002.tbl
#* │   └── genome_003.tbl
#* │
#* ├── geneB/
#* │   ├── genome_001.tbl
#* │   ├── genome_002.tbl
#* │   └── genome_004.tbl
BASE_DIR = Path("/home/luolintao/0-tmp/3-Bam_Tam/output") #* 包含所有tbl的文件夹，具体如上↑
OUTPUT_TSV = Path("/home/luolintao/0-tmp/3-Bam_Tam/output/presence_absence_matrix.tsv")

# 控制台输出的颜色样式 (ANSI 颜色代码)
COLOR_INFO = "\033[94m"
COLOR_SUCCESS = "\033[92m"
COLOR_WARNING = "\033[93m"
COLOR_RESET = "\033[0m"


def collect_tbl_files(base_dir: Path):
	"""
	遍历指定目录，查找所有以 `.tbl` 结尾的文件，
	并为每个文件提取基因名（文件所在文件夹名）、基因组名（文件名），
	以及该文件是否包含有效比对结果（hit）。
	
	返回:
		list[tuple[str, str, bool]]: 每个元素为 (gene_name, genome_name, has_hit)
	"""
	tbl_paths = [path for path in base_dir.rglob("*.tbl") if path.is_file()]
	total = len(tbl_paths)
	if total == 0:
		print(f"{COLOR_WARNING}[WARN]{COLOR_RESET} 未在 {base_dir} 找到 .tbl 文件")
		return []
	tbl_pairs = []
	for index, file_path in enumerate(sorted(tbl_paths), start=1):
		gene_name = file_path.parent.name  # 基因名称：文件所在的上一级目录名
		genome_name = file_path.stem       # 基因组名称：文件名（去掉扩展名）
		has_hit = tbl_has_hit(file_path)   # 检查文件是否有命中记录
		tbl_pairs.append((gene_name, genome_name, has_hit))
		render_progress(index, total, prefix="[COLLECT]")  # 渲染进度条
	if total:
		sys.stdout.write("\n")
		sys.stdout.flush()
	return tbl_pairs


def build_presence_absence_matrix(tbl_pairs):
	"""
	根据 (gene, genome, hit) 三元组，构建存在/缺失矩阵。
	矩阵以字典形式存储，键为基因组，值为对应的基因存在状态（0/1）。

	返回:
		tuple[list[str], list[str], dict[str, dict[str, int]]]:
			基因名列表、基因组名列表、presence矩阵字典
	"""
	gene_names = sorted({gene for gene, _, _ in tbl_pairs})
	genome_names = sorted({genome for _, genome, _ in tbl_pairs})

	# 初始化 presence 矩阵，所有值初始为 0
	presence = {genome: {gene: 0 for gene in gene_names} for genome in genome_names}
	for gene, genome, has_hit in tbl_pairs:
		presence[genome][gene] = 1 if has_hit else 0
	return gene_names, genome_names, presence


def tbl_has_hit(file_path: Path) -> bool:
	"""
	检测 .tbl 文件是否包含至少一条非注释的命中记录（hit）。
	只要出现一行非空且非注释的内容，就认为存在命中。

	返回:
		bool: 若文件中存在命中行则为 True，否则为 False
	"""
	try:
		with file_path.open("r", encoding="utf-8", errors="ignore") as handle:
			for line in handle:
				stripped = line.strip()
				if not stripped or stripped.startswith("#"):
					continue
				return True
	except OSError as exc:
		print(f"{COLOR_WARNING}[WARN]{COLOR_RESET} 无法读取文件 {file_path}: {exc}")
	return False


def render_progress(current, total, prefix=""):
	"""
	在命令行中渲染一个简单的进度条。

	参数:
		current (int): 当前进度位置
		total (int): 总任务数量
		prefix (str): 进度条前缀字符串
	"""
	if total == 0:
		sys.stdout.write(f"{COLOR_WARNING}{prefix} 没有可处理的项目{COLOR_RESET}\n")
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
	"""
	将 presence/absence 矩阵写入 TSV 文件。

	参数:
		output_path (Path): 输出文件路径
		gene_names (list[str]): 基因名列表
		genome_names (list[str]): 基因组名列表
		presence (dict): 存在矩阵字典
	"""
	output_path.parent.mkdir(parents=True, exist_ok=True)
	with output_path.open("w", encoding="utf-8") as handle:
		header = ["文件名", *gene_names]
		handle.write("\t".join(header) + "\n")
		for genome in genome_names:
			row = [genome] + [str(presence[genome][gene]) for gene in gene_names]
			handle.write("\t".join(row) + "\n")


def main():
	"""
	程序入口。
	扫描指定目录下的所有 .tbl 文件，解析并生成基因存在/缺失矩阵，
	最终将结果输出为 TSV 文件。
	"""
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
