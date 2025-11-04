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

# ==============================
# 硬编码配置
# ==============================
BASE_DIR = Path("/home/luolintao/0-tmp/3-Bam_Tam/output")  # 包含所有tbl的文件夹
OUTPUT_TSV = Path("/home/luolintao/0-tmp/3-Bam_Tam/output/presence_absence_matrix.tsv")
EVALUE_CUTOFF = 1e-5  # 新增：E-value 阈值，小于该值记为1，否则0

# 控制台颜色
COLOR_INFO = "\033[94m"
COLOR_SUCCESS = "\033[92m"
COLOR_WARNING = "\033[93m"
COLOR_RESET = "\033[0m"


def collect_tbl_files(base_dir: Path):
    """
    遍历指定目录下的所有 .tbl 文件。
    返回 (gene_name, genome_name, has_hit)
    """
    tbl_paths = [path for path in base_dir.rglob("*.tbl") if path.is_file()]
    total = len(tbl_paths)
    if total == 0:
        print(f"{COLOR_WARNING}[WARN]{COLOR_RESET} 未在 {base_dir} 找到 .tbl 文件")
        return []
    tbl_pairs = []
    for index, file_path in enumerate(sorted(tbl_paths), start=1):
        gene_name = file_path.parent.name
        genome_name = file_path.stem
        has_hit = tbl_has_hit(file_path)
        tbl_pairs.append((gene_name, genome_name, has_hit))
        render_progress(index, total, prefix="[COLLECT]")
    if total:
        sys.stdout.write("\n")
        sys.stdout.flush()
    return tbl_pairs


def build_presence_absence_matrix(tbl_pairs):
    """
    构建存在/缺失矩阵。
    """
    gene_names = sorted({gene for gene, _, _ in tbl_pairs})
    genome_names = sorted({genome for _, genome, _ in tbl_pairs})
    presence = {genome: {gene: 0 for gene in gene_names} for genome in genome_names}
    for gene, genome, has_hit in tbl_pairs:
        presence[genome][gene] = 1 if has_hit else 0
    return gene_names, genome_names, presence


def _parse_float_evalue(token: str):
    """安全解析 e-value 字符串为 float"""
    t = token.strip()
    if t in {"*", "-", "NA", "na"}:
        return None
    try:
        return float(t)
    except ValueError:
        return None


def tbl_has_hit(file_path: Path) -> bool:
    """
    检查 .tbl 文件是否有符合 EVALUE_CUTOFF 的命中：
    - 使用 “full sequence E-value” （第5列）
    - 若小于 EVALUE_CUTOFF，则返回 True，否则 False。
    """
    try:
        with file_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                cols = stripped.split()
                if len(cols) < 5:
                    continue
                e_full = _parse_float_evalue(cols[4])
                if e_full is not None and e_full < EVALUE_CUTOFF:
                    return True
    except OSError as exc:
        print(f"{COLOR_WARNING}[WARN]{COLOR_RESET} 无法读取文件 {file_path}: {exc}")
    return False


def render_progress(current, total, prefix=""):
    """简单进度条"""
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
    """写出 TSV 矩阵"""
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
    print(f"{COLOR_INFO}[INFO]{COLOR_RESET} 使用 e-value 阈值: {EVALUE_CUTOFF:g}")

    tbl_pairs = collect_tbl_files(BASE_DIR)
    total_files = len(tbl_pairs)
    print(f"{COLOR_INFO}[INFO]{COLOR_RESET} 共发现 {total_files} 个 .tbl 文件")

    gene_names, genome_names, presence = build_presence_absence_matrix(tbl_pairs)
    print(f"{COLOR_INFO}[INFO]{COLOR_RESET} 基因数量: {len(gene_names)} | 样本数量: {len(genome_names)}")

    write_matrix(OUTPUT_TSV, gene_names, genome_names, presence)
    print(f"{COLOR_SUCCESS}[DONE]{COLOR_RESET} 结果已写入: {OUTPUT_TSV}")


if __name__ == "__main__":
    main()
