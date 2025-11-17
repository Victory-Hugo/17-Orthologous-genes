#!/usr/bin/env python3
# 作者: BigLin
# 依赖: python3 (标准库: argparse, pathlib, sys)

import argparse
import sys
from pathlib import Path

# ==============================
# 默认配置，可被命令行/环境变量覆盖
# ==============================
DEFAULT_BASE_DIR = Path("/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output")
DEFAULT_OUTPUT_TSV = Path("/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/presence_absence_matrix.tsv")
DEFAULT_MERGED_TSV = Path("/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM/output/tbl_merge.tsv")
DEFAULT_EVALUE_CUTOFF = 1e-10  # E-value 阈值，小于该值记为1

# 控制台颜色
COLOR_INFO = "\033[94m"
COLOR_SUCCESS = "\033[92m"
COLOR_WARNING = "\033[93m"
COLOR_RESET = "\033[0m"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="读取 Hmmsearch 生成的 *.tbl，输出存在缺失矩阵和合并表。"
    )
    parser.add_argument(
        "--base-dir",
        default=str(DEFAULT_BASE_DIR),
        help="包含 *.tbl 结果的根目录（默认: %(default)s）",
    )
    parser.add_argument(
        "--output-tsv",
        default=str(DEFAULT_OUTPUT_TSV),
        help="存在缺失矩阵输出路径（默认: %(default)s）",
    )
    parser.add_argument(
        "--merged-tsv",
        default=str(DEFAULT_MERGED_TSV),
        help="合并后的 tbl 汇总输出路径（默认: %(default)s）",
    )
    parser.add_argument(
        "--evalue-cutoff",
        type=float,
        default=DEFAULT_EVALUE_CUTOFF,
        help="E-value 小于该值视为存在 (默认: %(default)s)",
    )
    return parser.parse_args()


def parse_tbl_file(file_path: Path, evalue_cutoff: float):
    """
    解析单个 .tbl 文件。

    返回:
        rows: List[List[str]] => 数据行, 每行长度为19
        has_hit: True/False/None => e-value < 阈值时为True; 文件可读但不满足为False; 读取失败为None
    """
    rows = []
    has_hit = False
    try:
        with file_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                parts = stripped.split(maxsplit=18)
                if len(parts) < 19:
                    missing = 19 - len(parts)
                    parts.extend(["" for _ in range(missing)])
                    print(f"{COLOR_WARNING}[WARN]{COLOR_RESET} 行格式异常, 已补齐空字段: {file_path}")
                rows.append(parts)
                e_full = _parse_float_evalue(parts[4])
                if e_full is not None and e_full < evalue_cutoff:
                    has_hit = True
        return rows, (True if has_hit else False)
    except OSError as exc:
        print(f"{COLOR_WARNING}[WARN]{COLOR_RESET} 无法读取文件 {file_path}: {exc}")
        return [], None


def collect_tbl_files(base_dir: Path, evalue_cutoff: float):
    """
    遍历指定目录下的所有 .tbl 文件。
    返回:
        tbl_pairs: List[(gene_name, genome_name, has_hit)]
        merged_rows: List[(genome_name, row_values)]

    其中:
    - gene_name = 基因目录名 (比如 geneA)
    - genome_name = 文件基名 (比如 genome_001)
    - has_hit:
        True  => 该tbl内存在 e-value 小于阈值的命中
        False => 该tbl存在，但没有任何满足阈值的命中
        None  => 文件存在但读取失败 (I/O错误等)
    """
    tbl_paths = [path for path in base_dir.rglob("*.tbl") if path.is_file()]
    total = len(tbl_paths)
    if total == 0:
        print(f"{COLOR_WARNING}[WARN]{COLOR_RESET} 未在 {base_dir} 找到 .tbl 文件")
        return [], []

    tbl_pairs = []
    merged_rows = []
    for index, file_path in enumerate(sorted(tbl_paths), start=1):
        gene_name = file_path.parent.name
        genome_name = file_path.stem
        rows, has_hit = parse_tbl_file(file_path, evalue_cutoff)
        for row in rows:
            merged_rows.append((genome_name, row))
        tbl_pairs.append((gene_name, genome_name, has_hit))
        render_progress(index, total, prefix="[COLLECT]")

    if total:
        sys.stdout.write("\n")
        sys.stdout.flush()
    return tbl_pairs, merged_rows


def build_presence_absence_matrix(tbl_pairs):
    """
    构建存在/缺失矩阵 (presence/absence matrix)。

    规则:
    - 初始化时，所有 [genome][gene] = "NA"
      含义: 这个组合没有对应的 .tbl 文件 => 没跑，不知道
    - 如果有该 .tbl 文件:
        has_hit is True  => 1
        has_hit is False => 0
        has_hit is None  => "NA" (读文件失败，视为未知)

    返回:
    gene_names (按字母序)
    genome_names (按字母序)
    presence[genome][gene] -> "NA" / 0 / 1
    """
    gene_names = sorted({gene for gene, _, _ in tbl_pairs})
    genome_names = sorted({genome for _, genome, _ in tbl_pairs})

    # 全部预填为 "NA"
    presence = {
        genome: {gene: "NA" for gene in gene_names}
        for genome in genome_names
    }

    # 用真实有文件的组合覆盖 "NA"
    for gene, genome, has_hit in tbl_pairs:
        if has_hit is True:
            presence[genome][gene] = 1
        elif has_hit is False:
            presence[genome][gene] = 0
        else:
            # has_hit is None，维持/设为 "NA"
            presence[genome][gene] = "NA"

    return gene_names, genome_names, presence


def _parse_float_evalue(token: str):
    """安全解析 e-value 字符串为 float。无法解析返回 None。"""
    t = token.strip()
    if t in {"*", "-", "NA", "na"}:
        return None
    try:
        return float(t)
    except ValueError:
        return None


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
    """
    写出 TSV 矩阵:
    第一列: 文件名(即genome名)
    后面一列列是各gene
    单元格为 1 / 0 / NA
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        header = ["文件名", *gene_names]
        handle.write("\t".join(header) + "\n")
        for genome in genome_names:
            row = [genome] + [str(presence[genome][gene]) for gene in gene_names]
            handle.write("\t".join(row) + "\n")


def write_merged_tsv(output_path: Path, merged_rows):
    """
    将所有 .tbl 数据汇总写入一个 TSV 文件。
    第一列: 文件名(即 genome 名)
    """
    header = [
        "文件名",
        "target name",
        "accession",
        "query name",
        "accession",
        "E-value",
        "score",
        "bias",
        "E-value",
        "score",
        "bias",
        "exp",
        "reg",
        "clu",
        "ov",
        "env",
        "dom",
        "rep",
        "inc",
        "description of target",
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(header) + "\n")
        for genome_name, row in merged_rows:
            line = [genome_name, *row]
            handle.write("\t".join(line) + "\n")


def main():
    args = parse_args()
    base_dir = Path(args.base_dir).expanduser()
    output_tsv = Path(args.output_tsv).expanduser()
    merged_tsv = Path(args.merged_tsv).expanduser()
    evalue_cutoff = args.evalue_cutoff

    print(f"{COLOR_INFO}[INFO]{COLOR_RESET} 开始扫描 {base_dir}")
    if not base_dir.exists():
        print(f"{COLOR_WARNING}[WARN]{COLOR_RESET} 目录不存在: {base_dir}")
        sys.exit(1)

    print(f"{COLOR_INFO}[INFO]{COLOR_RESET} 使用 e-value 阈值: {evalue_cutoff:g}")

    tbl_pairs, merged_rows = collect_tbl_files(base_dir, evalue_cutoff)
    total_files = len(tbl_pairs)
    total_entries = len(merged_rows)
    print(f"{COLOR_INFO}[INFO]{COLOR_RESET} 共发现 {total_files} 个 .tbl 文件")

    gene_names, genome_names, presence = build_presence_absence_matrix(tbl_pairs)
    print(f"{COLOR_INFO}[INFO]{COLOR_RESET} 基因数量: {len(gene_names)} | 样本数量: {len(genome_names)}")

    write_matrix(output_tsv, gene_names, genome_names, presence)
    print(f"{COLOR_SUCCESS}[DONE]{COLOR_RESET} 存在缺失矩阵已写入: {output_tsv}")

    write_merged_tsv(merged_tsv, merged_rows)
    print(f"{COLOR_INFO}[INFO]{COLOR_RESET} 合并条目数量: {total_entries}")
    print(f"{COLOR_SUCCESS}[DONE]{COLOR_RESET} 合并结果已写入: {merged_tsv}")


if __name__ == "__main__":
    main()
