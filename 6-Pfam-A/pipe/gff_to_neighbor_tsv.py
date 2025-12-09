#!/usr/bin/env python3
"""
从 GFF 生成邻域注释 TSV，供 rv0194_filter.py 使用。
逻辑：
- 读取目录下所有 .gff（或单个文件），提取 CDS 的 protein_id 和 product
- 按染色体/contig 的顺序，取每个 CDS 上下游固定窗口（默认 +/-5 个 CDS），汇总 product 字段
- 输出列：ProteinID\tNeighborProducts（分号分隔，已转为小写）
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple


def parse_gff_cds(path: Path) -> Dict[str, List[Tuple[int, str]]]:
    """
    返回 seqid -> [(start, protein_id, product), ...]，按 start 升序。
    """
    seq_map: Dict[str, List[Tuple[int, str, str]]] = {}
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _, feature, start_s, _end_s, _score, _strand, _phase, attrs = parts
            if feature != "CDS":
                continue
            start = int(start_s)
            attr_map = {}
            for field in attrs.split(";"):
                if not field:
                    continue
                if "=" in field:
                    key, val = field.split("=", 1)
                    attr_map[key] = val
            pid = attr_map.get("protein_id") or attr_map.get("Name")
            product = attr_map.get("product", "")
            if not pid:
                continue
            seq_map.setdefault(seqid, []).append((start, pid, product))
    for seqid in seq_map:
        seq_map[seqid].sort(key=lambda x: x[0])
    # 转换为 seqid -> [(start, pid, product)]
    return seq_map


def iter_neighbor_products(
    seq_map: Dict[str, List[Tuple[int, str, str]]],
    window: int,
    target_pids: Optional[Set[str]] = None,
) -> Iterable[Tuple[str, str]]:
    """
    生成器：对每个（或指定的）protein_id 汇总邻域 products，包含自身 product，范围为上下游 window 个 CDS。
    target_pids 为空则全量输出；否则仅输出目标集合。
    """
    # 建立 pid -> (seqid, idx) 映射，便于快速找到目标
    pid_index: Dict[str, Tuple[str, int]] = {}
    for seqid, entries in seq_map.items():
        for idx, (_start, pid, _prod) in enumerate(entries):
            pid_index[pid] = (seqid, idx)

    if target_pids is None:
        target_pids = set(pid_index.keys())

    for pid in target_pids:
        if pid not in pid_index:
            continue
        seqid, idx = pid_index[pid]
        entries = seq_map[seqid]
        start_idx = max(0, idx - window)
        end_idx = min(len(entries), idx + window + 1)
        products = [entries[i][2] for i in range(start_idx, end_idx) if entries[i][2]]
        yield pid, ";".join(products).lower()


def write_header(out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        f.write("ProteinID\tNeighborProducts\n")


def iter_gff_files(input_path: Path) -> List[Path]:
    if input_path.is_dir():
        return sorted(input_path.glob("*.gff"))
    if input_path.is_file():
        return [input_path]
    return []


def load_targets(target_file: Path) -> Dict[str, Set[str]]:
    """
    从候选列表加载目标：
    - 每行形如 Assembly|ProteinID，前缀作为组装 ID（文件名 stem），后半为 protein_id
    - 如果不含 '|'，则仅按 protein_id 匹配，不区分组装
    返回 assembly -> set(pid)
    """
    targets: Dict[str, Set[str]] = {}
    if not target_file:
        return targets
    with target_file.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if "|" in line:
                assembly, pid = line.split("|", 1)
                if assembly and pid:
                    targets.setdefault(assembly, set()).add(pid)
            else:
                targets.setdefault("", set()).add(line)
    return targets


def main() -> None:
    parser = argparse.ArgumentParser(description="从 GFF 生成邻域注释 TSV")
    parser.add_argument("--gff", required=True, type=Path, help="GFF 文件或目录（*.gff）")
    parser.add_argument("--output", required=True, type=Path, help="输出 TSV 路径")
    parser.add_argument("--window", type=int, default=5, help="上下游取多少个 CDS（默认 5）")
    parser.add_argument(
        "--targets",
        type=Path,
        help="候选列表文件，仅输出其中的蛋白（格式: Assembly|ProteinID 或单独 ProteinID）",
    )
    args = parser.parse_args()

    targets = load_targets(args.targets) if args.targets else {}
    gff_files = iter_gff_files(args.gff)
    if not gff_files:
        raise SystemExit(f"未找到 GFF：{args.gff}")

    write_header(args.output)
    total = 0
    with args.output.open("a", encoding="utf-8") as out_f:
        for gff in gff_files:
            assembly = gff.stem
            # 如果提供了 targets 且当前 assembly 不在其中，跳过
            if targets and assembly not in targets and "" not in targets:
                continue
            seq_map = parse_gff_cds(gff)
            target_pids = None
            if targets:
                # 优先使用对应 assembly 的集合，其次全局（"" 键）
                target_pids = set()
                if assembly in targets:
                    target_pids.update(targets[assembly])
                if "" in targets:
                    target_pids.update(targets[""])
            for pid, text in iter_neighbor_products(seq_map, args.window, target_pids):
                out_f.write(f"{pid}\t{text}\n")
                total += 1
            print(f"[INFO] 处理 {gff.name}，累积条目 {total}")

    print(f"[INFO] 完成，写入邻域 TSV：{args.output}，蛋白条目数：{total}")


if __name__ == "__main__":
    main()
