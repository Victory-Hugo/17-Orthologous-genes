#!/usr/bin/env python3
"""
基于 hmmscan (Pfam) CSV 与 TMHMM 结果筛选 Rv0194 型四域融合结构 (TMD–NBD–TMD–NBD)，可选开启 WalkerA/WalkerB/LSGGQ 保守基序与基因邻域关键词过滤。
"""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


ABC_PFAM_ID = "PF00005"


@dataclass
class Domain:
    pfam_id: str
    start: int
    end: int


@dataclass
class TMRegion:
    start: int
    end: int


def _coerce_int(value: str) -> Optional[int]:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _select_pfam_field(fieldnames: Iterable[str]) -> str:
    # 优先使用明确的 pfam_id，其次 query_accession（hmmsearch 模式下 Pfam ID 在 query 列）
    preferred = ["pfam_id", "query_accession", "query_name"]
    for name in preferred:
        if name in fieldnames:
            return name
    raise ValueError("hmmscan CSV 缺少 pfam_id/query_accession/query_name 列")


def parse_hmmscan_csv(csv_path: Path) -> Dict[str, List[Domain]]:
    """
    解析 hmmscan CSV，返回 protein_id -> ABC_tran 域列表。
    仅保留 pfam_id 以 PF00005 开头的行。
    
    注意：hmmsearch 模式下（HMM profile 搜索序列数据库）：
    - query_* 列是 HMM profile 信息（Pfam ID）
    - target_* 列是蛋白质序列信息（蛋白质 ID）
    """
    protein_to_domains: Dict[str, List[Domain]] = {}
    with csv_path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"{csv_path} 没有表头")
        pfam_field = _select_pfam_field(reader.fieldnames)
        # hmmsearch 模式下，蛋白质 ID 在 target_name 列
        name_field = "target_name" if "target_name" in reader.fieldnames else "query_name"
        for row in reader:
            pfam_raw = row.get(pfam_field, "")
            if not pfam_raw or not pfam_raw.startswith(ABC_PFAM_ID):
                continue
            start = _coerce_int(row.get("ali_from"))
            end = _coerce_int(row.get("ali_to"))
            protein_id = row.get(name_field, "")
            if not protein_id or start is None or end is None:
                continue
            protein_to_domains.setdefault(protein_id, []).append(
                Domain(pfam_id=pfam_raw.split(".")[0], start=start, end=end)
            )
    # 域按起点排序
    for domains in protein_to_domains.values():
        domains.sort(key=lambda d: d.start)
    return protein_to_domains


def load_fasta_seqs(fasta_path: Path) -> Dict[str, str]:
    """
    读取 FASTA，返回 id -> 序列。id 取第一个空格前的字段。
    """
    seqs: Dict[str, str] = {}
    current_id: Optional[str] = None
    buf: List[str] = []
    with fasta_path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id:
                    seqs[current_id] = "".join(buf)
                current_id = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if current_id:
            seqs[current_id] = "".join(buf)
    return seqs


def parse_tmhmm(tmhmm_path: Path) -> Dict[str, List[TMRegion]]:
    """
    解析 TMHMM 输出：
    - 短格式：ProteinID  TMhelix  start  end
    - 长格式：ProteinID  len=xxx  ExpAA=xxx  First60=xxx  PredHel=xxx  Topology=i40-61o76-95i...
    从 Topology 字段中提取 TM 区域（'i' 表示胞内侧，'o' 表示胞外侧）
    """
    tm_map: Dict[str, List[TMRegion]] = {}
    with tmhmm_path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("Temporary files"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            
            protein_id = parts[0]
            
            # 检查是否是长格式（含有 Topology= 字段）
            topology_field = None
            for part in parts:
                if part.startswith("Topology="):
                    topology_field = part.split("=", 1)[1]
                    break
            
            if topology_field:
                # 解析 Topology 字段：i40-61o76-95i144-166o...
                # 提取所有形如 数字-数字 的区域作为 TM 区域
                # Topology 格式：方向（i/o）+ 起点-终点
                pattern = re.compile(r'[io](\d+)-(\d+)')
                matches = pattern.findall(topology_field)
                for start_s, end_s in matches:
                    start = _coerce_int(start_s)
                    end = _coerce_int(end_s)
                    if start is not None and end is not None:
                        tm_map.setdefault(protein_id, []).append(TMRegion(start=start, end=end))
            elif len(parts) >= 4:
                # 尝试短格式：ProteinID  TMhelix  start  end
                feature, start_s, end_s = parts[1], parts[2], parts[3]
                if feature.lower() in {"tmhelix", "tm"}:
                    start = _coerce_int(start_s)
                    end = _coerce_int(end_s)
                    if start is not None and end is not None:
                        tm_map.setdefault(protein_id, []).append(TMRegion(start=start, end=end))
    
    for regions in tm_map.values():
        regions.sort(key=lambda r: r.start)
    return tm_map


def format_positions(regions: List[Tuple[int, int]]) -> str:
    return ";".join(f"{s}-{e}" for s, e in regions)


MOTIF_PATTERNS = {
    "walker_a": re.compile(r"G.{4}GK[ST]"),  # A/GxxxxGKS/T
    "walker_b": re.compile(r"[ILVFM]{3,5}DE"),
    "signature": re.compile(r"LSGGQ"),
}


def _check_motifs(
    protein_id: str,
    seq: Optional[str],
    domains: List[Domain],
    window: int,
) -> Tuple[bool, str]:
    if seq is None:
        return False, f"{protein_id}:No sequence"
    misses = []
    for idx, dom in enumerate(domains, 1):
        # 扩展窗口，避免 HMM 对齐截断 Walker A / C-loop
        start = max(0, dom.start - 1 - window)
        end = min(len(seq), dom.end + window)
        segment = seq[start:end]
        has_a = bool(MOTIF_PATTERNS["walker_a"].search(segment))
        has_sig = bool(MOTIF_PATTERNS["signature"].search(segment))
        has_b = bool(MOTIF_PATTERNS["walker_b"].search(segment))
        if not (has_a and has_sig and has_b):
            lacks = []
            if not has_a:
                lacks.append("WalkerA")
            if not has_sig:
                lacks.append("LSGGQ")
            if not has_b:
                lacks.append("WalkerB")
            misses.append(f"NBD{idx}缺失:" + ",".join(lacks))
    if misses:
        return False, ";".join(misses)
    return True, "Motifs_OK"


def parse_neighbor_products(neighbor_path: Path) -> Dict[str, str]:
    """
    解析邻域注释文件，需包含 ProteinID 列和 NeighborProducts 列（分号分隔描述或关键词）。
    """
    if not neighbor_path.exists():
        return {}
    neighbor_map: Dict[str, str] = {}
    with neighbor_path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{neighbor_path} 没有表头")
        if "ProteinID" not in reader.fieldnames:
            raise ValueError(f"{neighbor_path} 缺少 ProteinID 列")
        # 尝试多种列名兼容
        product_field: Optional[str] = None
        for name in ["NeighborProducts", "neighbor_products", "Neighbors", "neighbors"]:
            if name in reader.fieldnames:
                product_field = name
                break
        if product_field is None:
            raise ValueError(f"{neighbor_path} 缺少 NeighborProducts/Neighbors 列")

        for row in reader:
            pid = row.get("ProteinID", "")
            products = row.get(product_field, "")
            if pid:
                neighbor_map[pid] = products.lower()
    return neighbor_map


def _check_neighborhood(
    pid: str,
    neighbor_map: Dict[str, str],
    keywords: List[str],
) -> Tuple[bool, str]:
    if not neighbor_map or not keywords:
        return True, "neighbor_skip"
    if pid not in neighbor_map:
        return False, "No_neighbor_info"
    text = neighbor_map[pid]
    if any(kw.lower() in text for kw in keywords):
        return True, "Neighbor_OK"
    return False, "Neighbor_no_keyword"


def evaluate_protein(
    domains: List[Domain],
    tms: List[TMRegion],
    seq: Optional[str],
    require_motif: bool,
    motif_window: int,
    pid: str,
    neighbor_map: Dict[str, str],
    neighbor_keywords: List[str],
) -> Tuple[bool, str, str, str]:
    """
    针对 Rv0194 / ABCB 型四域融合蛋白的宽松筛选逻辑：
    1) 恰好 2 个 PF00005 域（两个 NBD）
    2) TM 总数 >= 6（至少是多跨膜蛋白）
    3) 至少 1 条 TM 在 NBD1 前（对应 TMD1）
    4) 至少 1 条 TM 在 NBD1 与 NBD2 之间（对应 TMD2）
    不再要求 NBD2 之后有 TM。
    """
    # 1) NBD 数量
    if len(domains) != 2:
        return False, f"PF00005 count = {len(domains)}", "motif_skip", "neighbor_skip"

    # 2) TM 总数
    if len(tms) < 6:
        return False, f"TM count < 6 ({len(tms)})", "motif_skip", "neighbor_skip"

    # 按起始位置排序，确保 nbd1 在前、nbd2 在后
    nbd1, nbd2 = sorted(domains, key=lambda d: d.start)

    # 3) 至少一条 TM 在 NBD1 之前
    tm_before_nbd1 = [tm for tm in tms if tm.end <= nbd1.start]

    # 4) 至少一条 TM 在 NBD1 与 NBD2 之间
    tm_between = [tm for tm in tms if tm.start >= nbd1.end and tm.end <= nbd2.start]

    if not tm_before_nbd1:
        return False, "No TM before NBD1", "motif_skip", "neighbor_skip"

    if not tm_between:
        return False, "No TM between NBD1 and NBD2", "motif_skip", "neighbor_skip"

    motif_reason = "motif_skip"
    neighbor_reason = "neighbor_skip"
    if require_motif:
        motifs_ok, motif_reason = _check_motifs(pid, seq, domains, motif_window)
        if not motifs_ok:
            return False, "Motif_fail", motif_reason, neighbor_reason

    neighbor_ok, neighbor_reason = _check_neighborhood(pid, neighbor_map, neighbor_keywords)
    if not neighbor_ok:
        return False, "Neighbor_fail", motif_reason, neighbor_reason

    return True, "Pass", motif_reason, neighbor_reason

def run(
    hmmscan_csv: Path,
    tmhmm_path: Path,
    fasta_path: Optional[Path],
    output_dir: Path,
    require_motif: bool,
    motif_window: int,
    neighbor_path: Optional[Path],
    neighbor_keywords: List[str],
) -> None:
    # 保留旧接口：单个 FASTA/单次输出
    domains_map = parse_hmmscan_csv(hmmscan_csv) if hmmscan_csv.exists() else {}
    tm_map = parse_tmhmm(tmhmm_path) if tmhmm_path.exists() else {}
    seqs = load_fasta_seqs(fasta_path) if (fasta_path and fasta_path.exists()) else {}
    neighbor_map = parse_neighbor_products(neighbor_path) if neighbor_path else {}
    protein_ids = sorted(set(domains_map.keys()) | set(tm_map.keys()) | set(seqs.keys()))
    write_reports(
        protein_ids=protein_ids,
        domains_map=domains_map,
        tm_map=tm_map,
        seqs=seqs,
        neighbor_map=neighbor_map,
        neighbor_keywords=neighbor_keywords,
        output_dir=output_dir,
        require_motif=require_motif,
        motif_window=motif_window,
    )


def write_reports(
    protein_ids: Iterable[str],
    domains_map: Dict[str, List[Domain]],
    tm_map: Dict[str, List[TMRegion]],
    seqs: Dict[str, str],
    neighbor_map: Dict[str, str],
    neighbor_keywords: List[str],
    output_dir: Path,
    require_motif: bool,
    motif_window: int,
) -> List[str]:
    """
    核心写文件逻辑，可被单样本或批量模式调用。
    返回通过筛选的蛋白 ID 列表。
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    abc_summary_path = output_dir / "abc_tran_summary.tsv"
    tm_summary_path = output_dir / "tm_summary.tsv"
    candidates_path = output_dir / "rv0194_like_candidates.txt"
    report_path = output_dir / "detailed_report.tsv"

    candidates: List[str] = []

    with abc_summary_path.open("w", encoding="utf-8", newline="") as abc_out, \
            tm_summary_path.open("w", encoding="utf-8", newline="") as tm_out, \
            report_path.open("w", encoding="utf-8", newline="") as report_out:
        abc_writer = csv.writer(abc_out, delimiter="\t")
        tm_writer = csv.writer(tm_out, delimiter="\t")
        report_writer = csv.writer(report_out, delimiter="\t")

        abc_writer.writerow(["ProteinID", "ABC_tran_count", "Positions"])
        tm_writer.writerow(["ProteinID", "TM_count", "TM_positions"])
        report_writer.writerow(
            [
                "ProteinID",
                "ABC_tran_count",
                "ABC_positions",
                "TM_count",
                "TM_positions",
                "pass_filter",
                "reason",
                "motif_reason",
                "neighbor_reason",
            ]
        )

        for pid in protein_ids:
            domains = domains_map.get(pid, [])
            tms = tm_map.get(pid, [])
            abc_positions = format_positions([(d.start, d.end) for d in domains]) if domains else ""
            tm_positions = format_positions([(tm.start, tm.end) for tm in tms]) if tms else ""

            seq = seqs.get(pid)
            passed, reason, motif_reason, neighbor_reason = evaluate_protein(
                domains,
                tms,
                seq=seq,
                require_motif=require_motif,
                motif_window=motif_window,
                pid=pid,
                neighbor_map=neighbor_map,
                neighbor_keywords=neighbor_keywords,
            )
            if passed:
                candidates.append(pid)

            abc_writer.writerow([pid, len(domains), abc_positions])
            tm_writer.writerow([pid, len(tms), tm_positions])
            report_writer.writerow(
                [
                    pid,
                    len(domains),
                    abc_positions,
                    len(tms),
                    tm_positions,
                    str(passed),
                    reason,
                    motif_reason,
                    neighbor_reason,
                ]
            )

    with candidates_path.open("w", encoding="utf-8") as f:
        for pid in candidates:
            f.write(f"{pid}\n")

    return candidates


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="基于 hmmscan (Pfam) CSV 与 TMHMM 结果筛选 Rv0194 型四域融合结构。"
    )
    parser.add_argument("--hmmscan-csv", required=True, type=Path, help="hmmscan 转换后的 CSV 文件路径")
    parser.add_argument("--tmhmm", required=True, type=Path, help="TMHMM 输出文件路径")
    fasta_group = parser.add_mutually_exclusive_group()
    fasta_group.add_argument(
        "--fasta",
        required=False,
        type=Path,
        help="单个 FASTA 文件（提供后可开启基序检查）",
    )
    fasta_group.add_argument(
        "--fasta-dir",
        type=Path,
        help="包含多个样本 FASTA 的目录（每个文件独立输出到子目录）",
    )
    parser.add_argument(
        "--fasta-glob",
        type=str,
        default="*.faa,*.fa,*.fasta",
        help="匹配 --fasta-dir 下文件的 glob，逗号分隔（默认: *.faa,*.fa,*.fasta）",
    )
    parser.add_argument("--output-dir", required=True, type=Path, help="结果输出目录")
    parser.add_argument(
        "--require-motif",
        action="store_true",
        help="开启 WalkerA/WalkerB/LSGGQ 保守基序检查（需要提供 --fasta）",
    )
    parser.add_argument(
        "--motif-window",
        type=int,
        default=100,
        help="基序检查时在 HMM 对齐范围两侧扩展的氨基酸数 (default: 35)",
    )
    parser.add_argument(
        "--neighbor-tsv",
        type=Path,
        help="基因邻域注释 TSV（需含 ProteinID 与 NeighborProducts/Neighbors 列，分号分隔）",
    )
    parser.add_argument(
        "--neighbor-keywords",
        type=str,
        default="abc transporter;permease;transporter",
        help="邻域关键词，分号分隔（默认: abc transporter;permease;transporter）",
    )
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    neighbor_keywords = [kw.strip() for kw in args.neighbor_keywords.split(";") if kw.strip()]

    # 批量模式：按目录中的 FASTA 逐个输出到子目录
    if args.fasta_dir:
        patterns = [p.strip() for p in args.fasta_glob.split(",") if p.strip()]
        fasta_files = []
        for pat in patterns:
            fasta_files.extend(sorted(args.fasta_dir.glob(pat)))
        if not fasta_files:
            raise SystemExit(f"未在 {args.fasta_dir} 下找到匹配 {patterns} 的 FASTA")

        domains_map = parse_hmmscan_csv(args.hmmscan_csv) if args.hmmscan_csv.exists() else {}
        tm_map = parse_tmhmm(args.tmhmm) if args.tmhmm.exists() else {}
        neighbor_map = parse_neighbor_products(args.neighbor_tsv) if args.neighbor_tsv else {}

        for fasta_file in fasta_files:
            seqs = load_fasta_seqs(fasta_file)
            protein_ids = sorted(set(seqs.keys()))
            out_dir = args.output_dir / fasta_file.stem
            print(f"[INFO] 处理 {fasta_file.name}，蛋白数：{len(protein_ids)}，输出目录：{out_dir}")
            write_reports(
                protein_ids=protein_ids,
                domains_map=domains_map,
                tm_map=tm_map,
                seqs=seqs,
                neighbor_map=neighbor_map,
                neighbor_keywords=neighbor_keywords,
                output_dir=out_dir,
                require_motif=args.require_motif,
                motif_window=args.motif_window,
            )
    else:
        run(
            hmmscan_csv=args.hmmscan_csv,
            tmhmm_path=args.tmhmm,
            fasta_path=args.fasta,
            output_dir=args.output_dir,
            require_motif=args.require_motif,
            motif_window=args.motif_window,
            neighbor_path=args.neighbor_tsv,
            neighbor_keywords=neighbor_keywords,
        )


if __name__ == "__main__":
    main()
