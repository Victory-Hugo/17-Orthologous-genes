#!/usr/bin/env python3
"""
基于 hmmscan (Pfam) CSV 与 TMHMM 结果筛选 Rv0194 型四域融合结构 (TMD–NBD–TMD–NBD)。
"""

from __future__ import annotations

import argparse
import csv
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
    # 优先使用明确的 pfam_id/target_accession，其次 target_name。
    preferred = ["pfam_id", "target_accession", "target_name"]
    for name in preferred:
        if name in fieldnames:
            return name
    raise ValueError("hmmscan CSV 缺少 pfam_id/target_accession/target_name 列")


def parse_hmmscan_csv(csv_path: Path) -> Dict[str, List[Domain]]:
    """
    解析 hmmscan CSV，返回 protein_id -> ABC_tran 域列表。
    仅保留 pfam_id 以 PF00005 开头的行。
    """
    protein_to_domains: Dict[str, List[Domain]] = {}
    with csv_path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"{csv_path} 没有表头")
        pfam_field = _select_pfam_field(reader.fieldnames)
        name_field = "query_name" if "query_name" in reader.fieldnames else "target_name"
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


def parse_tmhmm(tmhmm_path: Path) -> Dict[str, List[TMRegion]]:
    """
    解析 TMHMM 标准输出：ProteinID  TMhelix  start  end
    跳过非 TMhelix 行和注释行。
    """
    tm_map: Dict[str, List[TMRegion]] = {}
    with tmhmm_path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            protein_id, feature, start_s, end_s = parts[0], parts[1], parts[2], parts[3]
            if feature.lower() not in {"tmhelix", "tm"}:
                continue
            start = _coerce_int(start_s)
            end = _coerce_int(end_s)
            if start is None or end is None:
                continue
            tm_map.setdefault(protein_id, []).append(TMRegion(start=start, end=end))
    for regions in tm_map.values():
        regions.sort(key=lambda r: r.start)
    return tm_map


def format_positions(regions: List[Tuple[int, int]]) -> str:
    return ";".join(f"{s}-{e}" for s, e in regions)


def evaluate_protein(domains: List[Domain], tms: List[TMRegion]) -> Tuple[bool, str]:
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
        return False, f"PF00005 count = {len(domains)}"

    # 2) TM 总数
    if len(tms) < 6:
        return False, f"TM count < 6 ({len(tms)})"

    # 按起始位置排序，确保 nbd1 在前、nbd2 在后
    nbd1, nbd2 = sorted(domains, key=lambda d: d.start)

    # 3) 至少一条 TM 在 NBD1 之前
    tm_before_nbd1 = [tm for tm in tms if tm.end <= nbd1.start]

    # 4) 至少一条 TM 在 NBD1 与 NBD2 之间
    tm_between = [tm for tm in tms if tm.start >= nbd1.end and tm.end <= nbd2.start]

    if not tm_before_nbd1:
        return False, "No TM before NBD1"

    if not tm_between:
        return False, "No TM between NBD1 and NBD2"

    return True, "Pass"

def run(
    hmmscan_csv: Path,
    tmhmm_path: Path,
    fasta_path: Optional[Path],
    output_dir: Path,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    domains_map = parse_hmmscan_csv(hmmscan_csv) if hmmscan_csv.exists() else {}
    tm_map = parse_tmhmm(tmhmm_path) if tmhmm_path.exists() else {}

    protein_ids = sorted(set(domains_map.keys()) | set(tm_map.keys()))

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
            ]
        )

        for pid in protein_ids:
            domains = domains_map.get(pid, [])
            tms = tm_map.get(pid, [])
            abc_positions = format_positions([(d.start, d.end) for d in domains]) if domains else ""
            tm_positions = format_positions([(tm.start, tm.end) for tm in tms]) if tms else ""

            passed, reason = evaluate_protein(domains, tms)
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
                ]
            )

    with candidates_path.open("w", encoding="utf-8") as f:
        for pid in candidates:
            f.write(f"{pid}\n")

    # fasta_path 保留参数接口，如需后续处理可在此扩展；当前逻辑未使用该文件。


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="基于 hmmscan (Pfam) CSV 与 TMHMM 结果筛选 Rv0194 型四域融合结构。"
    )
    parser.add_argument("--hmmscan-csv", required=True, type=Path, help="hmmscan 转换后的 CSV 文件路径")
    parser.add_argument("--tmhmm", required=True, type=Path, help="TMHMM 输出文件路径")
    parser.add_argument(
        "--fasta",
        required=False,
        type=Path,
        help="候选蛋白序列 fasta（当前未使用，保留接口）",
    )
    parser.add_argument("--output-dir", required=True, type=Path, help="结果输出目录")
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    run(
        hmmscan_csv=args.hmmscan_csv,
        tmhmm_path=args.tmhmm,
        fasta_path=args.fasta,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
