#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从命中蛋白序列 (.hits.faa) 中提取对应的核酸序列。

核心思路：
1. 读取命中 fasta，解析 header 中的 `id=` 字段作为目标 ID。
2. 按样本名匹配对应的 GFF 与基因组序列（支持 .fna / .fasta / .fa 等）。
3. 过滤出包含这些 ID 的 GFF 子集，调用 gffread -g genome -x 输出对应 CDS。

支持命令行与 import 调用，主流程入口为 run()。
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import subprocess
import tempfile
import threading
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple


# 解析 hits header 用的简单正则
ID_RE = re.compile(r"id=([^|]+)")


@dataclass
class HitEntry:
    ids: Set[str]
    raw_headers: List[str]


def parse_hit_ids(hits_faa: Path) -> HitEntry:
    ids: Set[str] = set()
    headers: List[str] = []
    with hits_faa.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line.startswith(">"):
                continue
            headers.append(line[1:])
            match = ID_RE.search(line)
            if match:
                ids.add(match.group(1))
            else:
                token = line[1:].split("|", 1)[0].strip()
                if token:
                    ids.add(token)
    return HitEntry(ids=ids, raw_headers=headers)


def choose_file(base_dir: Path, sample: str, extensions: Sequence[str]) -> Optional[Path]:
    for ext in extensions:
        candidate = base_dir / f"{sample}{ext}"
        if candidate.exists():
            return candidate
    return None


def load_gff_features(gff_path: Path):
    features = []
    id_to_children: Dict[str, Set[int]] = defaultdict(set)
    opener = gzip.open if gff_path.suffix.endswith("gz") else open
    with opener(gff_path, "rt", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            attr_field = parts[8]
            attrs = {}
            for kv in attr_field.split(";"):
                if not kv or "=" not in kv:
                    continue
                key, val = kv.split("=", 1)
                attrs[key.strip()] = val.strip()
            feature_id = attrs.get("ID")
            parents = set()
            if "Parent" in attrs:
                parents = set(p.strip() for p in attrs["Parent"].split(",") if p.strip())
            idx_in_list = len(features)
            features.append(
                {
                    "index": idx_in_list,
                    "line": line,
                    "id": feature_id,
                    "parents": parents,
                    "attrs": attrs,
                }
            )
            for parent in parents:
                id_to_children[parent].add(idx_in_list)
    return features, id_to_children


def collect_related_indices(
    features: List[Dict[str, object]],
    id_to_children: Dict[str, Set[int]],
    target_ids: Set[str],
) -> Set[int]:
    matched_ids: Set[str] = set()
    matched_indices: Set[int] = set()

    # 初步匹配：属性值中包含目标 ID
    for feature in features:
        attrs: Dict[str, str] = feature["attrs"]  # type: ignore[assignment]
        fid = feature["id"]  # type: ignore[assignment]
        all_values: Set[str] = set()
        for val in attrs.values():
            for token in val.split(","):
                token = token.strip()
                if not token:
                    continue
                all_values.add(token)
                # 兼容形如 cds-NP_xxx 或 rna-NP_xxx 的前缀
                if "-" in token:
                    all_values.add(token.split("-", 1)[1])
        if all_values & target_ids:
            if fid:
                matched_ids.add(fid)
            matched_indices.add(feature["index"])  # type: ignore[arg-type]

    # 扩展：包含所有父节点与子节点
    queue: List[str] = list(matched_ids)
    seen: Set[str] = set(queue)
    id_to_parents: Dict[str, Set[str]] = defaultdict(set)
    for feature in features:
        fid = feature["id"]  # type: ignore[assignment]
        if fid:
            id_to_parents[fid].update(feature["parents"])  # type: ignore[arg-type]

    while queue:
        current = queue.pop()
        # 子节点
        for child_idx in id_to_children.get(current, set()):
            matched_indices.add(child_idx)
            child_id = features[child_idx]["id"]  # type: ignore[index]
            if child_id and child_id not in seen:
                seen.add(child_id)
                queue.append(child_id)
        # 父节点
        for parent_id in id_to_parents.get(current, set()):
            if parent_id not in seen:
                seen.add(parent_id)
                queue.append(parent_id)

    # 包含父/子但无 ID 的行（例如部分 GFF 记录）
    for idx, feature in enumerate(features):
        if feature["id"] is None and feature["parents"] and feature["parents"] & seen:  # type: ignore[operator]
            matched_indices.add(idx)

    return matched_indices


def write_filtered_gff(features, indices: Set[int], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [features[i]["line"] for i in sorted(indices)]
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_gffread(gffread_bin: str, genome: Path, gff_subset: Path, output_fasta: Path) -> None:
    cmd = [
        gffread_bin,
        "-g",
        str(genome),
        "-x",
        str(output_fasta),
        str(gff_subset),
    ]
    subprocess.run(cmd, check=True)


def parse_hits_metadata(hits_faa: Path) -> Dict[str, Dict[str, str]]:
    """
    从 .hits.faa header 解析元信息，返回 {target_id: {score/evalue/len/...}}
    """
    meta_map: Dict[str, Dict[str, str]] = {}
    with hits_faa.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            fields = header.split("|")
            target_token = fields[0]
            target_id = target_token
            for field in fields:
                if field.startswith("id="):
                    target_id = field.split("=", 1)[1]
                    break
            meta: Dict[str, str] = {}
            for field in fields[1:]:
                if "=" not in field:
                    continue
                k, v = field.split("=", 1)
                meta[k.strip()] = v.strip()
            meta_map[target_id] = meta
    return meta_map


def reheader_fasta(
    fasta_path: Path,
    sample: str,
    query: str,
    meta_lookup: Dict[str, Dict[str, str]],
) -> None:
    """
    重新写入 fasta，header 增加 sample/query/id/score/evalue/len 等信息，序列自动按 60 列换行。
    """
    if not fasta_path.exists():
        return
    records: List[Tuple[str, str]] = []
    with fasta_path.open("r", encoding="utf-8", errors="ignore") as handle:
        header = None
        seq_lines: List[str] = []
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_lines)))
                header = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header is not None:
            records.append((header, "".join(seq_lines)))

    with fasta_path.open("w", encoding="utf-8") as handle:
        for orig_id, seq in records:
            meta = meta_lookup.get(orig_id, meta_lookup.get(orig_id.split(".")[0], {}))
            parts = [
                f"sample={sample}",
                f"query={query}",
                f"id={orig_id}",
            ]
            for key in ("score", "evalue", "len", "cov"):
                if key in meta and meta[key]:
                    parts.append(f"{key}={meta[key]}")
            header = "|".join(parts)
            handle.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                handle.write(seq[i : i + 60] + "\n")


def run(
    hits_dir: Path,
    gff_dir: Path,
    genome_dir: Path,
    gffread_bin: str = "gffread",
    output_dir: Optional[Path] = None,
    log_dir: Optional[Path] = None,
    max_workers: int = 1,
) -> None:
    hits_dir = hits_dir.resolve()
    gff_dir = gff_dir.resolve()
    genome_dir = genome_dir.resolve()
    out_dir = output_dir.resolve() if output_dir else hits_dir
    log_root = (log_dir or Path("/mnt/d/5-NCBI-Reference/hmm分析示例/log")).resolve()

    if not hits_dir.is_dir():
        raise FileNotFoundError(f"hits 目录不存在: {hits_dir}")
    if not gff_dir.is_dir():
        raise FileNotFoundError(f"GFF 目录不存在: {gff_dir}")
    if not genome_dir.is_dir():
        raise FileNotFoundError(f"基因组目录不存在: {genome_dir}")

    fasta_exts = [".fna", ".fasta", ".fa", ".fa.gz", ".fasta.gz", ".fna.gz"]
    gff_exts = [".gff", ".gff3", ".gff.gz", ".gff3.gz"]

    hits_files = sorted(hits_dir.glob("*.hits.faa"))
    if not hits_files:
        print(f"⚠️  未找到 .hits.faa: {hits_dir}")
        return

    log_root.mkdir(parents=True, exist_ok=True)
    text_log = log_root / "hits_cds_extract.log"
    summary_tsv = log_root / "hits_cds_extract_summary.tsv"
    summary_json = log_root / "hits_cds_extract_summary.json"

    gff_cache: Dict[Path, Tuple[List[Dict[str, object]], Dict[str, Set[int]]]] = {}
    gff_cache_lock = threading.Lock()

    query_name = hits_dir.parent.name if hits_dir.parent.name else "unknown_query"
    results: List[Dict[str, object]] = []

    def log_line(msg: str) -> None:
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        with text_log.open("a", encoding="utf-8") as handle:
            handle.write(f"[{timestamp}] {msg}\n")

    def process_one(hits_file: Path) -> Dict[str, object]:
        start = time.time()
        sample_stem = hits_file.stem  # e.g., GCF_xxx.hits
        sample_name = sample_stem[:-5] if sample_stem.endswith(".hits") else sample_stem
        result = {
            "sample": sample_name,
            "query": query_name,
            "status": "skip",
            "reason": "",
            "output": "",
            "matched_ids": 0,
            "matched_features": 0,
            "elapsed_sec": 0.0,
        }

        gff_path = choose_file(gff_dir, sample_name, gff_exts)
        genome_path = choose_file(genome_dir, sample_name, fasta_exts)

        if not gff_path or not genome_path:
            reason = f"缺少GFF或基因组 gff={bool(gff_path)} genome={bool(genome_path)}"
            log_line(f"[SKIP] {sample_name}: {reason}")
            result["reason"] = reason
            result["elapsed_sec"] = time.time() - start
            return result

        hit_entry = parse_hit_ids(hits_file)
        if not hit_entry.ids:
            reason = f"{hits_file.name} 未解析到 id= 字段"
            log_line(f"[SKIP] {sample_name}: {reason}; header示例: {hit_entry.raw_headers[:2]}")
            result["reason"] = reason
            result["elapsed_sec"] = time.time() - start
            return result

        with gff_cache_lock:
            cached = gff_cache.get(gff_path)
        if cached:
            features, id_to_children = cached
        else:
            features, id_to_children = load_gff_features(gff_path)
            with gff_cache_lock:
                gff_cache[gff_path] = (features, id_to_children)

        matched_indices = collect_related_indices(features, id_to_children, hit_entry.ids)
        if not matched_indices:
            reason = f"在 {gff_path.name} 中未找到匹配 ID"
            log_line(f"[SKIP] {sample_name}: {reason}; ids示例: {list(hit_entry.ids)[:3]}")
            result["reason"] = reason
            result["elapsed_sec"] = time.time() - start
            return result

        with tempfile.TemporaryDirectory(prefix="gff_subset_") as tmpdir:
            subset_path = Path(tmpdir) / f"{sample_name}.gff"
            write_filtered_gff(features, matched_indices, subset_path)
            output_fasta = out_dir / f"{sample_name}.hits.cds.fna"
            try:
                run_gffread(gffread_bin, genome_path, subset_path, output_fasta)
            except subprocess.CalledProcessError as exc:
                reason = f"gffread失败: {exc}"
                log_line(f"[ERROR] {sample_name}: {reason}")
                result["reason"] = reason
                result["elapsed_sec"] = time.time() - start
                return result

            meta_lookup = parse_hits_metadata(hits_file)
            reheader_fasta(output_fasta, sample_name, query_name, meta_lookup)
            elapsed = time.time() - start
            log_line(f"[OK] {sample_name}: 输出 {output_fasta} ({len(matched_indices)} 条GFF记录，耗时 {elapsed:.2f}s)")
            result.update(
                {
                    "status": "ok",
                    "reason": "",
                    "output": str(output_fasta),
                    "matched_ids": len(hit_entry.ids),
                    "matched_features": len(matched_indices),
                    "elapsed_sec": elapsed,
                }
            )
            return result

    max_workers = max(1, int(max_workers))
    if max_workers == 1:
        for hits_file in hits_files:
            results.append(process_one(hits_file))
    else:
        from concurrent.futures import ThreadPoolExecutor

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            for res in executor.map(process_one, hits_files):
                results.append(res)

    # 写 summary
    if results:
        with summary_tsv.open("w", encoding="utf-8") as handle:
            handle.write("sample\tquery\tstatus\treason\toutput\tmatched_ids\tmatched_features\telapsed_sec\n")
            for r in results:
                handle.write(
                    f"{r['sample']}\t{r['query']}\t{r['status']}\t{r['reason']}\t"
                    f"{r['output']}\t{r['matched_ids']}\t{r['matched_features']}\t{r['elapsed_sec']:.3f}\n"
                )
        summary_json.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"✓ 处理完成，日志: {text_log}，汇总: {summary_tsv}")


def main():
    parser = argparse.ArgumentParser(description="根据命中蛋白序列提取对应核酸 CDS（使用 gffread）")
    parser.add_argument("--hits-dir", required=True, help="命中蛋白序列目录（*.hits.faa）")
    parser.add_argument("--gff-dir", required=True, help="GFF 目录（文件名需与样本名匹配）")
    parser.add_argument("--genome-dir", required=True, help="基因组目录（支持 .fna/.fasta/.fa）")
    parser.add_argument("--gffread-bin", default="gffread", help="gffread 可执行文件路径/命令名")
    parser.add_argument("--output-dir", default=None, help="输出目录（默认同 hits-dir）")
    parser.add_argument("--log-dir", default=None, help="日志输出目录，默认 /mnt/d/5-NCBI-Reference/hmm分析示例/log")
    parser.add_argument("--max-workers", type=int, default=1, help="并发 gffread 数量，默认 1")
    args = parser.parse_args()

    run(
        hits_dir=Path(args.hits_dir),
        gff_dir=Path(args.gff_dir),
        genome_dir=Path(args.genome_dir),
        gffread_bin=args.gffread_bin,
        output_dir=Path(args.output_dir) if args.output_dir else None,
        log_dir=Path(args.log_dir) if args.log_dir else None,
        max_workers=args.max_workers,
    )


if __name__ == "__main__":
    main()
