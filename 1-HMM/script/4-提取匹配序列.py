#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import math
import os
import sys
from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional


@dataclass
class HitRecord:
    query_name: str
    target_name: str
    evalue: float
    score: float


def parse_float(value: str, default: float) -> float:
    """安全地解析浮点数"""
    if value is None:
        return default
    value = value.strip()
    if not value or value == '-':
        return default
    try:
        return float(value)
    except ValueError:
        return default


def load_hits_for_sample(tbl_file: str, sample_name: str) -> "OrderedDict[str, List[HitRecord]]":
    """
    从合并后的 tbl 文件中提取指定样本、按 query name 分类的命中详情。
    返回 OrderedDict[str, list[HitRecord]]，保持 query 和 target 的原始顺序。
    """
    hits_by_query: "OrderedDict[str, List[HitRecord]]" = OrderedDict()
    seen_per_query: "defaultdict[str, set[str]]" = defaultdict(set)

    try:
        with open(tbl_file, 'r', encoding='utf-8', errors='ignore', newline='') as handle:
            reader = csv.reader(handle, delimiter='\t')
            header = next(reader, None)
            if header is None:
                return hits_by_query

            # 尝试定位样本和 target 列
            try:
                sample_idx = header.index('文件名')
            except ValueError:
                sample_idx = 0

            try:
                target_idx = header.index('target name')
            except ValueError:
                target_idx = 1 if len(header) > 1 else 0

            try:
                query_idx = header.index('query name')
            except ValueError:
                query_idx = 3 if len(header) > 3 else target_idx

            try:
                evalue_idx = header.index('E-value')
            except ValueError:
                evalue_idx = -1

            try:
                score_idx = header.index('score')
            except ValueError:
                score_idx = -1

            header_sample_value = header[sample_idx] if sample_idx < len(header) else ''
            header_target_value = header[target_idx] if target_idx < len(header) else ''

            for row in reader:
                if not row:
                    continue
                if row[0].startswith('#'):
                    continue
                # 如果合并时重复了表头，则跳过
                if sample_idx < len(row) and target_idx < len(row):
                    if row[sample_idx] == header_sample_value and row[target_idx] == header_target_value:
                        continue
                else:
                    continue

                if row[sample_idx] == sample_name:
                    query_name = row[query_idx] if query_idx < len(row) else ''
                    query_name = query_name.strip() or 'unknown_query'
                    target_name = row[target_idx]
                    evalue = parse_float(row[evalue_idx], float('inf')) if evalue_idx != -1 else float('inf')
                    score = parse_float(row[score_idx], float('-inf')) if score_idx != -1 else float('-inf')

                    if target_name and target_name not in seen_per_query[query_name]:
                        hits_by_query.setdefault(query_name, []).append(
                            HitRecord(
                                query_name=query_name,
                                target_name=target_name,
                                evalue=evalue,
                                score=score,
                            )
                        )
                        seen_per_query[query_name].add(target_name)

    except FileNotFoundError:
        raise
    except Exception as exc:
        raise RuntimeError(f"解析 {tbl_file} 失败: {exc}") from exc

    return hits_by_query


def extract_sequences_from_faa(faa_file, target_ids):
    """
    在 FAA FASTA 中提取目标蛋白序列，返回 {target_id: {'sequence': str, 'description': str}}。
    """
    sequences = {}
    targets = set(target_ids)

    try:
        with open(faa_file, 'r', encoding='utf-8', errors='ignore') as handle:
            current_id = None
            current_seq = []
            current_desc = ""

            for line in handle:
                line = line.rstrip('\n')
                if not line:
                    continue
                if line.startswith('>'):
                    if current_id in targets and current_seq:
                        sequences[current_id] = {
                            'sequence': ''.join(current_seq),
                            'description': current_desc,
                        }

                    header_line = line[1:].strip()
                    parts = header_line.split(maxsplit=1)
                    current_id = parts[0]
                    current_desc = parts[1] if len(parts) > 1 else ""
                    current_seq = []
                else:
                    current_seq.append(line)

            if current_id in targets and current_seq:
                sequences[current_id] = {
                    'sequence': ''.join(current_seq),
                    'description': current_desc,
                }

    except Exception as exc:
        raise RuntimeError(f"读取 FAA 文件失败 - {faa_file}: {exc}") from exc

    return sequences


def format_fasta_sequence(sequence, width=60):
    for i in range(0, len(sequence), width):
        yield sequence[i:i + width]


def extract_hit_sequences(
    faa_file: str,
    tbl_file: str,
    output_dir: str,
    max_evalue: Optional[float] = None,
    min_coverage: Optional[float] = None,
    min_length: Optional[int] = None,
    mode: str = 'all',
    model_length: Optional[float] = None,
):
    """
    从 FAA 文件中提取所有命中序列，每个命中以 >sample_编号 的形式输出。
    """
    if not os.path.isfile(faa_file):
        print(f"警告: FAA 文件不存在 - {faa_file}")
        return False

    if not os.path.isfile(tbl_file):
        print(f"警告: tbl 文件不存在 - {tbl_file}")
        return False

    filename_base = Path(faa_file).stem

    try:
        hits_by_query = load_hits_for_sample(tbl_file, filename_base)
    except FileNotFoundError:
        print(f"警告: tbl 文件不存在 - {tbl_file}")
        return False
    except RuntimeError as exc:
        print(exc)
        return False

    if not hits_by_query:
        print(f"❌ 无匹配结果: {filename_base}")
        return False

    ordered_hits = [hit.target_name for hits in hits_by_query.values() for hit in hits]
    if not ordered_hits:
        print(f"❌ 无匹配结果: {filename_base}")
        return False

    try:
        sequences = extract_sequences_from_faa(faa_file, ordered_hits)
    except RuntimeError as exc:
        print(exc)
        return False

    if not sequences:
        print(f"⚠️  警告: 在 FAA 文件中未找到任何命中序列 - {filename_base}")
        return False

    missing_overall = [hit for hit in ordered_hits if hit not in sequences]
    if missing_overall:
        print(f"⚠️  警告: {filename_base} 缺失 {len(missing_overall)} 条序列 (示例: {', '.join(missing_overall[:3])})")

    filter_summary: Dict[str, int] = {
        'total': 0,
        'passed': 0,
        'evalue': 0,
        'length': 0,
        'coverage': 0,
    }

    coverage_enabled = min_coverage is not None
    coverage_available = coverage_enabled and model_length and model_length > 0
    if coverage_enabled and not coverage_available:
        print(f"⚠️  覆盖度阈值 {min_coverage} 已设置，但未提供模型长度，已忽略覆盖度过滤。")
        coverage_enabled = False

    total_written = 0
    for query_name, query_hits in hits_by_query.items():
        if not query_hits:
            continue

        # 输出目录：query 名 + ".aln"（若尚未带后缀），并放入 sequence 子目录
        query_dir_name = query_name if query_name.endswith(".aln") else f"{query_name}.aln"
        query_dir = os.path.join(output_dir, query_dir_name)
        sequence_dir = os.path.join(query_dir, "sequence")
        os.makedirs(sequence_dir, exist_ok=True)
        output_file = os.path.join(sequence_dir, f"{filename_base}.hits.faa")

        try:
            with open(output_file, 'w', encoding='utf-8') as handle:
                count = 0
                filtered_hits = []

                for hit in query_hits:
                    sequence_entry = sequences.get(hit.target_name)
                    if not sequence_entry:
                        continue

                    sequence = sequence_entry.get('sequence', '')
                    description = sequence_entry.get('description', '')
                    desc_clean = description.replace('|', '/').replace('\t', ' ').strip()

                    filter_summary['total'] += 1

                    if max_evalue is not None and hit.evalue > max_evalue:
                        filter_summary['evalue'] += 1
                        continue

                    length = len(sequence)
                    if min_length is not None and length < min_length:
                        filter_summary['length'] += 1
                        continue

                    coverage = None
                    if coverage_enabled and coverage_available:
                        coverage = length / float(model_length)  # type: ignore[arg-type]
                        if coverage < float(min_coverage):
                            filter_summary['coverage'] += 1
                            continue

                    filtered_hits.append({
                        'hit': hit,
                        'sequence': sequence,
                        'length': length,
                        'coverage': coverage,
                        'description': desc_clean,
                    })
                    filter_summary['passed'] += 1

                if mode == 'only_one' and filtered_hits:
                    def _rank_key(item):
                        score = item['hit'].score if math.isfinite(item['hit'].score) else float('-inf')
                        evalue_component = -item['hit'].evalue if math.isfinite(item['hit'].evalue) else float('-inf')
                        return (
                            score,
                            evalue_component,
                            item['length'],
                        )

                    best_hit = max(filtered_hits, key=_rank_key)
                    filtered_hits = [best_hit]

                for idx, packaged in enumerate(filtered_hits, start=1):
                    count += 1
                    header_prefix = filename_base if mode == 'only_one' else f"{filename_base}_{idx}"
                    meta_parts = [
                        f"id={packaged['hit'].target_name}",
                        f"score={packaged['hit'].score:.2f}" if math.isfinite(packaged['hit'].score) else "score=NA",
                        f"evalue={packaged['hit'].evalue:.2e}" if math.isfinite(packaged['hit'].evalue) else "evalue=NA",
                        f"len={packaged['length']}",
                    ]
                    if packaged.get('coverage') is not None:
                        meta_parts.append(f"cov={packaged['coverage']:.3f}")
                    if packaged.get('description'):
                        meta_parts.append(f"desc={packaged['description']}")
                    header = f"{header_prefix}|" + "|".join(meta_parts)
                    handle.write(f">{header}\n")
                    for chunk in format_fasta_sequence(packaged['sequence']):
                        handle.write(chunk + '\n')

            if count == 0:
                print(f"⚠️  警告: {filename_base} 在 {query_name} 中未写入任何序列（全部命中未通过过滤）")
            else:
                total_written += count
                print(f"✓ 成功: {filename_base} - {query_name} -> 输出 {count} 条序列 (模式: {mode})")

        except Exception as exc:
            print(f"错误: 写入输出文件失败 - {output_file}: {exc}")

    if total_written == 0:
        print(f"❌ {filename_base} 未通过过滤，命中数: {filter_summary['total']}, "
              f"E-value 过滤: {filter_summary['evalue']}, "
              f"长度过滤: {filter_summary['length']}, "
              f"覆盖度过滤: {filter_summary['coverage']}, "
              f"最终通过: {filter_summary['passed']}")
        return False

    return True


def main():
    parser = argparse.ArgumentParser(
        description='从 HMM 合并结果中提取指定样本的全部命中序列',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例用法:
  # 处理单个文件
  python3 4-提取最佳匹配序列.py /path/to/sample.faa /path/to/tbl_merge.tsv /path/to/output
  
  # 使用 parallel 处理多个文件（假设 faa_files.txt 包含所有文件路径）
  cat faa_files.txt | parallel python3 4-提取最佳匹配序列.py {} /path/to/tbl_merge.tsv /path/to/output
        '''
    )

    parser.add_argument('faa_file', help='FAA 文件路径')
    parser.add_argument('tbl_file', help='合并的 tbl 结果文件路径')
    parser.add_argument('output_dir', help='输出目录')
    parser.add_argument('--max-evalue', type=float, default=None, help='允许的最大 E-value（全局命中）')
    parser.add_argument('--min-coverage', type=float, default=None, help='最小覆盖度（命中长度 / 模型长度）')
    parser.add_argument('--min-length', type=int, default=None, help='最小氨基酸长度')
    parser.add_argument('--mode', choices=['all', 'only_one'], default='all', help='all: 输出所有通过过滤的命中；only_one: 仅输出最高分命中')
    parser.add_argument('--model-length', type=float, default=None, help='HMM 模型长度（用于覆盖度计算）')

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    success = extract_hit_sequences(
        args.faa_file,
        args.tbl_file,
        args.output_dir,
        max_evalue=args.max_evalue,
        min_coverage=args.min_coverage,
        min_length=args.min_length,
        mode=args.mode,
        model_length=args.model_length,
    )

    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
