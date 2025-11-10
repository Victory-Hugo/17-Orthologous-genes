#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import os
import sys
import argparse
from collections import OrderedDict, defaultdict
from pathlib import Path


def load_hits_for_sample(tbl_file, sample_name):
    """
    从合并后的 tbl 文件中提取指定样本、按 query name 分类的命中蛋白 ID。
    返回 OrderedDict[str, list[str]]，保持 query name 与 target 的原始顺序。
    """
    hits_by_query: OrderedDict[str, list[str]] = OrderedDict()
    seen_per_query: defaultdict[str, set[str]] = defaultdict(set)

    try:
        with open(tbl_file, 'r', encoding='utf-8', errors='ignore', newline='') as handle:
            reader = csv.reader(handle, delimiter='\t')
            header = next(reader, None)
            if header is None:
                return hits

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

                    if target_name and target_name not in seen_per_query[query_name]:
                        hits_by_query.setdefault(query_name, []).append(target_name)
                        seen_per_query[query_name].add(target_name)

    except FileNotFoundError:
        raise
    except Exception as exc:
        raise RuntimeError(f"解析 {tbl_file} 失败: {exc}") from exc

    return hits_by_query


def extract_sequences_from_faa(faa_file, target_ids):
    """
    在 FAA FASTA 中提取目标蛋白序列，返回 {target_id: sequence}。
    """
    sequences = {}
    targets = set(target_ids)

    try:
        with open(faa_file, 'r', encoding='utf-8', errors='ignore') as handle:
            current_id = None
            current_seq = []

            for line in handle:
                line = line.rstrip('\n')
                if not line:
                    continue
                if line.startswith('>'):
                    if current_id in targets and current_seq:
                        sequences[current_id] = ''.join(current_seq)

                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)

            if current_id in targets and current_seq:
                sequences[current_id] = ''.join(current_seq)

    except Exception as exc:
        raise RuntimeError(f"读取 FAA 文件失败 - {faa_file}: {exc}") from exc

    return sequences


def format_fasta_sequence(sequence, width=60):
    for i in range(0, len(sequence), width):
        yield sequence[i:i + width]


def extract_hit_sequences(faa_file, tbl_file, output_dir):
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

    ordered_hits = [hit for hits in hits_by_query.values() for hit in hits]
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

    total_written = 0
    for query_name, query_hits in hits_by_query.items():
        if not query_hits:
            continue

        query_dir = os.path.join(output_dir, query_name)
        os.makedirs(query_dir, exist_ok=True)
        output_file = os.path.join(query_dir, f"{filename_base}.hits.faa")

        try:
            with open(output_file, 'w', encoding='utf-8') as handle:
                count = 0
                for hit in query_hits:
                    sequence = sequences.get(hit)
                    if not sequence:
                        continue
                    count += 1
                    header = f"{filename_base}_{count}"
                    handle.write(f">{header}\n")
                    for chunk in format_fasta_sequence(sequence):
                        handle.write(chunk + '\n')

            if count == 0:
                print(f"⚠️  警告: {filename_base} 在 {query_name} 中未写入任何序列")
            else:
                total_written += count
                print(f"✓ 成功: {filename_base} - {query_name} -> 输出 {count} 条序列")

        except Exception as exc:
            print(f"错误: 写入输出文件失败 - {output_file}: {exc}")

    if total_written == 0:
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

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    success = extract_hit_sequences(
        args.faa_file,
        args.tbl_file,
        args.output_dir,
    )

    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
