#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from pathlib import Path

def extract_best_hit_sequence(faa_file, results_dir, output_dir, results_prefix="hits_"):
    """
    从 FAA 文件中提取最佳匹配的序列
    
    参数:
        faa_file: FAA 文件路径
        results_dir: HMM 比对结果目录
        output_dir: 输出目录
        results_prefix: 结果文件前缀（默认 "hits_"）
    """
    
    # 检查文件是否存在
    if not os.path.isfile(faa_file):
        print(f"警告: FAA 文件不存在 - {faa_file}")
        return False
    
    # 获取文件名（不含扩展名）
    filename = os.path.basename(faa_file)
    filename_base = filename.replace('.faa', '')
    
    # 对应的结果文件
    result_file = os.path.join(results_dir, f"{results_prefix}{filename_base}.tbl")
    
    # 检查结果文件是否存在
    if not os.path.isfile(result_file):
        print(f"❌ 无匹配结果: {filename_base}")
        return False
    
    # 读取结果文件，获取最佳匹配的蛋白 ID（第一行就是最佳匹配）
    best_hit = None
    try:
        with open(result_file, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                # 跳过注释行和空行
                if line.startswith('#') or not line.strip():
                    continue
                # 第一列是蛋白 ID（target name）
                fields = line.split()
                if len(fields) > 0:
                    best_hit = fields[0]
                    break  # 只取第一行（最佳匹配）
    except Exception as e:
        print(f"错误: 读取结果文件失败 - {result_file}: {e}")
        return False
    
    # 如果没有找到任何匹配
    if best_hit is None:
        print(f"❌ 无匹配结果: {filename_base}")
        return False
    
    # 从 FAA 文件中提取最佳匹配的序列
    sequence = None
    try:
        with open(faa_file, 'r', encoding='utf-8', errors='ignore') as f:
            current_id = None
            current_seq = []
            
            for line in f:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    # 检查是否找到了目标序列
                    if current_id == best_hit and current_seq:
                        sequence = ''.join(current_seq)
                        break
                    
                    # 新的序列开始
                    # 提取序列 ID（第一列，空格分隔）
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # 检查最后一个序列
            if current_id == best_hit and current_seq:
                sequence = ''.join(current_seq)
    
    except Exception as e:
        print(f"错误: 读取 FAA 文件失败 - {faa_file}: {e}")
        return False
    
    # 如果没有找到序列
    if sequence is None:
        print(f"⚠️  警告: 在 FAA 文件中未找到蛋白 {best_hit} - {filename_base}")
        return False
    
    # 写入输出文件
    output_file = os.path.join(output_dir, f"{filename_base}.best_hit.faa")
    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            # 按要求格式化序列名称：使用原始文件名
            f.write(f">{filename_base}|lolD\n")
            # 按 60 个字符每行写入序列
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i+60] + '\n')
        
        print(f"✓ 成功: {filename_base} -> {best_hit}")
        return True
    
    except Exception as e:
        print(f"错误: 写入输出文件失败 - {output_file}: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description='从 HMM 比对结果中提取最佳匹配的蛋白序列',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例用法:
  # 处理单个文件
  python3 4-提取最佳匹配序列.py /path/to/file.faa /path/to/results /path/to/output
  
  # 使用 parallel 处理多个文件（假设 faa_files.txt 包含所有文件路径）
  cat faa_files.txt | parallel python3 4-提取最佳匹配序列.py {} /path/to/results /path/to/output
        '''
    )
    
    parser.add_argument('faa_file', help='FAA 文件路径')
    parser.add_argument('results_dir', help='HMM 比对结果目录')
    parser.add_argument('output_dir', help='输出目录')
    parser.add_argument('--results-prefix', default='hits_', 
                       help='结果文件前缀 (默认: hits_)')
    
    args = parser.parse_args()
    
    # 确保输出目录存在
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 处理文件
    success = extract_best_hit_sequence(
        args.faa_file,
        args.results_dir,
        args.output_dir,
        args.results_prefix
    )
    
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
