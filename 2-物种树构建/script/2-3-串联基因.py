#!/usr/bin/env python3
# ============================================================================
# 脚本名称: 2-3-串联基因.py
# 作者: BigLin
# 功能: 串联多个对齐的蛋白质序列文件，生成用于物种树构建的超矩阵
# ============================================================================
# 
# 依赖软件和包:
#   - Python 3.6+
#   - 标准库: os, sys, argparse, pathlib, collections, signal, atexit, shutil
#
# 使用方式:
#   python3 2-3-串联基因.py
#
# ============================================================================

import os
import sys
import argparse
import signal
import atexit
import shutil
from pathlib import Path
from collections import defaultdict

# ============================================================================
# 全局配置变量 (所有变量硬编码)
# ============================================================================

# 输入/输出路径配置
ALIGNMENT_DIR = "/mnt/f/15_Bam_Tam/5-补齐更多物种/output/alignments"
OUTPUT_DIR = "/mnt/f/15_Bam_Tam/5-补齐更多物种/output/concatenated"
OUTPUT_FILE = "/mnt/f/15_Bam_Tam/5-补齐更多物种/output/concatenated/supermatrix.faa"
TEMP_FILES = []  # 用于记录临时文件

# 颜色定义
class Colors:
    RESET = '\033[0m'
    RED = '\033[91m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    WHITE = '\033[97m'
    BOLD = '\033[1m'

def print_header(msg):
    """打印彩色标题"""
    print(f"\n{Colors.BOLD}{Colors.CYAN}{'='*70}{Colors.RESET}")
    print(f"{Colors.BOLD}{Colors.BLUE}▶ {msg}{Colors.RESET}")
    print(f"{Colors.BOLD}{Colors.CYAN}{'='*70}{Colors.RESET}\n")

def print_success(msg):
    """打印成功信息"""
    print(f"{Colors.GREEN}✓ {msg}{Colors.RESET}")

def print_error(msg):
    """打印错误信息"""
    print(f"{Colors.RED}✗ {msg}{Colors.RESET}", file=sys.stderr)

def print_warning(msg):
    """打印警告信息"""
    print(f"{Colors.YELLOW}⚠ {msg}{Colors.RESET}")

def print_info(msg):
    """打印信息"""
    print(f"{Colors.BLUE}ℹ {msg}{Colors.RESET}")

def print_progress_bar(current, total, width=50):
    """打印进度条"""
    percent = current / total
    filled = int(width * percent)
    bar = '█' * filled + '░' * (width - filled)
    print(f"\r{Colors.CYAN}[{bar}] {current}/{total} ({percent*100:.1f}%){Colors.RESET}", end='', flush=True)

def cleanup_on_exit():
    """脚本退出时清理未完成的临时文件"""
    if TEMP_FILES:
        print_warning("\n脚本中断，清理未完成的临时文件...")
        for temp_file in TEMP_FILES:
            if os.path.exists(temp_file):
                try:
                    if os.path.isdir(temp_file):
                        shutil.rmtree(temp_file)
                    else:
                        os.remove(temp_file)
                    print_info(f"已删除临时文件: {temp_file}")
                except Exception as e:
                    print_error(f"删除文件失败 {temp_file}: {e}")

def signal_handler(signum, frame):
    """处理 Ctrl+C 信号"""
    print_error("\n用户中断执行...")
    cleanup_on_exit()
    sys.exit(1)

# 注册信号处理器和退出处理器
signal.signal(signal.SIGINT, signal_handler)
atexit.register(cleanup_on_exit)

def parse_fasta_files(alignment_dir):
    """
    解析alignment目录中的所有FASTA文件
    返回字典：{species_id: concatenated_sequence}
    同时返回物种基因覆盖信息
    """
    sequences = defaultdict(str)
    gene_files = sorted(Path(alignment_dir).glob("*_aligned.faa"))
    
    if not gene_files:
        print_error(f"在 {alignment_dir} 中未找到 *_aligned.faa 文件")
        return None, None
    
    print_success(f"找到 {len(gene_files)} 个基因文件\n")
    
    species_gene_count = defaultdict(int)  # 记录每个物种找到的基因数
    
    for idx, gene_file in enumerate(gene_files, 1):
        print_progress_bar(idx, len(gene_files))
        
        with open(gene_file, 'r') as f:
            current_species = None
            current_seq = ""
            
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                if line.startswith('>'):
                    # 保存前一个序列
                    if current_species:
                        sequences[current_species] += current_seq
                        species_gene_count[current_species] += 1
                    
                    # 解析新的头部行，提取物种ID (GCA_XXXXX.X)
                    header = line[1:]  # 去掉 '>'
                    # 只提取第一个字段作为物种ID
                    current_species = header.split()[0]
                    current_seq = ""
                else:
                    # 序列行
                    current_seq += line
            
            # 保存最后一个序列
            if current_species:
                sequences[current_species] += current_seq
                species_gene_count[current_species] += 1
    
    print()  # 换行
    
    print_info("物种基因覆盖情况：")
    # 检查每个物种覆盖的基因数
    gene_counts = {}
    for species, count in species_gene_count.items():
        gene_counts[count] = gene_counts.get(count, 0) + 1
    
    for count in sorted(gene_counts.keys()):
        print(f"  {Colors.CYAN}•{Colors.RESET} {gene_counts[count]} 个物种包含 {count} 个基因")
    
    return dict(sequences), species_gene_count

def write_concatenated_fasta(sequences, output_file, species_gene_count=None, min_genes=None):
    """
    将串联后的序列写入FASTA文件
    如果指定min_genes，只保留在至少min_genes个基因中出现的物种
    """
    if min_genes:
        if not species_gene_count:
            print_warning("未提供物种基因覆盖信息，忽略 min_genes 过滤条件")
        else:
            print_info(f"筛选在至少 {min_genes} 个基因中出现的物种...")
            sequences = {
                species_id: seq
                for species_id, seq in sequences.items()
                if species_gene_count.get(species_id, 0) >= min_genes
            }
    
    if not sequences:
        print_error("没有物种满足筛选条件")
        return False

    total_species = len(sequences)
    species_list = sorted(sequences.keys())
    
    print_info(f"写入 {total_species} 个物种的序列...\n")
    
    try:
        # 添加到临时文件列表
        TEMP_FILES.append(output_file)
        
        with open(output_file, 'w') as f:
            for idx, species_id in enumerate(species_list, 1):
                print_progress_bar(idx, total_species)
                
                f.write(f">{species_id}\n")
                # 按80个字符一行写入序列（符合FASTA标准）
                seq = sequences[species_id]
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + "\n")
        
        print()  # 换行
        print_success(f"成功生成串联FASTA文件：{output_file}")
        
        # 完成后从临时文件列表移除
        if output_file in TEMP_FILES:
            TEMP_FILES.remove(output_file)
        
        return True
    except Exception as e:
        print_error(f"写入文件失败: {e}")
        return False

def main():
    """主程序函数"""
    print_header("蛋白质序列串联工具")
    print_info(f"作者: BigLin\n")
    
    # 使用全局硬编码配置
    alignment_dir = ALIGNMENT_DIR
    output_dir = OUTPUT_DIR
    output_file = OUTPUT_FILE
    
    # 检查输入目录是否存在
    if not os.path.isdir(alignment_dir):
        print_error(f"输入目录不存在: {alignment_dir}")
        return False
    
    print_success(f"输入目录: {alignment_dir}")
    
    # 检查输出目录是否存在
    if not os.path.isdir(output_dir):
        print_error(f"输出目录不存在: {output_dir}")
        return False
    
    print_success(f"输出目录: {output_dir}")
    
    # 解析FASTA文件
    print_header("第一步: 读取基因对齐文件")
    print_info(f"从 {alignment_dir} 读取基因对齐文件...\n")
    sequences, species_gene_count = parse_fasta_files(alignment_dir)
    
    if not sequences:
        return False
    
    # 获取基因总数
    total_genes = len([f for f in Path(alignment_dir).glob("*_aligned.faa")])
    print_success(f"总基因数: {total_genes}")
    
    # 过滤：只保留在所有基因中都出现的物种
    print_header("第二步: 筛选物种")
    print_info(f"筛选在所有 {total_genes} 个基因中都出现的物种...\n")
    
    filtered_sequences = {}
    for species_id, seq in sequences.items():
        if species_gene_count[species_id] == total_genes:
            filtered_sequences[species_id] = seq
    
    print_success(f"原始物种数: {len(sequences)}")
    print_success(f"筛选后物种数: {len(filtered_sequences)}\n")
    
    if not filtered_sequences:
        print_error("没有物种在所有基因中都出现！")
        return False
    
    # 统计串联长度
    first_species = list(filtered_sequences.keys())[0]
    total_length = len(filtered_sequences[first_species])
    print_success(f"每个物种的串联序列长度: {total_length}")
    
    # 写入输出文件
    print_header("第三步: 生成串联序列文件")
    success = write_concatenated_fasta(
        filtered_sequences,
        output_file,
        species_gene_count=species_gene_count,
        min_genes=None,
    )
    
    if not success:
        return False
    
    # 验证输出
    print_header("第四步: 验证输出文件")
    print_info("验证输出文件...\n")
    
    try:
        with open(output_file, 'r') as f:
            lines = f.readlines()
            num_species = len([l for l in lines if l.startswith('>')])
            total_aa = sum(len(l) for l in lines if not l.startswith('>') and l.strip())
            
        print_success(f"输出文件中的物种数: {num_species}")
        print_success(f"氨基酸总数: {total_aa}")
        print_success(f"输出文件大小: {os.path.getsize(output_file) / 1024:.2f} KB")
        
    except Exception as e:
        print_error(f"验证文件失败: {e}")
        return False
    
    print_header("完成")
    print_success("脚本执行成功！\n")

    return True

if __name__ == "__main__":
    try:
        success = main()
        if success:
            sys.exit(0)
        else:
            sys.exit(1)
    except KeyboardInterrupt:
        print_error("\n用户中断执行")
        cleanup_on_exit()
        sys.exit(1)
    except Exception as e:
        print_error(f"发生未预期的错误: {e}")
        cleanup_on_exit()
        sys.exit(1)
