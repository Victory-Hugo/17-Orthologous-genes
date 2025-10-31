#!/usr/bin/env python3
# ============================================================
# Author: BigLin
# Purpose: 去除 FASTA 序列中的 "-" 并保持序列名称不变；输出 line wrap = 80
# ============================================================
#
# 运行环境与依赖：
#   Python >= 3.8
#   需要安装的包：
#       biopython
#       tqdm
#       termcolor
#
# 安装方式：
#   pip install biopython tqdm termcolor
#
# ============================================================

from Bio import SeqIO
from termcolor import colored
from tqdm import tqdm
import textwrap
import os
import sys

# ------------------ 固定参数设置 ------------------
INPUT = "/data_ssd3/7-luolintao-ssd/0-GTDB-Database/GTDB_bac_MSA_aln.faa"  # 输入文件
OUTPUT = "/data_ssd3/7-luolintao-ssd/0-GTDB-Database/GTDB_bac.faa"         # 输出文件
LINE_WRAP = 80                                                              # 每行最大长度
# -------------------------------------------------

# ---------- 彩色打印 ----------
def print_info(msg):
    print(colored(f"[INFO] {msg}", "cyan"))

def print_ok(msg):
    print(colored(f"[OK]   {msg}", "green"))

def print_warn(msg):
    print(colored(f"[WARN] {msg}", "yellow"))

def print_fail(msg):
    print(colored(f"[FAIL] {msg}", "red"))

# ---------- 主功能函数 ----------
def remove_gaps_stream(input_path, output_path, line_wrap=80):
    """流式读取 FASTA，去除 '-' 并逐条写入输出文件"""
    if not os.path.exists(input_path):
        print_fail(f"找不到输入文件: {input_path}")
        sys.exit(1)

    # 尝试估计序列数（仅用于进度条显示）
    try:
        total = sum(1 for line in open(input_path, "r") if line.startswith(">"))
    except Exception:
        total = None

    if total and total == 0:
        print_fail("未检测到任何序列，请检查输入文件格式。")
        sys.exit(1)

    print_info(f"输入文件: {input_path}")
    print_info(f"输出文件: {output_path}")
    print_info(f"检测到约 {total if total else '?'} 条序列，开始流式处理...")
    print_info("===========================================")

    processed = 0
    with open(output_path, "w") as out_f:
        for record in tqdm(SeqIO.parse(input_path, "fasta"),
                           total=total,
                           desc=colored("Processing", "magenta"),
                           ncols=80):
            seq_str = str(record.seq).replace("-", "")
            wrapped_seq = "\n".join(textwrap.wrap(seq_str, line_wrap))
            out_f.write(f">{record.id}\n{wrapped_seq}\n")
            processed += 1

    print_ok(f"所有序列已处理完毕，共 {processed} 条。")
    print_info(f"输出文件位置: {output_path}")

# ---------- 主入口 ----------
if __name__ == "__main__":
    print_info("========== 去除 FASTA 序列中的 '-' ==========")
    print_info(f"输入文件: {INPUT}")
    print_info(f"输出文件: {OUTPUT}")
    print_info(f"换行长度: {LINE_WRAP}")
    print_info("===========================================")

    remove_gaps_stream(INPUT, OUTPUT, LINE_WRAP)

    print_ok("任务完成 ✅")
