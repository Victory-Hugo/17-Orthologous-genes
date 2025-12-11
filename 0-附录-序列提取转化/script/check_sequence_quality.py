#!/usr/bin/env python3
"""
序列质量检查工具
检查 FASTA 文件中的编码序列是否符合标准

支持两种调用方式:
1. 作为模块 import 使用: from check_sequence_quality import check_sequences
2. 作为命令行脚本调用: python check_sequence_quality.py <输入> <输出>
"""

import sys
import json
import argparse
from pathlib import Path
from typing import Dict


def parse_fasta(fasta_file: str) -> Dict[str, str]:
    """解析 FASTA 文件
    
    Args:
        fasta_file: FASTA 文件路径
        
    Returns:
        序列 ID 到序列字符串的字典
    """
    sequences = {}
    current_id = None
    current_seq = ""
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = current_seq
                current_id = line.lstrip('>')
                current_seq = ""
            else:
                current_seq += line
        
        if current_id is not None:
            sequences[current_id] = current_seq
    
    return sequences


def analyze_sequence(seq_id: str, seq: str) -> Dict:
    """分析单个序列的质量
    
    Args:
        seq_id: 序列 ID
        seq: 序列字符串
        
    Returns:
        包含分析结果的字典
    """
    if not seq:
        return {
            'id': seq_id,
            'length': 0,
            'status': 'empty',
            'issues': ['空序列']
        }
    
    length = len(seq)
    start_codon = seq[:3] if len(seq) >= 3 else ""
    end_codon = seq[-3:] if len(seq) >= 3 else ""
    
    issues = []
    
    # 检查长度是否能被3整除
    if length % 3 != 0:
        issues.append(f"长度不能被3整除 (长度={length}, 余数={length%3})")
    
    # 检查起始密码子
    start_ok = start_codon == 'ATG'
    if not start_ok:
        issues.append(f"起始密码子不标准 (实际: {start_codon}, 期望: ATG)")
    
    # 检查终止密码子
    stop_ok = end_codon in ['TAA', 'TAG', 'TGA']
    if not stop_ok:
        issues.append(f"终止密码子错误 (实际: {end_codon}, 期望: TAA/TAG/TGA)")
    
    # 检查非标准核苷酸
    invalid = set(seq.upper()) - set('ATGC')
    if invalid:
        issues.append(f"含有非标准核苷酸: {','.join(invalid)}")
    
    # 统计密码子
    atg_count = seq.upper().count('ATG')
    stop_count = (seq.upper().count('TAA') + 
                  seq.upper().count('TAG') + 
                  seq.upper().count('TGA'))
    
    # 计算 GC 含量
    gc_count = seq.upper().count('G') + seq.upper().count('C')
    gc_content = (gc_count / length * 100) if length > 0 else 0
    
    # 判断序列状态
    if not issues:
        status = 'normal'
    elif len(issues) == 1 and 'ATG' in issues[0]:
        # 只有起始密码子不标准，可能是生物学真实
        status = 'alternative_start'
    else:
        status = 'abnormal'
    
    return {
        'id': seq_id,
        'length': length,
        'start_codon': start_codon,
        'end_codon': end_codon,
        'atg_count': atg_count,
        'stop_count': stop_count,
        'gc_content': round(gc_content, 2),
        'status': status,
        'issues': issues
    }


def _write_text_report(report: Dict, text_report_file: str) -> None:
    """生成易读的文本格式报告
    
    Args:
        report: 报告字典
        text_report_file: 输出文件路径
    """
    normal = report['normal_sequences']
    alternative = report['alternative_start_sequences']
    abnormal = report['abnormal_sequences']
    
    with open(text_report_file, 'w', encoding='utf-8') as f:
        f.write("=" * 100 + "\n")
        f.write("序列质量检查报告\n")
        f.write("=" * 100 + "\n\n")
        
        f.write("【统计摘要】\n")
        f.write(f"总序列数: {report['summary']['total_sequences']}\n")
        f.write(f"正常序列: {report['summary']['normal']} 个\n")
        f.write(f"非标准起始密码子: {report['summary']['alternative_start']} 个\n")
        f.write(f"异常序列: {report['summary']['abnormal']} 个\n\n")
        
        if abnormal:
            f.write("=" * 100 + "\n")
            f.write("【异常序列详情】\n")
            f.write("=" * 100 + "\n\n")
            for r in abnormal:
                f.write(f"ID: {r['id']}\n")
                f.write(f"  长度: {r['length']} bp (余数={r['length']%3})\n")
                f.write(f"  起始密码子: {r['start_codon']}\n")
                f.write(f"  末端密码子: {r['end_codon']}\n")
                f.write("  问题: \n")
                for issue in r['issues']:
                    f.write(f"    - {issue}\n")
                f.write("\n")
        
        if alternative:
            f.write("=" * 100 + "\n")
            f.write("【非标准起始密码子序列】(可能是生物学真实)\n")
            f.write("=" * 100 + "\n\n")
            for r in alternative:
                f.write(f"ID: {r['id']:<30} 起始: {r['start_codon']}  长度: {r['length']:>5} bp  GC: {r['gc_content']:>6.2f}%\n")
            f.write("\n")
        
        if normal:
            f.write("=" * 100 + "\n")
            f.write(f"【正常序列列表】({len(normal)} 个)\n")
            f.write("=" * 100 + "\n\n")
            for i, r in enumerate(normal, 1):
                f.write(f"{i:3d}. {r['id']:<30} 长度: {r['length']:>5} bp  GC: {r['gc_content']:>6.2f}%\n")


def check_sequences(
    input_file: str,
    output_file: str,
    verbose: bool = True
) -> Dict:
    """检查序列质量并生成报告
    
    这是主要的 run() 函数，支持被其他 Python 脚本 import 调用
    
    Args:
        input_file: 输入 FASTA 文件路径
        output_file: 输出报告文件路径 (JSON 格式，.json 扩展名)
        verbose: 是否输出详细信息到标准输出 (默认为 True)
        
    Returns:
        包含完整分析结果的报告字典
        
    Raises:
        FileNotFoundError: 如果输入文件不存在
        IOError: 如果无法读取或写入文件
    """
    # 验证输入文件
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"输入文件不存在: {input_file}")
    
    # 解析 FASTA
    if verbose:
        print(f"正在读取: {input_file}")
    sequences = parse_fasta(input_file)
    if verbose:
        print(f"找到 {len(sequences)} 个序列")
    
    # 分析所有序列
    if verbose:
        print("正在分析序列...")
    results = []
    for seq_id, seq in sequences.items():
        result = analyze_sequence(seq_id, seq)
        results.append(result)
    
    # 分类统计
    normal = [r for r in results if r['status'] == 'normal']
    alternative = [r for r in results if r['status'] == 'alternative_start']
    abnormal = [r for r in results if r['status'] == 'abnormal']
    
    # 生成报告字典
    report = {
        'summary': {
            'total_sequences': len(results),
            'normal': len(normal),
            'alternative_start': len(alternative),
            'abnormal': len(abnormal)
        },
        'normal_sequences': normal,
        'alternative_start_sequences': alternative,
        'abnormal_sequences': abnormal
    }
    
    # 保存 JSON 报告
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(report, f, ensure_ascii=False, indent=2)
    
    # 同时生成易读的文本报告
    text_report_file = output_file.replace('.json', '.txt')
    _write_text_report(report, text_report_file)
    
    # 输出统计信息
    if verbose:
        print("\n【统计摘要】")
        print(f"总序列数: {report['summary']['total_sequences']}")
        print(f"正常序列: {report['summary']['normal']} 个")
        print(f"非标准起始密码子: {report['summary']['alternative_start']} 个")
        print(f"异常序列: {report['summary']['abnormal']} 个")
        
        if report['summary']['abnormal'] > 0:
            print(f"\n⚠️  发现 {report['summary']['abnormal']} 个异常序列:")
            for r in report['abnormal_sequences']:
                print(f"  - {r['id']}: {', '.join(r['issues'])}")
        
        print(f"\n✓ 报告已保存到: {output_file}")
        print(f"✓ 文本报告已保存到: {text_report_file}")
    
    return report


def main():
    """命令行入口点"""
    parser = argparse.ArgumentParser(
        description='序列质量检查工具 - 检查 FASTA 文件中的编码序列是否符合标准',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例:
  python check_sequence_quality.py input.fna output/report.json
  python check_sequence_quality.py /path/to/sequences.fna ./reports/report.json
        '''
    )
    
    parser.add_argument(
        'input_file',
        help='输入 FASTA 文件路径'
    )
    
    parser.add_argument(
        'output_file',
        help='输出报告文件路径 (建议使用 .json 扩展名)'
    )
    
    parser.add_argument(
        '-q', '--quiet',
        action='store_true',
        help='静默模式，不输出详细信息'
    )
    
    args = parser.parse_args()
    
    try:
        check_sequences(
            input_file=args.input_file,
            output_file=args.output_file,
            verbose=not args.quiet
        )
    except FileNotFoundError as e:
        print(f"错误: {e}", file=sys.stderr)
        sys.exit(1)
    except IOError as e:
        print(f"I/O 错误: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
