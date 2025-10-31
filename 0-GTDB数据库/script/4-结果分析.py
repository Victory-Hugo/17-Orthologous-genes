#!/usr/bin/env python3

"""
GTDB-Tk 结果分析脚本
用法: python3 4-结果分析.py <输出目录>
"""

import os
import sys
import pandas as pd
from pathlib import Path

def analyze_gtdbtk_results(output_dir):
    """分析 GTDB-Tk 分类结果"""
    
    output_path = Path(output_dir)
    
    if not output_path.exists():
        print(f"❌ 错误: 输出目录不存在: {output_dir}")
        sys.exit(1)
    
    print("=" * 60)
    print("GTDB-Tk 结果分析")
    print("=" * 60)
    print(f"结果目录: {output_dir}")
    print()
    
    # 查找摘要文件
    classify_summary = output_path / "classify" / "classify_wf.summary.tsv"
    markers_summary = output_path / "classify" / "classify_wf.markers_summary.tsv"
    
    if not classify_summary.exists():
        print("❌ 找不到分类摘要文件")
        return
    
    print("📊 分类结果摘要")
    print("-" * 60)
    
    try:
        df_classify = pd.read_csv(classify_summary, sep='\t')
        print(f"基因组总数: {len(df_classify)}")
        print()
        
        # 分析分类结果
        print("分类类型分布:")
        classification_counts = df_classify['classification'].value_counts()
        for clas_type, count in classification_counts.items():
            print(f"  {clas_type}: {count} ({count/len(df_classify)*100:.1f}%)")
        
        # 检查质量指标
        if 'fastq_error_rate' in df_classify.columns:
            print()
            print("质量统计:")
            print(f"  平均错误率: {df_classify['fastq_error_rate'].mean():.4f}")
        
        # 显示前几行详细信息
        print()
        print("详细分类结果 (前10个基因组):")
        print("-" * 60)
        display_cols = ['user_genome', 'classification', 'fastq_error_rate'] 
        available_cols = [col for col in display_cols if col in df_classify.columns]
        print(df_classify[available_cols].head(10).to_string(index=False))
        
    except Exception as e:
        print(f"❌ 读取分类文件出错: {e}")
    
    # 分析标记基因
    if markers_summary.exists():
        print()
        print("📊 标记基因摘要")
        print("-" * 60)
        
        try:
            df_markers = pd.read_csv(markers_summary, sep='\t')
            print(f"基因组总数: {len(df_markers)}")
            print()
            
            # 分析完整性
            if 'completeness' in df_markers.columns:
                print("完整性统计 (Completeness):")
                print(f"  平均值: {df_markers['completeness'].mean():.2f}%")
                print(f"  最小值: {df_markers['completeness'].min():.2f}%")
                print(f"  最大值: {df_markers['completeness'].max():.2f}%")
            
            # 分析污染度
            if 'contamination' in df_markers.columns:
                print()
                print("污染度统计 (Contamination):")
                print(f"  平均值: {df_markers['contamination'].mean():.2f}%")
                print(f"  最小值: {df_markers['contamination'].min():.2f}%")
                print(f"  最大值: {df_markers['contamination'].max():.2f}%")
            
            # 质量评级
            if 'genome_quality' in df_markers.columns:
                print()
                print("基因组质量等级:")
                quality_counts = df_markers['genome_quality'].value_counts()
                for quality, count in quality_counts.items():
                    print(f"  {quality}: {count}")
            
            # 显示详细信息
            print()
            print("标记基因详情 (前10个基因组):")
            print("-" * 60)
            display_cols = ['genome', 'completeness', 'contamination', 'genome_quality']
            available_cols = [col for col in display_cols if col in df_markers.columns]
            print(df_markers[available_cols].head(10).to_string(index=False))
            
        except Exception as e:
            print(f"❌ 读取标记基因文件出错: {e}")
    
    # 统计输出文件
    print()
    print("📂 输出文件统计")
    print("-" * 60)
    
    output_files = list(output_path.glob("**/*"))
    file_types = {}
    for f in output_files:
        if f.is_file():
            ext = f.suffix or "no_ext"
            file_types[ext] = file_types.get(ext, 0) + 1
    
    for ext, count in sorted(file_types.items(), key=lambda x: x[1], reverse=True):
        print(f"  {ext}: {count} 个文件")
    
    print()
    print("=" * 60)
    print("✅ 分析完成")
    print("=" * 60)


def main():
    if len(sys.argv) < 2:
        print("用法: python3 4-结果分析.py <输出目录>")
        print("示例: python3 4-结果分析.py ./gtdbtk_output")
        sys.exit(1)
    
    output_dir = sys.argv[1]
    analyze_gtdbtk_results(output_dir)


if __name__ == "__main__":
    main()
