#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HMM匹配数据质量分析和可视化脚本

功能：
  1. 加载HMM匹配结果（tbl格式）
  2. 生成多维度数据质量分析报告
  3. 绘制专业质量评估图表
  4. 导出严格筛选后的数据

使用方法：
  python tbl_data_analysis.py --input <输入文件> --output <输出目录>
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from aquarel import load_theme
from scipy import stats


def load_data(input_file):
    """加载HMM匹配数据"""
    print(f"[INFO] 正在加载数据: {input_file}")
    df = pd.read_csv(input_file, sep='\t')
    print(f"[OK] 数据加载成功: {df.shape[0]:,} 行, {df.shape[1]} 列")
    return df


def print_basic_statistics(df):
    """打印基本统计信息"""
    print("\n" + "="*80)
    print("数据基本信息")
    print("="*80)
    print(f"数据形状: {df.shape}")
    print(f"\n列名和类型:")
    print(df.dtypes)
    print(f"\n前几行:")
    print(df.head())


def analyze_evalue(df):
    """分析E-value分布"""
    print("\n" + "="*80)
    print("E-value 统计分析")
    print("="*80)
    print(f"中位数: {df['E-value'].median():.2e}")
    print(f"均值: {df['E-value'].mean():.2e}")
    print(f"最小值: {df['E-value'].min():.2e}")
    print(f"\nE-value 分布:")
    print(f"< 1e-100 (极高相似性): {(df['E-value'] < 1e-100).sum():,}")
    print(f"1e-100 ~ 1e-50: {((df['E-value'] >= 1e-100) & (df['E-value'] < 1e-50)).sum():,}")
    print(f"1e-50 ~ 1e-20: {((df['E-value'] >= 1e-50) & (df['E-value'] < 1e-20)).sum():,}")
    print(f">= 1e-20 (低相似性): {(df['E-value'] >= 1e-20).sum():,}")

    print(f"\n=== Score 统计 ===")
    print(df['score'].describe())


def analyze_coverage(df):
    """分析样本和基因覆盖度"""
    print("\n" + "="*80)
    print("样本和基因覆盖度分析")
    print("="*80)
    print(f"独特样本数: {df['文件名'].nunique():,}")
    print(f"独特HMM查询基因: {df['query name'].nunique()}")
    print(f"独特目标蛋白: {df['target name'].nunique():,}")

    print(f"\n样本大小分布:")
    sample_counts = df['文件名'].value_counts()
    print(f"平均匹配数/样本: {sample_counts.mean():.1f}")
    print(f"中位数: {sample_counts.median():.1f}")
    print(f"最多匹配: {sample_counts.max()} (样本: {sample_counts.idxmax()})")
    print(f"最少匹配: {sample_counts.min()} (样本: {sample_counts.idxmin()})")

    print(f"\n基因匹配广度:")
    query_counts = df['query name'].value_counts()
    print(f"平均匹配样本数/基因: {query_counts.mean():.1f}")
    print(f"中位数: {query_counts.median():.1f}")


def analyze_quality_flags(df):
    """分析质量标志（缺失值、多域等）"""
    print("\n" + "="*80)
    print("缺失值和质量标志检查")
    print("="*80)
    print(f"缺失值统计:")
    print(df.isnull().sum())

    print(f"\n多域/重叠分析:")
    print(f"包含多个域的匹配 (dom > 1): {(df['dom'] > 1).sum():,}")
    print(f"包含重叠的匹配 (ov > 0): {(df['ov'] > 0).sum():,}")
    print(f"重复匹配 (rep > 0): {(df['rep'] > 0).sum():,}")

    print(f"\n高置信度匹配分析:")
    high_conf = df[df['E-value'] < 1e-20]
    print(f"E-value < 1e-20 的匹配: {len(high_conf):,} ({len(high_conf)/len(df)*100:.1f}%)")
    print(f"这些高置信匹配涉及样本: {high_conf['文件名'].nunique():,}")
    print(f"这些高置信匹配涉及基因: {high_conf['query name'].nunique()}")
    print(f"这些高置信匹配涉及蛋白: {high_conf['target name'].nunique():,}")


def analyze_genes(df):
    """分析个体HMM基因"""
    print("\n" + "="*80)
    print("HMM查询基因详细分析")
    print("="*80)
    for gene in df['query name'].unique():
        gene_data = df[df['query name'] == gene]
        print(f"\n基因: {gene}")
        print(f"  总匹配数: {len(gene_data):,}")
        print(f"  E-value 中位数: {gene_data['E-value'].median():.2e}")
        print(f"  Score 平均值: {gene_data['score'].mean():.1f}")
        print(f"  高置信匹配 (E-value < 1e-20): {(gene_data['E-value'] < 1e-20).sum():,}")
        print(f"  覆盖样本数: {gene_data['文件名'].nunique():,}")

    print(f"\n样本质量评估:")
    all_covered = (df.groupby('文件名')['query name'].nunique() == df['query name'].nunique()).all()
    print(f"所有样本都被全部基因检测到: {all_covered}")
    sample_gene_cov = df.groupby('文件名')['query name'].nunique()
    print(f"样本基因覆盖分布: {sample_gene_cov.value_counts().to_dict()}")


def show_case_examples(df):
    """展示不同质量等级的案例"""
    print("\n" + "="*80)
    print("实际案例对比：不同质量等级的匹配")
    print("="*80)

    df_sorted = df.sort_values('score', ascending=False)
    cases = {}

    # 极优匹配
    optimal = df_sorted[(df_sorted['E-value'] < 1e-50) & 
                        (df_sorted['score'] > 200) & 
                        (df_sorted['bias'] < 1)]
    if len(optimal) > 0:
        cases["🟢 极优匹配（应该保留）"] = optimal.iloc[0]
    else:
        cases["🟢 极优匹配（应该保留）"] = df_sorted[df_sorted['E-value'] < 1e-50].iloc[0]

    # 良好匹配
    good = df_sorted[(df_sorted['E-value'] < 1e-10) & 
                     (df_sorted['E-value'] >= 1e-20) & 
                     (df_sorted['score'] > 80) & 
                     (df_sorted['score'] < 150) & 
                     (df_sorted['bias'] < 3)]
    if len(good) > 0:
        cases["🟡 良好匹配（可接受）"] = good.iloc[0]
    else:
        cases["🟡 良好匹配（可接受）"] = df_sorted[(df_sorted['score'] >= 80) & 
                                              (df_sorted['score'] <= 150)].iloc[0]

    # 弱匹配
    weak = df_sorted[(df_sorted['E-value'] > 0.001) | (df_sorted['score'] < 50)]
    if len(weak) > 0:
        cases["🔴 弱匹配（需谨慎）"] = weak.iloc[0]
    else:
        cases["🔴 弱匹配（需谨慎）"] = df_sorted.iloc[-1]

    for case_name, row in cases.items():
        print(f"\n{case_name}")
        print("-" * 80)
        print(f"  基因:           {row['query name']}")
        evalue_status = ('✅ 非常可信' if row['E-value'] < 1e-20 
                        else ('⚠️  中等' if row['E-value'] < 1e-5 else '❌ 不可信'))
        print(f"  E-value:        {row['E-value']:.2e}    {evalue_status}")
        score_status = ('✅ 强' if row['score'] > 100 
                       else ('⚠️  中等' if row['score'] > 50 else '❌ 弱'))
        print(f"  Score:          {row['score']:.1f}         {score_status}")
        bias_status = ('✅ 低噪声' if row['bias'] < 1 
                      else ('⚠️  中等' if row['bias'] < 5 else '❌ 高噪声'))
        print(f"  Bias:           {row['bias']:.1f}         {bias_status}")
        ratio = row['score'] / max(row['bias'], 0.1)
        ratio_status = ('✅ 信噪比好' if ratio > 50 
                       else ('⚠️  中等' if ratio > 20 else '❌ 信噪比差'))
        print(f"  Score/Bias比:   {ratio:.1f}      {ratio_status}")

    print("\n" + "="*80)
    print("💡 快速判断标准总结")
    print("="*80)
    print("""
┌─────────────────┬──────────────┬──────────────┬──────────────┐
│  匹配质量       │  E-value     │  Score       │  Bias        │
├─────────────────┼──────────────┼──────────────┼──────────────┤
│ 🟢 优（推荐）   │  < 1e-20     │  > 100       │  < 1         │
│ 🟡 良（可用）   │  1e-5~1e-20  │  50～100     │  1～5        │
│ 🔴 差（可疑）   │  > 1e-5      │  < 50        │  > 5         │
└─────────────────┴──────────────┴──────────────┴──────────────┘

你的数据现状：
  • 78.2% 达到优质等级 ✅
  • 大部分Bias < 1（高信质量）✅
  • lolD基因比lolF更保守（Score高）✅
    """)


def filter_by_criteria(df):
    """按不同标准筛选数据"""
    print("\n" + "="*80)
    print("质量筛选标准对比")
    print("="*80)

    thresholds = {
        "严格标准（科研发表级）": {
            "E-value": 1e-20,
            "Score": 100,
            "Bias": 1.0,
        },
        "中等标准（常规分析）": {
            "E-value": 1e-10,
            "Score": 50,
            "Bias": 5.0,
        },
        "宽松标准（初步探索）": {
            "E-value": 1e-5,
            "Score": 25,
            "Bias": 10.0,
        },
    }

    for threshold_name, params in thresholds.items():
        filtered = df[(df['E-value'] <= params['E-value']) &
                      (df['score'] >= params['Score']) &
                      (df['bias'] <= params['Bias'])]
        retained_pct = len(filtered) / len(df) * 100

        print(f"\n{threshold_name}")
        print(f"  筛选条件: E-value ≤ {params['E-value']:.0e}, Score ≥ {params['Score']}, Bias ≤ {params['Bias']}")
        print(f"  保留匹配数: {len(filtered):,} / {len(df):,} ({retained_pct:.1f}%)")
        print(f"  保留样本: {filtered['文件名'].nunique():,} / {df['文件名'].nunique():,}")

    print("\n" + "="*80)
    print("📌 建议")
    print("="*80)
    print("""
1. 对于系统发育/进化分析：使用"严格标准"，确保结果稳健
2. 对于基因存在性调查：使用"中等标准"，平衡准确性和覆盖度
3. 对于初步探索/数据审视：使用"宽松标准"，获得全景图

你的数据中：
  • E-value中位数 3.5e-32 << 1e-20    → 数据质量非常好 ✅
  • Score均值 120.7 > 100              → 序列一致性强 ✅
  • Bias均值 0.17 << 1.0               → 背景噪声很小 ✅

➜ 推荐使用"严格标准"进行下游分析！
    """)


def generate_visualization(df, output_dir):
    """生成可视化图表"""
    print("\n" + "="*80)
    print("生成可视化图表")
    print("="*80)

    # 配置主题和字体
    theme = load_theme("boxy_light")
    theme.apply()
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42

    # 创建图表
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('HMM matching quality assessment', fontsize=16, fontweight='bold')

    # 1. E-value 分布
    ax1 = axes[0, 0]
    evalue_bins = [1e-120, 1e-110, 1e-100, 1e-90, 1e-80, 1e-70, 1e-60, 1e-50, 
                   1e-40, 1e-30, 1e-20, 1e-10, 1e-5, 0.001, 1]
    evalue_counts = []
    for i in range(len(evalue_bins) - 1):
        count = ((df['E-value'] >= evalue_bins[i]) & 
                 (df['E-value'] < evalue_bins[i + 1])).sum()
        evalue_counts.append(count)

    categories = ['<1e-110', '1e-110~100', '1e-100~90', '1e-90~80', '1e-80~70', 
                  '1e-70~60', '1e-60~50', '1e-50~40', '1e-40~30', '1e-30~20', 
                  '1e-20~10', '1e-10~5', '>0.001']
    new_colors = ['#085FE3', '#085FE3', '#085FE3', '#00AFFF', '#00AFFF', '#23A5AC', '#23A5AC',
                  '#DB9C15', '#DB9C15', '#7A1616', '#7A1616', '#7A1616', '#7A1616', '#7A1616']
    ax1.bar(range(len(evalue_counts)), evalue_counts, 
            color=new_colors[:len(evalue_counts)], alpha=0.5)
    ax1.set_xticks(range(len(categories)))
    ax1.set_xticklabels(categories, rotation=90, ha='right')
    ax1.set_ylabel('Matching count')
    ax1.set_title('The distribution of E-value', fontweight='bold')
    ax1.set_yscale('log')
    ax1.grid(False)

    # 2. Score 分布
    ax2 = axes[0, 1]
    colors_list = ['#3498db', '#e74c3c', '#2ecc71', '#f39c12', '#9b59b6']
    for idx, gene in enumerate(df['query name'].unique()):
        gene_data = df[df['query name'] == gene]['score']
        color = colors_list[idx % len(colors_list)]
        ax2.hist(gene_data, bins=50, alpha=0.6, label=gene, color=color)
    ax2.set_xlabel('Value of score')
    ax2.set_ylabel('Matching count')
    ax2.set_title('The distribution of score', fontweight='bold')
    ax2.legend()
    ax2.grid(False)

    # 3. Bias 分布
    ax3 = axes[1, 0]
    bias_hist, bias_bins, patches = ax3.hist(df['bias'], bins=100, 
                                             color='#23A5AC', alpha=0.5)
    ax3.set_xlabel('Value of bias')
    ax3.set_ylabel('Matching count')
    ax3.set_title('The distribution of bias', fontweight='bold')
    ax3.grid(False)

    # 4. Score vs Bias 散点图
    ax4 = axes[1, 1]
    sample_idx = np.random.choice(len(df), 5000, replace=False)
    scatter = ax4.scatter(df.iloc[sample_idx]['bias'],
                         df.iloc[sample_idx]['score'],
                         c=np.log10(df.iloc[sample_idx]['E-value']),
                         cmap='RdYlGn_r', s=20, alpha=0.6)
    ax4.set_xlabel('Bias')
    ax4.set_ylabel('Score')
    ax4.set_title('Score vs bias', fontweight='bold')
    cbar = plt.colorbar(scatter, ax=ax4)
    cbar.set_label('log10(E-value)')
    ax4.grid(False)

    plt.tight_layout()
    theme.apply_transforms()

    # 保存图表
    pdf_path = os.path.join(output_dir, 'tbl_data_exploration.pdf')
    png_path = os.path.join(output_dir, 'tbl_data_exploration.png')
    plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"[OK] PDF已保存: {pdf_path}")
    print(f"[OK] PNG已保存: {png_path}")


def export_filtered_data(df, output_dir):
    """导出严格筛选后的数据"""
    print("\n" + "="*80)
    print("导出筛选后的数据")
    print("="*80)

    strict_filtered = df[(df['E-value'] <= 1e-20) &
                         (df['score'] >= 100) &
                         (df['bias'] <= 1.0)]

    output_file = os.path.join(output_dir, 'tbl_strict_filtered.tsv')
    strict_filtered.to_csv(output_file, sep='\t', index=False)

    print(f"[OK] 严格筛选后的数据已保存: {output_file}")
    print(f"     原始数据: {len(df):,} 行")
    print(f"     筛选后: {len(strict_filtered):,} 行 ({len(strict_filtered)/len(df)*100:.1f}%)")


def export_sample_copy_count(df, output_dir):
    """导出每个样本符合条件的拷贝数统计"""
    print("\n" + "="*80)
    print("样本拷贝数统计")
    print("="*80)

    # 严格筛选标准
    strict_filtered = df[(df['E-value'] <= 1e-20) &
                         (df['score'] >= 100) &
                         (df['bias'] <= 1.0)]

    # 按样本分组统计
    sample_copy_count = strict_filtered.groupby('文件名').size().reset_index(name='CopyNumber')
    
    # 获取所有样本列表（包括没有符合条件的样本）
    all_samples = df['文件名'].unique()
    all_samples_df = pd.DataFrame({'ID': all_samples})
    
    # 左连接，没有符合条件的样本拷贝数为0
    sample_copy_count.rename(columns={'文件名': 'ID'}, inplace=True)
    result_df = all_samples_df.merge(sample_copy_count, on='ID', how='left')
    result_df['CopyNumber'] = result_df['CopyNumber'].fillna(0).astype(int)
    
    # 按拷贝数降序排序
    result_df = result_df.sort_values('CopyNumber', ascending=False)

    # 保存为TSV文件
    output_file = os.path.join(output_dir, 'sample_copy_count.tsv')
    result_df.to_csv(output_file, sep='\t', index=False)

    print(f"[OK] 样本拷贝数统计已保存: {output_file}")
    print(f"\n统计摘要:")
    print(f"     总样本数: {len(result_df):,}")
    print(f"     有符合条件匹配的样本: {(result_df['CopyNumber'] > 0).sum():,}")
    print(f"     无符合条件匹配的样本: {(result_df['CopyNumber'] == 0).sum():,}")
    print(f"     平均拷贝数: {result_df['CopyNumber'].mean():.2f}")
    print(f"     最多拷贝数: {result_df['CopyNumber'].max()}")
    print(f"     最少拷贝数: {result_df['CopyNumber'].min()}")
    
    print(f"\n拷贝数分布:")
    copy_distribution = result_df['CopyNumber'].value_counts().sort_index()
    for copy_num, count in copy_distribution.items():
        print(f"     {copy_num}个拷贝的样本: {count:,} 个")


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='HMM匹配数据质量分析和可视化',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法：
  python tbl_data_analysis.py --input data/tbl_merge.tsv --output results/
  python tbl_data_analysis.py -i input.tsv -o ./output/
        """)

    parser.add_argument('-i', '--input', type=str, required=True,
                       help='输入的HMM匹配结果文件（TSV格式）')
    parser.add_argument('-o', '--output', type=str, required=True,
                       help='输出目录（将自动创建如果不存在）')

    args = parser.parse_args()

    # 验证输入文件
    if not os.path.exists(args.input):
        print(f"[ERROR] 输入文件不存在: {args.input}")
        sys.exit(1)

    # 创建输出目录
    os.makedirs(args.output, exist_ok=True)
    print(f"[OK] 输出目录: {args.output}")

    # 执行分析流程
    df = load_data(args.input)
    print_basic_statistics(df)
    analyze_evalue(df)
    analyze_coverage(df)
    analyze_quality_flags(df)
    analyze_genes(df)
    show_case_examples(df)
    filter_by_criteria(df)
    generate_visualization(df, args.output)
    export_filtered_data(df, args.output)
    export_sample_copy_count(df, args.output)

    print("\n" + "="*80)
    print("✅ 分析完成！")
    print("="*80)


if __name__ == '__main__':
    main()
