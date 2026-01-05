#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Rv1819c 基因分类统计分析脚本

功能：
1. 读取NCBI物种分类信息和过滤后的结果
2. 数据合并和预处理
3. 生成多种可视化图表：
   - 热图：展示各分类级别的丰度分布
   - 堆叠柱状图：展示分类级别丰度分布
   - 总体统计图表
   - 分类比例对比图
4. 生成详细的统计数据表格

作者: 
创建时间: 2026-01-05
修改时间: 2026-01-05
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

# ==================== 配置参数 ====================
META_FILE = "/mnt/d/5-NCBI-Reference/1-Bac/meta/final_meta.csv"  # NCBI 物种分类信息文件
FILTER_FILE = "/mnt/l/19-Rv1819c-Gene/1-BLAST/output/Final/rv1819c_filter_all.tsv"  # 过滤后的最终结果
OUT_PATH = "/mnt/l/19-Rv1819c-Gene/1-BLAST/output/Final"  # 输出文件夹

# ==================== 数据读取和预处理 ====================
print("开始读取数据文件...")

# 读取数据文件
df_META = pd.read_csv(META_FILE)
df_FILTER = pd.read_csv(FILTER_FILE, sep="\t")

print("数据读取完成!")
print(f"META文件包含 {len(df_META)} 条记录")
print(f"FILTER文件包含 {len(df_FILTER)} 条记录")

# 数据预处理
print("\n开始数据预处理...")

# 提取Assembly信息
df_FILTER['Assembly'] = df_FILTER['seq_id'].str.split('|').str[0]
df_FILTER.loc[df_FILTER['Assembly']=='Rv1819c', 'Assembly'] = 'GCF_000195955.2'  # 将 Rv1819c 替换为其对应的 Assembly 编号

# 数据合并
df_merge = df_FILTER.loc[:,['Assembly','class']].merge(df_META, on='Assembly', how='outer')

# 将class==NaN的命名为C
df_merge['class'] = df_merge['class'].fillna('C')

print(f"合并后数据包含 {len(df_merge)} 条记录")

# 保存合并后的数据
df_merge.to_csv(OUT_PATH + "/rv1819c_filter_meta.csv", sep=",", index=False)
print(f"合并数据已保存到: {OUT_PATH}/rv1819c_filter_meta.csv")

# ==================== 可视化准备 ====================
print("\n开始准备可视化...")

# 确保输出目录存在
os.makedirs(OUT_PATH, exist_ok=True)

# 设置图形样式
plt.style.use('default')
sns.set_palette("husl")

# 定义分类级别
taxonomic_levels = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

# 为每个分类级别创建丰度统计表
print("\n生成各分类级别的丰度统计表...")
abundance_data = {}
for level in taxonomic_levels:
    # 计算每个分类单元在每个类别中的丰度
    level_counts = df_merge.groupby([level, 'class']).size().unstack(fill_value=0)
    abundance_data[level] = level_counts
    
    # 保存数据表格
    level_counts.to_csv(f"{OUT_PATH}/abundance_{level.lower()}.csv")
    print(f"已保存 {level} 级别丰度数据到: {OUT_PATH}/abundance_{level.lower()}.csv")

# ==================== 可视化函数定义 ====================

def create_heatmap(data_dict, output_path):
    """创建热图展示分类丰度分布"""
    
    print("\n正在创建热图...")
    
    # 计算需要的子图数量
    n_levels = len(data_dict)
    fig, axes = plt.subplots(n_levels, 1, figsize=(12, 4*n_levels))
    
    if n_levels == 1:
        axes = [axes]
    
    for i, (level, data) in enumerate(data_dict.items()):
        # 只显示前20个最丰富的分类单元（如果超过20个）
        if len(data) > 20:
            # 按总丰度排序，取前20个
            data_sorted = data.loc[data.sum(axis=1).nlargest(20).index]
        else:
            data_sorted = data
        
        # 创建热图
        sns.heatmap(data_sorted, 
                   annot=True, 
                   fmt='d', 
                   cmap='Blues',
                   ax=axes[i],
                   cbar_kws={'shrink': 0.8})
        
        axes[i].set_title(f'{level} 级别丰度分布热图', fontsize=14, fontweight='bold')
        axes[i].set_xlabel('分类 (Class)', fontsize=12)
        axes[i].set_ylabel(f'{level}', fontsize=12)
        
        # 旋转y轴标签以便阅读
        axes[i].tick_params(axis='y', rotation=0)
        axes[i].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}/taxonomic_abundance_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()  # 关闭图形，节省内存
    print(f"热图已保存到: {output_path}/taxonomic_abundance_heatmap.png")


def create_stacked_barplot(data_dict, output_path):
    """创建堆叠柱状图展示分类丰度分布"""
    
    print("\n正在创建堆叠柱状图...")
    
    # 计算需要的子图数量
    n_levels = len(data_dict)
    fig, axes = plt.subplots(n_levels, 1, figsize=(15, 5*n_levels))
    
    if n_levels == 1:
        axes = [axes]
    
    # 定义颜色
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    
    for i, (level, data) in enumerate(data_dict.items()):
        # 只显示前15个最丰富的分类单元（如果超过15个）
        if len(data) > 15:
            data_sorted = data.loc[data.sum(axis=1).nlargest(15).index]
        else:
            data_sorted = data
        
        # 创建堆叠柱状图
        data_sorted.plot(kind='bar', 
                        stacked=True, 
                        ax=axes[i],
                        color=colors,
                        figsize=(15, 5))
        
        axes[i].set_title(f'{level} 级别丰度分布堆叠柱状图', fontsize=14, fontweight='bold')
        axes[i].set_xlabel(f'{level}', fontsize=12)
        axes[i].set_ylabel('物种数量', fontsize=12)
        axes[i].legend(title='分类 (Class)', loc='upper right')
        
        # 旋转x轴标签以便阅读
        axes[i].tick_params(axis='x', rotation=45, labelsize=10)
        axes[i].grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}/taxonomic_abundance_stacked_barplot.png", dpi=300, bbox_inches='tight')
    plt.close()  # 关闭图形，节省内存
    print(f"堆叠柱状图已保存到: {output_path}/taxonomic_abundance_stacked_barplot.png")


def create_overall_statistics(df, output_path):
    """创建总体统计图表"""
    
    print("\n正在创建总体统计图表...")
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # 3.1 各类别总体分布饼图
    class_counts = df['class'].value_counts()
    axes[0, 0].pie(class_counts.values, 
                   labels=class_counts.index, 
                   autopct='%1.1f%%',
                   colors=['#ff9999', '#66b3ff', '#99ff99'])
    axes[0, 0].set_title('各类别 (Class) 总体分布', fontsize=14, fontweight='bold')
    
    # 3.2 各门 (Phylum) 分布条形图
    phylum_counts = df['Phylum'].value_counts().head(10)
    sns.barplot(x=phylum_counts.values, 
                y=phylum_counts.index, 
                ax=axes[0, 1],
                palette='viridis')
    axes[0, 1].set_title('前10个门 (Phylum) 分布', fontsize=14, fontweight='bold')
    axes[0, 1].set_xlabel('物种数量')
    
    # 3.3 各纲 (Class) 分布条形图
    class_tax_counts = df['Class'].value_counts().head(10)
    sns.barplot(x=class_tax_counts.values, 
                y=class_tax_counts.index, 
                ax=axes[1, 0],
                palette='plasma')
    axes[1, 0].set_title('前10个纲 (Class) 分布', fontsize=14, fontweight='bold')
    axes[1, 0].set_xlabel('物种数量')
    
    # 3.4 各科 (Family) 分布条形图
    family_counts = df['Family'].value_counts().head(10)
    sns.barplot(x=family_counts.values, 
                y=family_counts.index, 
                ax=axes[1, 1],
                palette='cividis')
    axes[1, 1].set_title('前10个科 (Family) 分布', fontsize=14, fontweight='bold')
    axes[1, 1].set_xlabel('物种数量')
    
    plt.tight_layout()
    plt.savefig(f"{output_path}/overall_statistics.png", dpi=300, bbox_inches='tight')
    plt.close()  # 关闭图形，节省内存
    print(f"总体统计图表已保存到: {output_path}/overall_statistics.png")


def create_summary_tables(df, abundance_data, output_path):
    """生成详细的汇总统计表"""
    
    print("\n正在生成汇总统计表...")
    
    # 4.1 总体汇总统计
    summary_stats = pd.DataFrame({
        '分类级别': taxonomic_levels,
        '总类群数': [df[level].nunique() for level in taxonomic_levels],
        'A类物种数': [len(abundance_data[level][abundance_data[level]['A'] > 0]) if 'A' in abundance_data[level].columns else 0 for level in taxonomic_levels],
        'B类物种数': [len(abundance_data[level][abundance_data[level]['B'] > 0]) if 'B' in abundance_data[level].columns else 0 for level in taxonomic_levels],
        'C类物种数': [len(abundance_data[level][abundance_data[level]['C'] > 0]) if 'C' in abundance_data[level].columns else 0 for level in taxonomic_levels]
    })
    
    # 4.2 各分类级别的详细统计
    detailed_stats = {}
    for level in taxonomic_levels:
        level_summary = []
        for taxa in df[level].value_counts().head(20).index:
            taxa_data = df[df[level] == taxa]['class'].value_counts()
            level_summary.append({
                f'{level}': taxa,
                'A类数量': taxa_data.get('A', 0),
                'B类数量': taxa_data.get('B', 0),
                'C类数量': taxa_data.get('C', 0),
                '总数': taxa_data.sum()
            })
        detailed_stats[level] = pd.DataFrame(level_summary)
    
    # 保存汇总统计表
    summary_stats.to_csv(f"{output_path}/summary_statistics.csv", index=False)
    print(f"总体汇总统计已保存到: {output_path}/summary_statistics.csv")
    
    # 保存各分类级别详细统计
    for level, df_detail in detailed_stats.items():
        df_detail.to_csv(f"{output_path}/detailed_statistics_{level.lower()}.csv", index=False)
        print(f"{level} 详细统计已保存到: {output_path}/detailed_statistics_{level.lower()}.csv")
    
    return summary_stats, detailed_stats


def create_proportion_comparison(df, output_path):
    """创建各分类级别中A、B、C类比例对比图"""
    
    print("\n正在创建分类比例对比图...")
    
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.flatten()
    
    for i, level in enumerate(taxonomic_levels):
        # 计算前10个最丰富分类单元的比例
        top_taxa = df[level].value_counts().head(10).index
        proportions = []
        labels = []
        
        for taxa in top_taxa:
            taxa_data = df[df[level] == taxa]['class'].value_counts()
            total = taxa_data.sum()
            prop_A = taxa_data.get('A', 0) / total * 100
            prop_B = taxa_data.get('B', 0) / total * 100
            prop_C = taxa_data.get('C', 0) / total * 100
            
            proportions.append([prop_A, prop_B, prop_C])
            # 截断过长的标签
            label = taxa if len(taxa) <= 20 else taxa[:17] + "..."
            labels.append(label)
        
        # 创建堆叠条形图
        proportions = np.array(proportions)
        x = np.arange(len(labels))
        width = 0.8
        
        axes[i].bar(x, proportions[:, 0], width, label='A类', color='#ff9999')
        axes[i].bar(x, proportions[:, 1], width, bottom=proportions[:, 0], label='B类', color='#66b3ff')
        axes[i].bar(x, proportions[:, 2], width, bottom=proportions[:, 0] + proportions[:, 1], label='C类', color='#99ff99')
        
        axes[i].set_title(f'{level} 级别各类比例', fontsize=12, fontweight='bold')
        axes[i].set_ylabel('百分比 (%)')
        axes[i].set_xticks(x)
        axes[i].set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
        axes[i].legend(loc='upper right')
        axes[i].grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}/proportion_comparison.png", dpi=300, bbox_inches='tight')
    plt.close()  # 关闭图形，节省内存
    print(f"分类比例对比图已保存到: {output_path}/proportion_comparison.png")


# ==================== 主要执行流程 ====================

if __name__ == "__main__":
    print("\n" + "="*60)
    print("              Rv1819c 基因分类统计分析")
    print("="*60)
    
    # 1. 创建热图
    create_heatmap(abundance_data, OUT_PATH)
    
    # 2. 创建堆叠柱状图
    create_stacked_barplot(abundance_data, OUT_PATH)
    
    # 3. 创建总体统计图表
    create_overall_statistics(df_merge, OUT_PATH)
    
    # 4. 生成汇总统计表
    summary_stats, detailed_stats = create_summary_tables(df_merge, abundance_data, OUT_PATH)
    
    # 5. 创建分类比例对比图
    create_proportion_comparison(df_merge, OUT_PATH)
    
    # 6. 创建可视化报告总结
    print("\n" + "="*60)
    print("              可视化分析报告总结")
    print("="*60)
    print(f"数据集总规模: {len(df_merge):,} 个物种记录")
    print(f"分类级别: {len(taxonomic_levels)} 个 ({', '.join(taxonomic_levels)})")

    # 类别分布统计
    class_dist = df_merge['class'].value_counts()
    print("\n类别分布:")
    for cls, count in class_dist.items():
        percentage = count / len(df_merge) * 100
        print(f"  {cls}类: {count:,} ({percentage:.1f}%)")

    # 分类多样性统计
    print("\n分类多样性:")
    for level in taxonomic_levels:
        unique_count = df_merge[level].nunique()
        print(f"  {level}: {unique_count:,} 个不同类群")

    print("\n输出文件列表:")
    print(f"  合并数据: {OUT_PATH}/rv1819c_filter_meta.csv")
    print(f"  热图: {OUT_PATH}/taxonomic_abundance_heatmap.png")
    print(f"  堆叠柱状图: {OUT_PATH}/taxonomic_abundance_stacked_barplot.png")
    print(f"  总体统计图: {OUT_PATH}/overall_statistics.png")
    print(f"  比例对比图: {OUT_PATH}/proportion_comparison.png")
    print(f"  总体统计表: {OUT_PATH}/summary_statistics.csv")
    for level in taxonomic_levels:
        print(f"  {level}丰度数据: {OUT_PATH}/abundance_{level.lower()}.csv")
        print(f"  {level}详细统计: {OUT_PATH}/detailed_statistics_{level.lower()}.csv")

    # 显示总体统计（前几行）
    print("\n=== 总体汇总统计 ===")
    print(summary_stats.to_string(index=False))

    print("\n所有可视化图表和数据表格已成功生成！")
    print("="*60)