#!/bin/bash

# 分析基因覆盖率和物种分布
SINGLE_COPY_DIR="/mnt/f/15_Bam_Tam/2-物种树/output/all_single_copy"
OUTPUT_DIR="/mnt/f/15_Bam_Tam/2-物种树/output/coverage_analysis"

mkdir -p "$OUTPUT_DIR"

echo "开始分析基因和物种的覆盖率..."

# 1. 统计每个基因的物种数量
echo "基因ID,物种数量" > "$OUTPUT_DIR/gene_species_count.csv"
for gene_dir in "$SINGLE_COPY_DIR"/*; do
    if [ -d "$gene_dir" ]; then
        gene_id=$(basename "$gene_dir")
        if [ "$gene_id" != "README.txt" ]; then
            species_count=$(ls "$gene_dir"/*.faa 2>/dev/null | wc -l)
            echo "$gene_id,$species_count" >> "$OUTPUT_DIR/gene_species_count.csv"
        fi
    fi
done

# 2. 获取所有物种列表
echo "收集所有物种..."
all_species_file="$OUTPUT_DIR/all_species.list"
find "$SINGLE_COPY_DIR" -name "*.faa" -exec basename {} .faa \; | sort | uniq > "$all_species_file"
total_species=$(wc -l < "$all_species_file")
echo "总物种数: $total_species"

# 3. 统计每个物种在多少个基因中出现
echo "物种名称,基因数量,覆盖率%" > "$OUTPUT_DIR/species_gene_count.csv"
while IFS= read -r species; do
    gene_count=0
    for gene_dir in "$SINGLE_COPY_DIR"/*; do
        if [ -d "$gene_dir" ]; then
            gene_id=$(basename "$gene_dir")
            if [ "$gene_id" != "README.txt" ] && [ -f "$gene_dir/${species}.faa" ]; then
                gene_count=$((gene_count + 1))
            fi
        fi
    done
    
    total_genes=$(ls "$SINGLE_COPY_DIR" | grep -v README.txt | wc -l)
    coverage=$(echo "scale=2; $gene_count * 100 / $total_genes" | bc)
    echo "$species,$gene_count,$coverage" >> "$OUTPUT_DIR/species_gene_count.csv"
    echo "  $species: $gene_count/$total_genes 基因 ($coverage%)"
done < "$all_species_file"

# 4. 生成基因-物种矩阵
echo "生成基因-物种存在矩阵..."
matrix_file="$OUTPUT_DIR/gene_species_matrix.csv"

# 创建表头
echo -n "Gene_ID" > "$matrix_file"
while IFS= read -r species; do
    echo -n ",$species" >> "$matrix_file"
done < "$all_species_file"
echo "" >> "$matrix_file"

# 填充矩阵
for gene_dir in "$SINGLE_COPY_DIR"/*; do
    if [ -d "$gene_dir" ]; then
        gene_id=$(basename "$gene_dir")
        if [ "$gene_id" != "README.txt" ]; then
            echo -n "$gene_id" >> "$matrix_file"
            while IFS= read -r species; do
                if [ -f "$gene_dir/${species}.faa" ]; then
                    echo -n ",1" >> "$matrix_file"
                else
                    echo -n ",0" >> "$matrix_file"
                fi
            done < "$all_species_file"
            echo "" >> "$matrix_file"
        fi
    fi
done

echo ""
echo "================ 分析完成 ================"
echo "总物种数: $total_species"
echo "总基因数: $(ls "$SINGLE_COPY_DIR" | grep -v README.txt | wc -l)"
echo ""
echo "输出文件:"
echo "  - 基因物种数统计: $OUTPUT_DIR/gene_species_count.csv"
echo "  - 物种基因数统计: $OUTPUT_DIR/species_gene_count.csv"
echo "  - 基因-物种矩阵: $OUTPUT_DIR/gene_species_matrix.csv"
echo "  - 所有物种列表: $OUTPUT_DIR/all_species.list"