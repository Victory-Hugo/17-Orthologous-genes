#!/bin/bash

# 定义输入和输出目录
BUSCO_OUTPUT_DIR="/mnt/f/15_Bam_Tam/5-补齐更多物种/output/busco_out"
OUTPUT_DIR="/mnt/f/15_Bam_Tam/5-补齐更多物种/output/all_single_copy"

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 获取所有单拷贝基因的ID（从第一个样本提取）
FIRST_SAMPLE=$(ls "$BUSCO_OUTPUT_DIR" | head -n 1)
SINGLE_COPY_DIR="$BUSCO_OUTPUT_DIR/$FIRST_SAMPLE/run_bacteria_odb10/busco_sequences/single_copy_busco_sequences"

if [ ! -d "$SINGLE_COPY_DIR" ]; then
    echo "ERROR: Cannot find single copy sequences directory in $SINGLE_COPY_DIR"
    exit 1
fi

# 获取所有单拷贝基因的ID
GENE_IDS=$(ls "$SINGLE_COPY_DIR" | sed 's/.faa$//' | sort)

echo "Found $(echo "$GENE_IDS" | wc -w) single copy genes"
echo "Creating directories and collecting sequences..."

# 为每个基因创建目录并收集序列
for gene_id in $GENE_IDS; do
    gene_dir="$OUTPUT_DIR/$gene_id"
    mkdir -p "$gene_dir"
    
    # 遍历所有样本，收集该基因的序列
    for sample_dir in "$BUSCO_OUTPUT_DIR"/*; do
        if [ -d "$sample_dir" ]; then
            sample_name=$(basename "$sample_dir")
            single_copy_file="$sample_dir/run_bacteria_odb10/busco_sequences/single_copy_busco_sequences/${gene_id}.faa"
            
            # 如果该样本有这个单拷贝基因，复制到目标目录并重命名为样本名
            if [ -f "$single_copy_file" ]; then
                cp "$single_copy_file" "$gene_dir/${sample_name}.faa"
                echo "  ✓ $gene_id: $sample_name"
            fi
        fi
    done
done

echo ""
echo "Done! Single copy sequences collected in: $OUTPUT_DIR"
echo "Statistics:"
for gene_dir in "$OUTPUT_DIR"/*; do
    gene_id=$(basename "$gene_dir")
    count=$(ls "$gene_dir"/*.faa 2>/dev/null | wc -l)
    echo "  $gene_id: $count sequences"
done
