#!/bin/bash

# 脚本功能：将单拷贝基因文件夹从源目录合并到目标目录
# 说明：将源目录中的所有基因ID文件夹及其内的序列文件合并到目标目录中
# 作者：自动生成
# 日期：2025-11-04

# 定义源目录和目标目录
SOURCE_DIR="/mnt/f/15_Bam_Tam/2-物种树/output/all_single_copy"
TARGET_DIR="/mnt/f/15_Bam_Tam/5-补齐更多物种/output/all_single_copy"

# 颜色输出定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${YELLOW}=== 单拷贝基因文件夹合并脚本 ===${NC}"
echo ""

# 检查源目录是否存在
if [ ! -d "$SOURCE_DIR" ]; then
    echo -e "${RED}错误：源目录不存在 - $SOURCE_DIR${NC}"
    exit 1
fi

# 创建目标目录（如果不存在）
if [ ! -d "$TARGET_DIR" ]; then
    echo -e "${YELLOW}目标目录不存在，正在创建...${NC}"
    mkdir -p "$TARGET_DIR"
fi

# 显示源目录信息
SOURCE_GENE_COUNT=$(find "$SOURCE_DIR" -mindepth 1 -maxdepth 1 -type d | wc -l)
SOURCE_FILE_COUNT=$(find "$SOURCE_DIR" -type f | wc -l)
echo -e "${GREEN}源目录信息：${NC}"
echo "  路径：$SOURCE_DIR"
echo "  基因ID文件夹数：$SOURCE_GENE_COUNT"
echo "  序列文件数：$SOURCE_FILE_COUNT"
echo ""

# 显示目标目录信息
TARGET_GENE_COUNT_BEFORE=$(find "$TARGET_DIR" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l)
TARGET_FILE_COUNT_BEFORE=$(find "$TARGET_DIR" -type f 2>/dev/null | wc -l)
echo -e "${GREEN}目标目录合并前：${NC}"
echo "  路径：$TARGET_DIR"
echo "  基因ID文件夹数：$TARGET_GENE_COUNT_BEFORE"
echo "  序列文件数：$TARGET_FILE_COUNT_BEFORE"
echo ""

# 确认操作
read -p "是否继续合并? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo -e "${YELLOW}操作已取消${NC}"
    exit 0
fi

# 合并文件和文件夹
echo -e "${YELLOW}开始合并...${NC}"

# 使用 cp -r 将源目录中的所有基因ID文件夹复制到目标目录
# -u 选项：仅在源文件较新或目标文件不存在时复制
if cp -r "$SOURCE_DIR"/* "$TARGET_DIR/" 2>/dev/null; then
    echo -e "${GREEN}✓ 成功合并所有文件和文件夹${NC}"
else
    echo -e "${RED}✗ 合并文件时出错${NC}"
    exit 1
fi

# 统计合并后的文件数量
TARGET_GENE_COUNT_AFTER=$(find "$TARGET_DIR" -mindepth 1 -maxdepth 1 -type d | wc -l)
TARGET_FILE_COUNT_AFTER=$(find "$TARGET_DIR" -type f | wc -l)

echo ""
echo -e "${GREEN}=== 合并完成 ===${NC}"
echo -e "${GREEN}目标目录合并后：${NC}"
echo "  基因ID文件夹数：$TARGET_GENE_COUNT_AFTER （增加 $((TARGET_GENE_COUNT_AFTER - TARGET_GENE_COUNT_BEFORE))）"
echo "  序列文件数：$TARGET_FILE_COUNT_AFTER （增加 $((TARGET_FILE_COUNT_AFTER - TARGET_FILE_COUNT_BEFORE))）"
echo ""

# 询问是否删除源目录
read -p "是否删除源目录? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo -e "${YELLOW}正在删除源目录...${NC}"
    rm -rf "$SOURCE_DIR"
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ 源目录已删除${NC}"
    else
        echo -e "${RED}✗ 删除源目录失败${NC}"
    fi
else
    echo -e "${YELLOW}源目录保留${NC}"
fi

echo ""
echo -e "${GREEN}脚本执行完毕！${NC}"
