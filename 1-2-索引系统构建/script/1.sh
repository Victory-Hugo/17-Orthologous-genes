#!/usr/bin/env bash
# 常用命令示例（在项目根目录运行：${BASE}）

set -euo pipefail
BASE="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-2-索引系统构建"
BUILDPY="${BASE}/python/build_index.py"
EXTRACTPY="${BASE}/python/extract.py"
CONFIG="${BASE}/python/config.json"
INDEX="${BASE}/index"
TRANS="${BASE}/python/trans.py"

# 1) 增量构建索引（默认只处理新组装）
python ${BUILDPY} --config ${CONFIG} --workers 16

# 2) 重建指定组装（保留其它组装）
# python ${BUILDPY} --config ${CONFIG} --overwrite-existing --assemblies ${BASE}/database/GCA_009268105.2/..

# 3) 全量重建（清空 index 后重建，慎用）
# python ${BUILDPY} --config ${CONFIG} --force

# 4) 单条查询示例
# 默认 faa+fna 输出到终端
# python ${EXTRACTPY} --assembly GCA_009268105.2 --query dnaA --index ${INDEX}
# 只要蛋白，并写文件
# python ${EXTRACTPY} --assembly GCA_009268105.2 --query dnaA --format faa --out output.faa --index ${INDEX}

# 5) 批量查询示例
# 准备一个 CSV/TSV，列：Assembly,Query,Format[,OutFile]（表头可有可无）
# 将每行写入独立文件（自动按 组装号_query.后缀）
# python ${EXTRACTPY} --batch ${BASE}/script/queries_1.csv --index ${INDEX}
# 将多行追加到指定文件（按 OutFile 字段）
# python ${EXTRACTPY} --batch ${BASE}/script/queries_2.csv --index ${INDEX}

# python ${TRANS} \
#     --input ${BASE}/input/GCF_000006765.1.faa \
#     --output ${BASE}/input/GCF_000006765.1.fna

# python ${TRANS} \
#     --batch ${BASE}/script/queries_3.csv