#!/bin/bash

# GTDB-Tk 快速配置脚本
# 使用方法: bash 2-快速配置.sh

set -e

echo "========================================"
echo "GTDB-Tk 快速配置脚本"
echo "========================================"

# 定义路径
DATA_DIR="/mnt/d/3-GTDB-Database/data"
ARCHIVE_FILE="$DATA_DIR/gtdbtk_data.tar.gz"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 检查压缩文件是否存在
if [ ! -f "$ARCHIVE_FILE" ]; then
    echo "❌ 错误: 找不到 $ARCHIVE_FILE"
    exit 1
fi

echo ""
echo "📦 步骤 1: 检查压缩文件"
echo "---"
ls -lh "$ARCHIVE_FILE"
echo "✅ 文件存在"

echo ""
echo "📦 步骤 2: 解压数据库（这可能需要几分钟）"
echo "---"
cd "$DATA_DIR"
echo "正在解压到: $DATA_DIR"

# 检查是否已解压
if ls release*/files 1> /dev/null 2>&1; then
    echo "⚠️  数据库已解压"
    RELEASE_DIR=$(ls -d release*/ | head -n 1)
else
    echo "正在解压..."
    tar -xzf "$ARCHIVE_FILE"
    RELEASE_DIR=$(ls -d release*/ | head -n 1)
    echo "✅ 解压完成"
fi

RELEASE_PATH="${DATA_DIR}/${RELEASE_DIR%/}"
echo "📍 数据库路径: $RELEASE_PATH"

echo ""
echo "📦 步骤 3: 验证数据库完整性"
echo "---"

# 检查必要的文件
required_files=(
    "files/bac120.ms"
    "files/ar122.ms"
    "metadata/bac120_metadata.tsv"
    "metadata/ar122_metadata.tsv"
)

all_exists=true
for file in "${required_files[@]}"; do
    if [ -f "$RELEASE_PATH/$file" ]; then
        echo "✅ $file"
    else
        echo "❌ 缺失: $file"
        all_exists=false
    fi
done

if [ "$all_exists" = false ]; then
    echo ""
    echo "⚠️  警告: 某些文件可能缺失"
fi

echo ""
echo "📦 步骤 4: 设置环境变量"
echo "---"

# 创建环境变量设置脚本
ENV_SCRIPT="$SCRIPT_DIR/set_gtdbtk_env.sh"
cat > "$ENV_SCRIPT" << EOF
#!/bin/bash
# GTDB-Tk 环境变量设置脚本

export GTDBTK_DATA_PATH="$RELEASE_PATH"

echo "✅ GTDBTK_DATA_PATH 已设置为:"
echo "   \$GTDBTK_DATA_PATH = \$GTDBTK_DATA_PATH"
EOF

chmod +x "$ENV_SCRIPT"
echo "✅ 环境脚本已创建: $ENV_SCRIPT"

echo ""
echo "📦 步骤 5: 永久配置（可选）"
echo "---"
echo "将以下内容添加到 ~/.bashrc 或 ~/.bash_profile："
echo ""
echo "export GTDBTK_DATA_PATH=\"$RELEASE_PATH\""
echo ""
echo "然后执行: source ~/.bashrc"
echo ""

echo ""
echo "📦 步骤 6: 安装 GTDB-Tk (如果尚未安装)"
echo "---"

# 检查是否已安装 gtdbtk
if command -v gtdbtk &> /dev/null; then
    echo "✅ GTDB-Tk 已安装"
    gtdbtk --version
else
    echo "❌ 未检测到 GTDB-Tk"
    echo ""
    echo "请运行以下命令安装:"
    echo ""
    echo "  conda create -n gtdbtk python=3.10 -y"
    echo "  conda activate gtdbtk"
    echo "  conda install -c bioconda gtdbtk hmmer pplacer fasttree -y"
    echo ""
fi

echo ""
echo "========================================"
echo "✅ 配置完成！"
echo "========================================"
echo ""
echo "📋 后续使用步骤:"
echo ""
echo "1. 激活 GTDB-Tk 环境:"
echo "   conda activate gtdbtk"
echo ""
echo "2. 设置环境变量 (每次使用前):"
echo "   source $ENV_SCRIPT"
echo ""
echo "3. 验证安装:"
echo "   gtdbtk check_install_dir --install_dir \$GTDBTK_DATA_PATH"
echo ""
echo "4. 运行分类:"
echo "   gtdbtk classify_wf --genome_dir ./genomes --out_dir ./results --cpus 4"
echo ""
