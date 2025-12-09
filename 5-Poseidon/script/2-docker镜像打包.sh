#!/bin/bash

###############################################################################
# PoSeiDon Docker 镜像离线打包脚本
# 功能：下载 poseidon.nf 工作流所需的所有 Docker 镜像并打包
# 使用方法：./download_docker_images.sh
###############################################################################

set -e

# 定义输出目录
OUTPUT_DIR="docker_images"
ARCHIVE_NAME="poseidon_docker_images.tar.gz"
MANIFEST_FILE="docker_images_manifest.txt"

# 定义所有需要的 Docker 镜像
declare -a DOCKER_IMAGES=(
    "nanozoo/translatorx:1.1--1f1af23"
    "nanozoo/bioruby:2.0.1--1f8a188"
    "nanozoo/texlive-core:3.14159265--bb5506c"
    "nanozoo/raxml:8.2.12--27d10cf"
    "nanozoo/newick_utils:1.6--014d613"
    "nanozoo/hyphy:2.5.14--ed56f29"
    "nanozoo/hyphy:2.2.7--9dc0fe5"
    "nanozoo/paml:4.9--cff12e2"
)

echo "=================================="
echo "PoSeiDon Docker 镜像离线打包工具"
echo "=================================="
echo ""

# 检查 Docker 是否安装
if ! command -v docker &> /dev/null; then
    echo "❌ 错误: Docker 未安装或不在 PATH 中"
    echo "请先安装 Docker"
    exit 1
fi

# 检查 Docker daemon 是否运行
if ! docker ps &> /dev/null; then
    echo "❌ 错误: Docker daemon 未运行"
    echo "请先启动 Docker"
    exit 1
fi

echo "✓ Docker 已安装并运行"
echo ""

# 创建输出目录
mkdir -p "$OUTPUT_DIR"
echo "✓ 创建输出目录: $OUTPUT_DIR"
echo ""

# 初始化清单文件
> "$MANIFEST_FILE"
echo "开始下载 Docker 镜像..." > "$MANIFEST_FILE"
echo "生成时间: $(date)" >> "$MANIFEST_FILE"
echo "" >> "$MANIFEST_FILE"
echo "镜像列表:" >> "$MANIFEST_FILE"

echo "开始下载 Docker 镜像..."
echo ""

# 下载每个镜像并保存为 tar 文件
for image in "${DOCKER_IMAGES[@]}"; do
    echo "正在处理: $image"
    
    # 将镜像名称转换为文件名（替换特殊字符）
    filename=$(echo "$image" | tr '/:' '_' | tr '.' '_')
    tar_file="$OUTPUT_DIR/${filename}.tar"
    
    # 拉取镜像
    echo "  [1/3] 拉取镜像..."
    if docker pull "$image"; then
        echo "  ✓ 镜像拉取成功"
    else
        echo "  ✗ 镜像拉取失败: $image"
        continue
    fi
    
    # 保存镜像为 tar 文件
    echo "  [2/3] 保存镜像为 tar 文件..."
    if docker save -o "$tar_file" "$image"; then
        echo "  ✓ 镜像已保存: $tar_file"
    else
        echo "  ✗ 镜像保存失败: $image"
        rm -f "$tar_file"
        continue
    fi
    
    # 获取文件大小
    file_size=$(du -h "$tar_file" | cut -f1)
    echo "  [3/3] 文件大小: $file_size"
    echo "  ✓ 完成: $image"
    
    # 添加到清单
    echo "  - $image (文件: ${filename}.tar, 大小: $file_size)" >> "$MANIFEST_FILE"
    echo ""
done

echo ""
echo "=================================="
echo "镜像下载完成"
echo "=================================="
echo ""

# 显示已下载的镜像
echo "已下载的镜像文件:"
ls -lh "$OUTPUT_DIR"/*.tar 2>/dev/null | awk '{print "  " $9 " (" $5 ")"}'
echo ""

# 计算总大小
total_size=$(du -sh "$OUTPUT_DIR" | cut -f1)
echo "总大小: $total_size"
echo ""

# 打包所有镜像
echo "正在打包所有镜像到 $ARCHIVE_NAME..."
tar -czf "$ARCHIVE_NAME" "$OUTPUT_DIR/" "$MANIFEST_FILE"
archive_size=$(du -h "$ARCHIVE_NAME" | cut -f1)
echo "✓ 打包完成: $ARCHIVE_NAME ($archive_size)"
echo ""

# 生成导入脚本
cat > load_docker_images.sh << 'EOF'
#!/bin/bash

###############################################################################
# Docker 镜像导入脚本
# 功能：在离线环境中导入 Docker 镜像
# 使用方法：./load_docker_images.sh
###############################################################################

set -e

# 检查 Docker 是否安装
if ! command -v docker &> /dev/null; then
    echo "❌ 错误: Docker 未安装或不在 PATH 中"
    exit 1
fi

# 检查 Docker daemon 是否运行
if ! docker ps &> /dev/null; then
    echo "❌ 错误: Docker daemon 未运行"
    exit 1
fi

echo "=================================="
echo "Docker 镜像导入工具"
echo "=================================="
echo ""

OUTPUT_DIR="docker_images"

if [ ! -d "$OUTPUT_DIR" ]; then
    echo "❌ 错误: 找不到目录 $OUTPUT_DIR"
    echo "请确保你在正确的目录中运行此脚本"
    exit 1
fi

echo "开始导入 Docker 镜像..."
echo ""

# 计数器
count=0
success=0

for tar_file in "$OUTPUT_DIR"/*.tar; do
    if [ -f "$tar_file" ]; then
        count=$((count + 1))
        filename=$(basename "$tar_file")
        echo "[$count] 导入: $filename"
        
        if docker load -i "$tar_file"; then
            echo "  ✓ 导入成功"
            success=$((success + 1))
        else
            echo "  ✗ 导入失败"
        fi
        echo ""
    fi
done

echo "=================================="
echo "导入完成"
echo "=================================="
echo "成功导入: $success / $count 个镜像"
echo ""

# 列出已导入的镜像
echo "已加载的镜像:"
docker images | grep nanozoo || echo "  (未找到 nanozoo 镜像)"

EOF

chmod +x load_docker_images.sh
echo "✓ 已生成导入脚本: load_docker_images.sh"
echo ""

# 生成 README 文档
cat > OFFLINE_SETUP_README.md << 'EOF'
# PoSeiDon 离线环境配置指南

## 概述

本指南说明如何在没有网络连接的服务器上使用 PoSeiDon 工作流。

## 文件说明

- `poseidon_docker_images.tar.gz` - 所有 Docker 镜像的打包文件
- `load_docker_images.sh` - Docker 镜像导入脚本
- `docker_images_manifest.txt` - 镜像清单和元数据
- `docker_images/` - Docker 镜像目录（解压后）

## 所需工具

在离线环境的服务器上需要安装：
- Docker（已配置并可正常运行）
- Nextflow（>=23.01.0）
- 足够的磁盘空间（至少 50GB 用于镜像）

## 步骤 1: 准备镜像

### 在有网络的机器上：
```bash
# 下载镜像（耗时较长，取决于网络速度）
./download_docker_images.sh

# 生成的文件：
# - poseidon_docker_images.tar.gz  (~30-40GB)
# - docker_images/                 (镜像文件目录)
# - load_docker_images.sh          (导入脚本)
# - docker_images_manifest.txt     (清单文件)
```

### 传输到离线环境：
```bash
# 方式1: 使用 U 盘或移动硬盘
# 方式2: 使用 rsync 通过中间机器传输
# 方式3: 使用其他文件传输方式

# 建议：只传输 poseidon_docker_images.tar.gz 文件 (较小)
# 然后在目标服务器上解压
```

## 步骤 2: 在离线服务器上导入镜像

```bash
# 1. 解压镜像文件（如果还未解压）
tar -xzf poseidon_docker_images.tar.gz

# 2. 导入所有镜像
chmod +x load_docker_images.sh
./load_docker_images.sh

# 3. 验证镜像是否已导入
docker images | grep nanozoo
```

预期输出应该列出以下镜像：
```
REPOSITORY                  TAG
nanozoo/translatorx         1.1--1f1af23
nanozoo/bioruby             2.0.1--1f8a188
nanozoo/texlive-core        3.14159265--bb5506c
nanozoo/raxml               8.2.12--27d10cf
nanozoo/newick_utils        1.6--014d613
nanozoo/hyphy               2.5.14--ed56f29
nanozoo/hyphy               2.2.7--9dc0fe5
nanozoo/paml                4.9--cff12e2
```

## 步骤 3: 修改 Nextflow 配置

在你的 Nextflow 工作目录中，创建或修改 `nextflow.config`，使用本地镜像而不从网络拉取：

```groovy
// 在 configs/container.config 中确保以下配置存在
process {   
    withLabel: translatorx  { container = 'nanozoo/translatorx:1.1--1f1af23' } 
    withLabel: bioruby      { container = 'nanozoo/bioruby:2.0.1--1f8a188' }
    withLabel: tex          { container = 'nanozoo/texlive-core:3.14159265--bb5506c' }
    withLabel: raxml        { container = 'nanozoo/raxml:8.2.12--27d10cf' }
    withLabel: newick_utils { container = 'nanozoo/newick_utils:1.6--014d613' }
    withLabel: hyphy        { container = 'nanozoo/hyphy:2.5.14--ed56f29' }
    withLabel: gard         { container = 'nanozoo/hyphy:2.2.7--9dc0fe5' }
    withLabel: codeml       { container = 'nanozoo/paml:4.9--cff12e2' }
}
```

或者在命令行中禁用镜像拉取：

```bash
# 确保 Docker 守护进程运行
sudo systemctl start docker

# 运行 PoSeiDon（仅使用本地镜像）
nextflow run poseidon.nf \
    --fasta '*.fasta' \
    -profile local,docker \
    -with-docker
```

## 步骤 4: 运行工作流

```bash
# 确保 Docker 运行
sudo systemctl start docker

# 运行 PoSeiDon 工作流
nextflow run poseidon.nf \
    --fasta 'input/*.fasta' \
    --output results \
    -profile standard \
    -with-docker
```

## 常见问题

### Q1: Docker 镜像大小太大，无法传输
**A:** 你可以选择性地只传输某些镜像。例如：
```bash
# 只传输特定镜像
tar -czf minimal_images.tar.gz docker_images/nanozoo_translatorx*.tar
```

### Q2: 导入镜像时出现 "docker daemon" 错误
**A:** 启动 Docker 服务：
```bash
sudo systemctl start docker
# 或者
sudo service docker start
```

### Q3: 镜像导入后仍然无法找到
**A:** 检查 Docker 是否完全导入：
```bash
docker images | grep nanozoo
docker image inspect nanozoo/translatorx:1.1--1f1af23
```

### Q4: 运行工作流时仍然尝试拉取镜像
**A:** 确保在 Nextflow 命令中使用了正确的 profile：
```bash
# 检查配置
cat nextflow.config | grep -A 5 "docker"

# 尝试禁用网络拉取（编辑 docker 配置）
```

## 性能建议

1. **磁盘空间**: 确保有至少 50GB 的可用空间用于镜像
2. **内存**: PoSeiDon 对某些步骤需要较多内存，建议至少 16GB
3. **CPU**: 使用 `-profile` 指定本地 CPU 核心数

## 技术细节

### Docker 镜像来源
所有镜像来自 [NanoZoo](https://github.com/NanoZoo) 项目，一个为生物信息学提供的 Docker 镜像库。

### 镜像大小参考
- translatorx: ~1.2GB
- bioruby: ~700MB
- texlive-core: ~3GB
- raxml: ~800MB
- newick_utils: ~300MB
- hyphy: ~1.5GB
- paml: ~500MB

### 磁盘空间计算
- 单个镜像：总计 ~8-10GB
- 未压缩的 docker_images/ 目录：~8-10GB
- 压缩后的 tar.gz：~3-5GB

## 获取帮助

如有问题，请参考：
1. PoSeiDon GitHub: https://github.com/hoelzer/poseidon
2. Nextflow 文档: https://www.nextflow.io/docs/latest/
3. Docker 文档: https://docs.docker.com/

EOF

echo "✓ 已生成配置指南: OFFLINE_SETUP_README.md"
echo ""

# 最终总结
echo "=================================="
echo "✓ 所有准备工作已完成！"
echo "=================================="
echo ""
echo "接下来的步骤："
echo "1. 传输以下文件到离线服务器："
echo "   - $ARCHIVE_NAME"
echo "   - load_docker_images.sh"
echo "   - OFFLINE_SETUP_README.md"
echo ""
echo "2. 在离线服务器上运行："
echo "   tar -xzf $ARCHIVE_NAME"
echo "   chmod +x load_docker_images.sh"
echo "   ./load_docker_images.sh"
echo ""
echo "3. 详细步骤请查看: OFFLINE_SETUP_README.md"
echo ""
