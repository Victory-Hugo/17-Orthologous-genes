#!/bin/bash

###############################################################################
# 快速参考：如何在离线环境配置 PoSeiDon
###############################################################################

# ============================================================
# 第一步：在有网络的机器上下载 Docker 镜像
# ============================================================

# 进入 poseidon 目录
cd /mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/5-Poseidon

# 给脚本执行权限并运行
chmod +x download_docker_images.sh
./download_docker_images.sh

# 这会生成：
# - poseidon_docker_images.tar.gz (~3-5GB 压缩文件)
# - docker_images/ (解压后的镜像目录)
# - load_docker_images.sh (导入脚本)
# - OFFLINE_SETUP_README.md (详细说明)
# - docker_images_manifest.txt (清单)


# ============================================================
# 第二步：传输文件到离线服务器
# ============================================================

# 选项 1：传输压缩文件（推荐，文件较小）
# scp poseidon_docker_images.tar.gz user@server:/path/to/
# scp load_docker_images.sh user@server:/path/to/
# scp OFFLINE_SETUP_README.md user@server:/path/to/

# 选项 2：传输整个目录
# scp -r docker_images load_docker_images.sh OFFLINE_SETUP_README.md user@server:/path/to/


# ============================================================
# 第三步：在离线服务器上导入镜像
# ============================================================

# SSH 连接到服务器
ssh user@server

# 进入文件所在目录
cd /path/to/

# 如果传输的是压缩文件，先解压
tar -xzf poseidon_docker_images.tar.gz

# 确保脚本有执行权限
chmod +x load_docker_images.sh

# 启动 Docker（如果还未启动）
sudo systemctl start docker

# 导入所有镜像（耗时较长，取决于磁盘 I/O）
./load_docker_images.sh

# 验证镜像是否导入成功
docker images | grep nanozoo


# ============================================================
# 第四步：运行 PoSeiDon 工作流
# ============================================================

# 进入 poseidon 工作目录
cd /path/to/poseidon

# 方法 1：使用 Docker profile（推荐）
nextflow run poseidon.nf \
    --fasta 'input/*.fasta' \
    --output results \
    -profile standard

# 方法 2：使用本地 + Docker
nextflow run poseidon.nf \
    --fasta 'input/*.fasta' \
    --output results \
    -profile local,docker \
    --cores 8 \
    --max_cores 16

# 方法 3：使用条件和参数
nextflow run poseidon.nf \
    --fasta 'input/*.fasta' \
    --output results \
    --reference 'reference_species' \
    --outgroup 'outgroup_species' \
    -profile local,docker \
    --cores 6


# ============================================================
# 故障排除
# ============================================================

# 问题 1：Docker daemon 未运行
# 解决：
sudo systemctl start docker
sudo systemctl enable docker

# 问题 2：导入镜像失败
# 解决：检查磁盘空间和权限
df -h
docker ps  # 验证 Docker 可以正常工作

# 问题 3：Nextflow 仍然尝试从网络拉取镜像
# 解决：检查 nextflow.config 的 docker 配置
cat configs/container.config

# 问题 4：权限被拒绝
# 解决：
sudo usermod -aG docker $USER
newgrp docker
# 然后重新登录


# ============================================================
# 镜像大小参考（预期下载量）
# ============================================================

# 镜像名称                                  大小
# nanozoo/translatorx:1.1--1f1af23          ~1.2GB
# nanozoo/bioruby:2.0.1--1f8a188           ~700MB
# nanozoo/texlive-core:3.14159265--bb5506c ~3GB
# nanozoo/raxml:8.2.12--27d10cf            ~800MB
# nanozoo/newick_utils:1.6--014d613        ~300MB
# nanozoo/hyphy:2.5.14--ed56f29            ~1.5GB
# nanozoo/hyphy:2.2.7--9dc0fe5             ~1.5GB
# nanozoo/paml:4.9--cff12e2                ~500MB
#
# 总计（解压后）：~9GB
# 压缩后：~3-5GB


# ============================================================
# 额外说明
# ============================================================

# 1. 所有镜像来自 NanoZoo 项目
#    官方网站: https://github.com/NanoZoo/docker

# 2. 此方案不需要网络连接即可运行工作流
#    前提是：Docker 镜像已离线导入

# 3. 建议定期更新镜像（在有网络时）
#    这样可以获得最新的 bug fixes 和 improvements

# 4. 如果只需要某些镜像，可以修改脚本
#    只下载需要的部分，以节省时间和空间

