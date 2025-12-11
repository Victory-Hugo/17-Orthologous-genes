#!/bin/bash

# PoSeiDon Docker 镜像加载脚本
# 用于在没有网络的服务器上加载所有Docker镜像

set -e

echo "开始加载 PoSeiDon 所需的 Docker 镜像..."

# 检查docker_images目录是否存在
if [[ ! -d "docker_images" ]]; then
    echo "错误: 找不到 docker_images 目录"
    echo "请确保已将打包的镜像文件上传到此目录"
    exit 1
fi

# 加载所有tar文件
for tar_file in docker_images/*.tar; do
    if [[ -f "$tar_file" ]]; then
        echo "正在加载镜像: $tar_file"
        docker load -i "$tar_file"
        echo "已加载: $tar_file"
        echo "---"
    fi
done

echo "所有镜像加载完成！"
echo ""
echo "验证已加载的镜像:"
docker images | grep nanozoo

echo ""
echo "Docker 镜像部署完成！"
echo "现在可以运行 Nextflow 流水线了。"