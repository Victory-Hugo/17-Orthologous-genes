#!/bin/bash

# PoSeiDon Docker 镜像打包脚本
# 用于在有网络的环境中拉取并保存所有必需的Docker镜像

set -e

echo "开始拉取和保存 PoSeiDon 所需的 Docker 镜像..."

# 创建保存镜像的目录
mkdir -p docker_images

# 读取镜像列表并拉取保存
while IFS= read -r image; do
    if [[ -n "$image" && ! "$image" =~ ^#.* ]]; then
        echo "正在拉取镜像: $image"
        docker pull "$image"
        
        # 将镜像名称中的特殊字符替换为下划线，用作文件名
        safe_name=$(echo "$image" | sed 's/[\/:]/_/g' | sed 's/--/_/g')
        
        echo "正在保存镜像到文件: docker_images/${safe_name}.tar"
        docker save "$image" -o "docker_images/${safe_name}.tar"
        
        echo "已完成: $image"
        echo "---"
    fi
done < docker_images_list.txt

echo "所有镜像已成功拉取并保存到 docker_images/ 目录"
echo "打包完成的文件列表:"
ls -lh docker_images/

# 创建镜像映射文件，用于部署时的参考
echo "# 镜像名称到文件名的映射" > docker_images/image_mapping.txt
while IFS= read -r image; do
    if [[ -n "$image" && ! "$image" =~ ^#.* ]]; then
        safe_name=$(echo "$image" | sed 's/[\/:]/_/g' | sed 's/--/_/g')
        echo "$image -> ${safe_name}.tar" >> docker_images/image_mapping.txt
    fi
done < docker_images_list.txt

echo ""
echo "镜像映射文件已创建: docker_images/image_mapping.txt"