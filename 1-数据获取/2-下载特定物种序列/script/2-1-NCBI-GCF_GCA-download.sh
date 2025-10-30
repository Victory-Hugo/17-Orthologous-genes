: '
脚本功能：
    本脚本用于下载 NCBI 的组装（Assembly）数据。
    通过调用指定的 Python 脚本，并传递配置文件和下载路径参数，实现自动化下载。

参数说明：
    PYTHON_SCRIPT  - 指定用于下载的 Python 脚本路径。
    TXT_FILE       - 包含下载信息的配置文件路径。
    PYTHON_PATH    - Python 解释器路径（建议使用虚拟环境）。
    DOWNLOAD_PATH  - 数据下载保存的目标路径。

使用说明：
    运行本脚本前，请确保上述路径均已正确设置，且所需 Python 环境及依赖已安装。
    执行脚本后，将自动根据配置文件内容下载 NCBI 组装数据至指定目录。
'
#! /bin/bash
set -euo pipefail
unset http_proxy
unset https_proxy
#* 下载 NCBI 数据的组装数据

PYTHON_SCRIPT="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-数据获取/2-下载特定物种序列/2-2-下载NCBI的Assembly.py"
TXT_FILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-数据获取/2-下载特定物种序列/conf/balanced_accessions_ID.txt"
PYTHON_PATH="/home/luolintao/miniconda3/envs/pyg/bin/python3"
DOWNLOAD_PATH="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-数据获取/2-下载特定物种序列/download/data/"
mkdir -p "${DOWNLOAD_PATH}"

${PYTHON_PATH} \
  ${PYTHON_SCRIPT} \
  --base-path "${DOWNLOAD_PATH}" \
  --file-path "${TXT_FILE}" \
  -w 5 --resume



