#!/bin/bash


NEXT_FLOW="nextflow"                                                                                                   #! 版本 25.10.2
MAIN_SRC="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/5-Poseidon/poseidon.nf"                        #? 主脚本位置
INPUT_FILE="/mnt/c/Users/Administrator/Desktop/cleaned.fna"    #? 输入文件位置
OUT_BASE_DIR="/mnt/c/Users/Administrator/Desktop/"
OUTPUT_DIR="${OUT_BASE_DIR}/output"                                  #? 输出文件目录
WORKDIR="${OUT_BASE_DIR}/workdir"                                    #? 中间文件目录
CACHE_DIR="${OUT_BASE_DIR}/cachedir"                                 #? 缓存文件目录
ROOT="NA"                                                            #? 系统发育树定根
REFERENCE="NA"                                                       #? 与哪个序列比较

mkdir -p ${CACHE_DIR}
cd ${CACHE_DIR}
# #? 获取帮助
# nextflow run ${MAIN_SRC} --help

#? 第一次运行很慢要下载容器
 ${NEXT_FLOW} run \
    ${MAIN_SRC} \
    --fasta ${INPUT_FILE} \
    --output ${OUTPUT_DIR} \
    --cores 4 \
    --root "${ROOT}" \
    --reference "${REFERENCE}" \
    --memory "16" \
    --max_cores "16" \
    --translatorx_code "1" \
    --codeml_code "0" \
    --skip_gard "false" \
    --workdir "${WORKDIR}" \
    --cachedir "${CACHE_DIR}" \
    -resume


# 删除所有以 .nextflow.log 开头的文件
rm -f ./*nextflow.log*
