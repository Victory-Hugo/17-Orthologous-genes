#!/bin/bash
set -Eeuo pipefail

########################################
# 配置区（保持功能不变，纯顺序执行）
########################################
SCRIPT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/1-HMM"     # 项目根目录
LIST_FILE="${SCRIPT_DIR}/conf/1-prodigal_faa_list.txt"                                  # 样本蛋白列表
HMM_FILE_TXT="${SCRIPT_DIR}/conf/2-hmm库_list.txt"                                      # HMM 库列表
OUTPUT_DIR="/mnt/d/5-NCBI-Reference/hmm分析示例/output/"                                # 输出目录

CPU_PER_JOB="2"                                                                         # hmmsearch CPU
E_VALUE="1e-5"                                                                          # hmmsearch e-value 阈值
MATRIX_EVALUE_CUTOFF="1e-10"                                                            # 矩阵筛选阈值

MATRIX_SCRIPT="${SCRIPT_DIR}/script/3-生成存在缺失矩阵.py"                               # 矩阵生成脚本
HMM_SEARCH_SCRIPT="${SCRIPT_DIR}/script/2-备选-HMM比对.sh"                              # 顺序 hmmsearch 脚本
EXTRACT_SCRIPT="${SCRIPT_DIR}/script/4-提取匹配序列.py"                                  # 提取匹配序列
HITS_PYTHON_SCRIPT="${SCRIPT_DIR}/script/4-1-整理HMM命中列表.py"                         # 整理命中清单
HITS_QUERIES="FimH.aln"                                                                 # 目标 query 名

PRESENCE_MATRIX_TSV="${OUTPUT_DIR}/presence_absence_matrix.tsv"                         # 矩阵输出
MERGED_TBL_TSV="${OUTPUT_DIR}/tbl_merge.tsv"                                            # 合并 tbl 输出
HITS_OUTPUT_TSV="${OUTPUT_DIR}/${HITS_QUERIES}/hits_manifest.tsv"                       # 命中清单输出

########################################
# 1) hmmsearch
########################################
export LIST_FILE HMM_FILE_TXT OUTPUT_DIR CPU_PER_JOB E_VALUE MATRIX_SCRIPT PRESENCE_MATRIX_TSV MERGED_TBL_TSV MATRIX_EVALUE_CUTOFF
echo "[INFO] 调用 hmmsearch 顺序脚本: $HMM_SEARCH_SCRIPT"
if [[ -x "$HMM_SEARCH_SCRIPT" ]]; then
  "$HMM_SEARCH_SCRIPT"
else
  bash "$HMM_SEARCH_SCRIPT"
fi

########################################
# 2) 提取匹配序列 (直接 Python)
#    使用 only_one 模式，保留最高分
########################################
echo "[INFO] 提取匹配序列到: $OUTPUT_DIR"
while IFS= read -r faa_path || [[ -n "$faa_path" ]]; do
  [[ -z "$faa_path" ]] && continue
  cmd=(python3 "$EXTRACT_SCRIPT" "$faa_path" "$MERGED_TBL_TSV" "$OUTPUT_DIR" --mode only_one)
  echo "运行命令: ${cmd[*]}"
  if "${cmd[@]}"; then
    :
  else
    echo "[WARN] 提取失败或无命中: $faa_path"
  fi
done < "$LIST_FILE"

########################################
# 3) 整理命中清单 (直接 Python)
########################################
echo "[INFO] 整理命中清单: $HITS_OUTPUT_TSV"
python3 "$HITS_PYTHON_SCRIPT" \
  --hits-dir "$OUTPUT_DIR" \
  --output-tsv "$HITS_OUTPUT_TSV" \
  --queries "$HITS_QUERIES"
