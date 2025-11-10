#!/bin/bash
#? 统计多个FASTA文件中序列的长度（含'-'和不含'-'），并输出到汇总文件中，用来计算基因在不同物种中的长度变化情况。
FASTA_LIST_TXT="/home/luolintao/0_Github/17-Orthologous-genes/2-HMM/conf/4-搜索完_list.txt" 
SUMMARY_FILE="/home/luolintao/5-AB-Baoman/2-Lol家族基因搜索/data/序列长度统计.txt"

# 统计总数
TOTAL=$(grep -cve '^\s*$' "${FASTA_LIST_TXT}")
COUNT=0

# 清空旧结果并添加表头
echo -e "Sequence_ID\tLength_all\tLength_no_dash" > "${SUMMARY_FILE}"

# 定义进度条函数
progress_bar() {
    local current=$1
    local total=$2
    local width=40
    local progress=$((current * width / total))
    local percent=$((current * 100 / total))

    # 构建进度条，不再使用 cut
    local filled=$(printf "%0.s#" $(seq 1 $progress))
    local empty=$(printf "%0.s " $(seq $((progress + 1)) $width))
    printf "\r[%s%s] %3d%% (%d/%d)" "$filled" "$empty" "$percent" "$current" "$total"
}

# 主循环
while IFS= read -r FASTA; do
    ((COUNT++))
    # 跳过空行
    [[ -z "$FASTA" ]] && continue

    awk '
    /^>/ {
        if (seq) {
            len_all = length(seq)
            gsub(/-/, "", seq)
            len_no_dash = length(seq)
            print id "\t" len_all "\t" len_no_dash
        }
        id = substr($0, 2)
        seq = ""
        next
    }
    /^[^>]/ { seq = seq $0 }
    END {
        if (seq) {
            len_all = length(seq)
            gsub(/-/, "", seq)
            len_no_dash = length(seq)
            print id "\t" len_all "\t" len_no_dash
        }
    }' "${FASTA}" >> "${SUMMARY_FILE}"

    progress_bar "$COUNT" "$TOTAL"
done < "${FASTA_LIST_TXT}"

# 输出换行防止最后一行覆盖
echo -e "\nDone! Summary written to: ${SUMMARY_FILE}"
