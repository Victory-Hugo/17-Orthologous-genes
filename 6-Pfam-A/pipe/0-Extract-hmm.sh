#!/bin/bash

# 从 Pfam-A.hmm 中提取 PF00005 域
# Extract PF00005 HMM profile from Pfam-A.hmm

set -e

# 设置文件路径
HMM_DB="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/6-Pfam-A/data/Pfam-A.hmm"
OUTPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/6-Pfam-A/data"
OUTPUT_FILE="${OUTPUT_DIR}/PF00005.hmm"
PFAM_ACC="PF00005.34"  # PF00005 的完整 accession 号

# 检查输入文件是否存在
if [ ! -f "$HMM_DB" ]; then
    echo "Error: HMM database file not found: $HMM_DB"
    exit 1
fi

# 使用 hmmfetch 提取 PF00005（推荐方法，需要 HMMER 工具集）
if command -v hmmfetch &> /dev/null; then
    echo "Using hmmfetch to extract ${PFAM_ACC}..."
    hmmfetch "$HMM_DB" "$PFAM_ACC" > "$OUTPUT_FILE"
    if [ $? -eq 0 ]; then
        echo "Successfully extracted ${PFAM_ACC} to: $OUTPUT_FILE"
    else
        echo "Error: hmmfetch failed"
        exit 1
    fi
else
    # 备选方案：手动提取（如果没有安装 HMMER）
    echo "hmmfetch not found. Using manual extraction..."
    
    awk -v acc="${PFAM_ACC}" '
    BEGIN {
        found = 0
        in_entry = 0
    }
    /^HMMER3\/f/ {
        in_entry = 1
    }
    in_entry && /^ACC   / {
        if ($0 ~ "ACC   " acc) {
            found = 1
        } else if (found) {
            # 已经找完了目标条目
            exit
        }
    }
    found {
        print $0
        if (/^\/\//) {
            # 到达末尾标记，提取完成
            exit
        }
    }
    ' "$HMM_DB" > "$OUTPUT_FILE"
    
    echo "Successfully extracted ${PFAM_ACC} to: $OUTPUT_FILE"
fi

# 验证输出文件是否包含 PF00005
if [ -s "$OUTPUT_FILE" ] && grep -q "PF00005" "$OUTPUT_FILE"; then
    echo "Verification: ${PFAM_ACC} found in output file"
    echo "Output file size: $(wc -l < "$OUTPUT_FILE") lines"
    
    # 自动编译 HMM 数据库以供 hmmscan 使用
    if command -v hmmpress &> /dev/null; then
        echo ""
        echo "[INFO] 编译 HMM 数据库以供 hmmscan 使用..."
        hmmpress "$OUTPUT_FILE"
        if [ $? -eq 0 ]; then
            echo "[INFO] HMM 数据库编译完成"
        else
            echo "[WARNING] HMM 数据库编译失败，请手动运行: hmmpress $OUTPUT_FILE"
        fi
    fi
else
    echo "Warning: ${PFAM_ACC} not found in output file"
    exit 1
fi
