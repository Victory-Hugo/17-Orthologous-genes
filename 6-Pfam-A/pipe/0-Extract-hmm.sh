#!/bin/bash

set -e

CONF_FILE="/mnt/f/onedrive/文档（科研）/脚本/Download/17-Orthologous-genes/6-Pfam-A/conf/1-提取hmm.conf"

if [ ! -f "$CONF_FILE" ]; then
    echo "Error: config file not found: $CONF_FILE"
    exit 1
fi

# shellcheck source=/dev/null
source "$CONF_FILE"

if [ -z "${HMM_DB:-}" ] || [ -z "${OUTPUT_DIR:-}" ] || [ -z "${OUTPUT_FILE:-}" ]; then
    echo "Error: HMM_DB, OUTPUT_DIR, OUTPUT_FILE must be set in config"
    exit 1
fi

if [ -z "${PFAM_IDS:-}" ]; then
    echo "Error: PFAM_IDS must be set in config"
    exit 1
fi

# 检查输入文件是否存在
if [ ! -f "$HMM_DB" ]; then
    echo "Error: HMM database file not found: $HMM_DB"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"
: > "$OUTPUT_FILE"

extract_with_hmmfetch() {
    local pfam_id="$1"
    hmmfetch "$HMM_DB" "$pfam_id" >> "$OUTPUT_FILE"
}

extract_with_awk() {
    local pfam_id="$1"
    awk -v id="${pfam_id}" '
    BEGIN { found = 0; in_entry = 0 }
    /^HMMER3\/f/ { in_entry = 1 }
    in_entry && /^ACC   / {
        if ($0 ~ "ACC   " id) {
            found = 1
        } else if (found) {
            exit
        }
    }
    in_entry && /^NAME  / {
        if ($0 ~ "NAME  " id) {
            found = 1
        } else if (found) {
            exit
        }
    }
    found {
        print $0
        if (/^\/\//) {
            exit
        }
    }
    ' "$HMM_DB" >> "$OUTPUT_FILE"
}

if command -v hmmfetch &> /dev/null; then
    echo "Using hmmfetch to extract selected Pfam IDs..."
    for pfam_id in ${PFAM_IDS}; do
        extract_with_hmmfetch "$pfam_id"
    done
else
    echo "hmmfetch not found. Using manual extraction..."
    for pfam_id in ${PFAM_IDS}; do
        extract_with_awk "$pfam_id"
    done
fi

if [ -s "$OUTPUT_FILE" ]; then
    echo "Output file created: $OUTPUT_FILE"
    echo "Output file size: $(wc -l < "$OUTPUT_FILE") lines"

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
    echo "Warning: output file is empty"
    exit 1
fi
