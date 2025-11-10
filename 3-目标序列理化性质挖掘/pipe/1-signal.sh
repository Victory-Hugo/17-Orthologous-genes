#!/bin/bash
# Author: BigLin
# Requires: bash, GNU parallel, signalp4(core 4.0+), coreutils

set -euo pipefail

FAA_LIST="/home/luolintao/0_Github/17-Orthologous-genes/3-目标序列理化性质挖掘/conf/signal_faa.txt"
SOFTWARE="signalp4"
JOBS=8

TEMP_DIR="/home/luolintao/5-AB-Baoman/2-Lol家族基因搜索/tmp"
OUTPUT_DIR="/home/luolintao/5-AB-Baoman/2-Lol家族基因搜索/output/1-signal"
JOB_LOG="/home/luolintao/5-AB-Baoman/2-Lol家族基因搜索/log/signalp_parallel.log"

mkdir -p "$TEMP_DIR" "$OUTPUT_DIR"

if [[ ! -s "$FAA_LIST" ]]; then
    echo "FAA list not found or empty: $FAA_LIST" >&2
    exit 1
fi

process_faa() {
    local faa_path="$1"
    local base_name
    base_name=$(basename "$faa_path")
    local stem="${base_name%.*}"
    local final_output="$OUTPUT_DIR/${stem}.txt"
    local temp_output
    temp_output=$(mktemp "$TEMP_DIR/${stem}.XXXXXX")

    cleanup() {
        rm -f "$temp_output"
    }
    trap cleanup EXIT

    "$SOFTWARE" \
        -f short \
        -t gram- \
        -T "$TEMP_DIR" \
        "$faa_path" > "$temp_output"

    mv -f "$temp_output" "$final_output"
    trap - EXIT
}
export SOFTWARE TEMP_DIR OUTPUT_DIR
export -f process_faa

parallel \
    --jobs "$JOBS" \
    --joblog "$JOB_LOG" \
    --resume-failed \
    process_faa {} :::: "$FAA_LIST"

#* 完成之后的清理工作
rm -rf "$TEMP_DIR"