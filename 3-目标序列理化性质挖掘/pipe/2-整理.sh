#!/bin/bash
# Author: BigLin
# Requirements: bash, gawk (or awk supporting asorti), coreutils (printf, cat, wc)
# Description: Merge physicochemical TSV outputs into one table and drop columns that are entirely NA.

set -euo pipefail
OUTPUT_BASE="/mnt/l/18-Rv0194-Gene/3-SingleGene-tree/output/理化性质/"
LENGTH_TSV="${OUTPUT_BASE}/1-序列长度/length.tsv"
GC_TSV="${OUTPUT_BASE}/2-GC含量/gc_content.tsv"
PI_TSV="${OUTPUT_BASE}/3-等电点/isoelectric_point.tsv"
HYDRO_TSV="${OUTPUT_BASE}/4-疏水性指数/hydrophobicity.tsv"
CAI_TSV="${OUTPUT_BASE}/5-CAI/cai.tsv"
MERGED_TSV="${OUTPUT_BASE}/merged_physicochemical.tsv"

for file in "${LENGTH_TSV}" "${GC_TSV}" "${PI_TSV}" "${HYDRO_TSV}" "${CAI_TSV}"; do
    if [[ ! -f "${file}" ]]; then
        echo "缺少输入文件: ${file}" >&2
        exit 1
    fi
done

TOTAL_RECORDS=0
for file in "${LENGTH_TSV}" "${GC_TSV}" "${PI_TSV}" "${HYDRO_TSV}" "${CAI_TSV}"; do
    line_count=$(wc -l < "${file}")
    if (( line_count > 1 )); then
        TOTAL_RECORDS=$((TOTAL_RECORDS + line_count - 1))
    fi
done

PROGRESS_WIDTH=40

awk -F $'\t' -v OFS=$'\t' \
    -v LEN_FILE="${LENGTH_TSV}" \
    -v GC_FILE="${GC_TSV}" \
    -v PI_FILE="${PI_TSV}" \
    -v HYDRO_FILE="${HYDRO_TSV}" \
    -v CAI_FILE="${CAI_TSV}" \
    -v TOTAL_RECORDS="${TOTAL_RECORDS}" \
    -v PROGRESS_WIDTH="${PROGRESS_WIDTH}" \
    '
    function store(key, column_label, value) {
        if (value == "" || value == "NA") {
            value = "NA"
        } else {
            has_value[column_label] = 1
        }
        data[key SUBSEP column_label] = value
    }

    function render_progress(processed, total, width) {
        if (total <= 0) {
            return
        }
        progress_ratio = processed / total
        if (progress_ratio > 1) {
            progress_ratio = 1
        }
        percent = int(progress_ratio * 100)
        filled = int(progress_ratio * width + 0.5)
        if (filled > width) filled = width
        bar = ""
        for (i = 1; i <= width; i++) {
            bar = bar ((i <= filled) ? "#" : " ")
        }
        printf "\r[%s] %3d%% (%d/%d)", bar, percent, processed, total > "/dev/stderr"
        fflush("/dev/stderr")
    }

    FNR == 1 { next }  # skip headers

    {
        key = $1 OFS $2
        keys[key] = 1

        if (FILENAME == LEN_FILE) {
            store(key, "Length_all", $3)
            store(key, "Length_no_gap", $4)
        } else if (FILENAME == GC_FILE) {
            store(key, "GC_percent", $3)
        } else if (FILENAME == PI_FILE) {
            store(key, "Isoelectric_point", $3)
        } else if (FILENAME == HYDRO_FILE) {
            store(key, "Hydrophobicity", $3)
        } else if (FILENAME == CAI_FILE) {
            store(key, "CAI", $3)
        }

        processed++
        if (processed == TOTAL_RECORDS || processed % 100 == 0) {
            render_progress(processed, TOTAL_RECORDS, PROGRESS_WIDTH)
        }
    }

    END {
        metrics[1] = "Length_all"
        metrics[2] = "Length_no_gap"
        metrics[3] = "GC_percent"
        metrics[4] = "Isoelectric_point"
        metrics[5] = "Hydrophobicity"
        metrics[6] = "CAI"

        printf "Source_Fasta\tSequence_ID"
        metric_count = 0
        for (i = 1; i <= 6; i++) {
            column_label = metrics[i]
            if (has_value[column_label]) {
                metric_order[++metric_count] = column_label
                printf "\t%s", column_label
            }
        }
        printf "\n"

        n = asorti(keys, sorted_keys)
        for (i = 1; i <= n; i++) {
            split(sorted_keys[i], parts, OFS)
            printf "%s\t%s", parts[1], parts[2]
            for (j = 1; j <= metric_count; j++) {
                column_label = metric_order[j]
                value = data[sorted_keys[i] SUBSEP column_label]
                if (value == "") {
                    value = "NA"
                }
                printf "\t%s", value
            }
            printf "\n"
        }

        if (TOTAL_RECORDS > 0) {
            render_progress(TOTAL_RECORDS, TOTAL_RECORDS, PROGRESS_WIDTH)
            printf "\n" > "/dev/stderr"
        }
    }
    ' \
    "${LENGTH_TSV}" "${GC_TSV}" "${PI_TSV}" "${HYDRO_TSV}" "${CAI_TSV}" > "${MERGED_TSV}"

echo "合并完成: ${MERGED_TSV}"
