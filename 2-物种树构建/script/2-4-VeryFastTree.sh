#!/bin/bash

INPUT_ALN="/mnt/f/15_Bam_Tam/5-补齐更多物种/output/concatenated/supermatrix.faa"
OUTPUT_TREE="/mnt/f/15_Bam_Tam/5-补齐更多物种/output/concatenated/tree"
THREADS="16"
SOFTWARE="VeryFastTree"



$SOFTWARE \
    -threads "$THREADS" \
    < "$INPUT_ALN" >\
    "$OUTPUT_TREE"


python3 /home/luolintao/test_mail.py \
    "2-物种树构建/script/2-4-VeryFastTree.sh任务完成通知"\
    "<p>2-物种树构建/script/2-4-VeryFastTree.sh分析已完成，请查看结果目录。</p>"
