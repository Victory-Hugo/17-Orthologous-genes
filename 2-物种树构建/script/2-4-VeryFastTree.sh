#!/bin/bash

INPUT_ALN="/home/luolintao/0-tmp/3-Bam_Tam/output/物种树/supermatrix.faa"
OUTPUT_TREE="/home/luolintao/0-tmp/3-Bam_Tam/output/物种树/tree"
THREADS="64"
SOFTWARE="VeryFastTree"



$SOFTWARE \
    -threads "$THREADS" \
    < "$INPUT_ALN" >\
    "$OUTPUT_TREE"


python3 /home/luolintao/test_mail.py "2-物种树构建/script/2-4-VeryFastTree.sh任务完成通知" "<p>2-物种树构建/script/2-4-VeryFastTree.sh分析已完成，请查看结果目录。</p>"
