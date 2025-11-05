#!/bin/bash

INPUT_ALN="/mnt/f/15_Bam_Tam/5-补齐更多物种/output/concatenated/supermatrix.faa"
OUTPUT_TREE="/mnt/f/15_Bam_Tam/5-补齐更多物种/output/concatenated/tree"
THREADS="16"
SOFTWARE="VeryFastTree"



$SOFTWARE \
    -threads "$THREADS" \
    < "$INPUT_ALN" >\
    "$OUTPUT_TREE"