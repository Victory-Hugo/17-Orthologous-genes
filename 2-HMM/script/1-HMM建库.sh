#!/bin/bash

INPUT_FASTA="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/2-HMM/data/lolD_merge.faa"
ALN_FASTA="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/2-HMM/data/lolD.aln.faa"
OUTPUT_HMM="/mnt/f/OneDrive/文档（科研）/脚本/Download/17-Orthologous-genes/2-HMM/data/LolD.hmm"



# 定义彩色输出函数
echo_cyan()   { echo -e "\033[1;36m$*\033[0m"; } 

mafft --auto --reorder --anysymbol \
    $INPUT_FASTA > $ALN_FASTA
    
trimal -gt 0.8 -in $ALN_FASTA -out ${ALN_FASTA%.faa}.trim.faa

hmmbuild "${OUTPUT_HMM}" "${ALN_FASTA%.faa}.trim.faa"

echo_cyan "HMM library LolD.hmm has been built successfully."



