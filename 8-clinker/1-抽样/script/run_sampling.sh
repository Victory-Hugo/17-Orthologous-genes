#!/usr/bin/env bash
set -euo pipefail

# --n_clade 表示“在目标在 clade 内至少抽取条代表序列”的目标数量。

Rscript "$(dirname "$0")/sample_representatives.R" \
  --tree "/mnt/f/onedrive/文档（科研）/脚本/Download/17-Orthologous-genes/8-clinker/1-抽样/conf/Rv1819c_gene.aln.trimal.tree" \
  --meta "/mnt/f/onedrive/文档（科研）/脚本/Download/17-Orthologous-genes/8-clinker/1-抽样/conf/meta.tsv" \
  --outdir "/mnt/l/19-Rv1819c-Gene/4-Clinker/input" \
  --n 100 \
  --n_clade 50 \
  --seed 123 \
  --rv_tip "Rv1819c"
