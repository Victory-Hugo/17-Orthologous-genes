#!/bin/bash

PY_SRC="/mnt/l/19-Rv1819c-Gene/4-Clinker/script/compute_taxonomy_heatmaps.p"
MATRIX="/mnt/l/19-Rv1819c-Gene/4-Clinker/output/run_20260106_124347/clinker/similarity_matrix.tsv"
META="/mnt/l/19-Rv1819c-Gene/4-Clinker/output/run_20260106_124347/clinker/2-存在Rv1819c_META.csv"

python3 ${PY_SRC} \
  --similarity ${MATRIX} \
  --metadata ${META}
