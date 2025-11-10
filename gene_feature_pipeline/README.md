# Gene Feature Pipeline

This module contains a minimal pipeline for computing descriptive statistics for bacterial genes and their encoded proteins. Input coding sequences are provided as FASTA files in `input/`, and derived outputs (including example results) are placed inside `data/`.

## Directory layout

- `data/` – reference codon usage tables and generated feature summaries.
- `input/` – nucleotide FASTA files containing coding sequences to analyse.
- `script/` – processing logic implemented in Python.
- `pipe/` – lightweight orchestration scripts or workflow entry points.

## Usage

```bash
python pipe/run_pipeline.py
```

The pipeline invokes `script/compute_gene_features.py`, which can also be executed directly for custom datasets:

```bash
python script/compute_gene_features.py input/example_genes.fasta \
  --codon-usage data/reference_codon_usage.json \
  --output data/custom_features.csv
```

The resulting CSV summarises gene length, GC content, predicted protein length, isoelectric point, GRAVY hydropathy score, and Codon Adaptation Index (CAI).
