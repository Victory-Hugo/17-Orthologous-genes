"""Pipeline entry point for computing gene feature summaries."""
from __future__ import annotations

import subprocess
from pathlib import Path


def run_pipeline() -> None:
    project_root = Path(__file__).resolve().parents[1]
    script_path = project_root / "script" / "compute_gene_features.py"
    input_fasta = project_root / "input" / "example_genes.fasta"
    codon_usage = project_root / "data" / "reference_codon_usage.json"
    output_csv = project_root / "data" / "example_gene_features.csv"

    command = [
        "python", str(script_path),
        str(input_fasta),
        "--codon-usage", str(codon_usage),
        "--output", str(output_csv),
    ]
    subprocess.run(command, check=True)


if __name__ == "__main__":
    run_pipeline()
