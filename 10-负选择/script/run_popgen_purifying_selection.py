#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Run purifying selection popgen metrics for full and fragments."""

import os
import re
import subprocess

BASE_DIR = "/mnt/f/onedrive/文档（科研）/脚本/Download/17-Orthologous-genes/10-负选择"
INPUT_DIR = os.path.join(BASE_DIR, "input")
SCRIPT = os.path.join(BASE_DIR, "script", "popgen_purifying_selection.py")
OUTPUT_DIR = os.path.join(BASE_DIR, "output")
OUTGROUP_MD = os.path.join(INPUT_DIR, "外群.md")

ALIGNMENTS = {
    "full": os.path.join(INPUT_DIR, "Rv1819c_merge_inter_rename_full_nt_ali.fasta"),
    "frag1": os.path.join(INPUT_DIR, "Rv1819c_merge_inter_rename_frag1_nt_ali.fasta"),
    "frag2": os.path.join(INPUT_DIR, "Rv1819c_merge_inter_rename_frag2_nt_ali.fasta"),
}


def read_outgroup(md_path):
    with open(md_path, "r", encoding="utf-8") as handle:
        text = handle.read()
    # Extract backticked accession
    m = re.search(r"`([^`]+)`", text)
    if not m:
        raise ValueError("Outgroup accession not found in 外群.md")
    return m.group(1).strip()


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    outgroup = read_outgroup(OUTGROUP_MD)

    summary_path = os.path.join(OUTPUT_DIR, "Rv1819c_popgen_summary.md")
    with open(summary_path, "w", encoding="utf-8") as summary:
        summary.write("# Rv1819c 群体遗传净化选择结果汇总\n\n")
        summary.write(f"Outgroup: `{outgroup}`\n\n")

        for label, fasta in ALIGNMENTS.items():
            if not os.path.exists(fasta):
                summary.write(f"- {label}: missing alignment {fasta}\n")
                continue
            out_tsv = os.path.join(OUTPUT_DIR, f"Rv1819c_popgen_{label}.tsv")
            cmd = ["python", SCRIPT, "--fasta", fasta, "--outgroup", outgroup, "--output", out_tsv]
            subprocess.run(cmd, check=True)

            summary.write(f"## {label}\n\n")
            summary.write(f"Output: `{out_tsv}`\n\n")
            # Read key metrics for quick summary
            metrics = {}
            with open(out_tsv, "r", encoding="utf-8") as handle:
                for line in handle:
                    if line.strip() and not line.startswith("metric"):
                        k, v = line.rstrip().split("\t", 1)
                        metrics[k] = v
            summary.write("- piN: {0}\n".format(metrics.get("pi_n")))
            summary.write("- piS: {0}\n".format(metrics.get("pi_s")))
            summary.write("- piN/piS: {0}\n".format(metrics.get("pi_n_over_pi_s")))
            summary.write("- MK fixed (N/S): {0}/{1}\n".format(metrics.get("mk_fixed_nonsyn"), metrics.get("mk_fixed_syn")))
            summary.write("- MK poly (N/S): {0}/{1}\n".format(metrics.get("mk_poly_nonsyn"), metrics.get("mk_poly_syn")))
            summary.write("- MK p-value (Fisher, greater): {0}\n\n".format(metrics.get("mk_pvalue")))


if __name__ == "__main__":
    main()
