#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Population genetics metrics for purifying selection from codon alignments.

Inputs: codon-aligned nucleotide FASTA with an outgroup sequence.
Outputs: summary metrics (piN, piS, piN/piS) and MK counts/p-value.
"""

import argparse
import itertools
import math
from collections import Counter

try:
    from scipy.stats import fisher_exact
except Exception as exc:  # pragma: no cover
    fisher_exact = None


CODON_TABLE = {
    # Standard genetic code
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

VALID_NT = set("ACGT")


def read_fasta(path):
    records = []
    name = None
    seq_parts = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq_parts).upper()))
                name = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line)
    if name is not None:
        records.append((name, "".join(seq_parts).upper()))
    return records


def is_valid_codon(codon):
    if len(codon) != 3:
        return False
    if any(nt not in VALID_NT for nt in codon):
        return False
    if CODON_TABLE.get(codon) == "*":
        return False
    return True


def codon_syn_sites(codon):
    """Number of synonymous sites in a codon (NG method)."""
    aa = CODON_TABLE.get(codon)
    if aa is None or aa == "*":
        return None
    syn_sites = 0.0
    for i in range(3):
        syn_changes = 0
        for nt in "ACGT":
            if nt == codon[i]:
                continue
            new_codon = codon[:i] + nt + codon[i + 1 :]
            aa2 = CODON_TABLE.get(new_codon)
            if aa2 is None or aa2 == "*":
                continue
            if aa2 == aa:
                syn_changes += 1
        syn_sites += syn_changes / 3.0
    return syn_sites


def codon_diff(c1, c2):
    """Expected synonymous/nonsynonymous differences between codons.

    Uses NG-style averaging across all shortest paths, skipping paths
    with stop-codon intermediates.
    """
    if c1 == c2:
        return 0.0, 0.0
    if not is_valid_codon(c1) or not is_valid_codon(c2):
        return None
    diffs = [i for i in range(3) if c1[i] != c2[i]]
    ndiff = len(diffs)
    if ndiff == 0:
        return 0.0, 0.0

    def step_count(path):
        codon = list(c1)
        syn = 0
        nonsyn = 0
        for pos in path:
            codon[pos] = c2[pos]
            new_codon = "".join(codon)
            aa_old = CODON_TABLE.get("".join(codon[:pos] + [c1[pos]] + codon[pos + 1 :]))
            aa_new = CODON_TABLE.get(new_codon)
            if aa_old is None or aa_new is None or aa_new == "*":
                return None
            if aa_new == aa_old:
                syn += 1
            else:
                nonsyn += 1
        return syn, nonsyn

    syn_total = 0.0
    nonsyn_total = 0.0
    valid_paths = 0

    for path in itertools.permutations(diffs):
        codon = list(c1)
        syn = 0
        nonsyn = 0
        valid = True
        for pos in path:
            old_codon = "".join(codon)
            codon[pos] = c2[pos]
            new_codon = "".join(codon)
            aa_old = CODON_TABLE.get(old_codon)
            aa_new = CODON_TABLE.get(new_codon)
            if aa_old is None or aa_new is None or aa_new == "*":
                valid = False
                break
            if aa_new == aa_old:
                syn += 1
            else:
                nonsyn += 1
        if valid:
            syn_total += syn
            nonsyn_total += nonsyn
            valid_paths += 1

    if valid_paths == 0:
        return None
    return syn_total / valid_paths, nonsyn_total / valid_paths


def calc_pi(seq_records, outgroup_id):
    outgroup = None
    ingroup = []
    for name, seq in seq_records:
        if name == outgroup_id:
            outgroup = (name, seq)
        else:
            ingroup.append((name, seq))

    if outgroup is None:
        raise ValueError(f"Outgroup {outgroup_id} not found in alignment.")

    if not ingroup:
        raise ValueError("No ingroup sequences found.")

    length = len(ingroup[0][1])
    for _, seq in ingroup:
        if len(seq) != length:
            raise ValueError("Sequences are not the same length.")
    if length % 3 != 0:
        raise ValueError("Alignment length is not a multiple of 3.")

    codon_count = length // 3

    # Pairwise piN/piS
    total_syn_sites = 0.0
    total_nonsyn_sites = 0.0
    total_syn_diffs = 0.0
    total_nonsyn_diffs = 0.0

    for (name1, seq1), (name2, seq2) in itertools.combinations(ingroup, 2):
        for i in range(codon_count):
            c1 = seq1[i * 3 : i * 3 + 3]
            c2 = seq2[i * 3 : i * 3 + 3]
            if not is_valid_codon(c1) or not is_valid_codon(c2):
                continue
            s1 = codon_syn_sites(c1)
            s2 = codon_syn_sites(c2)
            if s1 is None or s2 is None:
                continue
            syn_sites = (s1 + s2) / 2.0
            nonsyn_sites = 3.0 - syn_sites
            diff = codon_diff(c1, c2)
            if diff is None:
                continue
            syn_d, nonsyn_d = diff
            total_syn_sites += syn_sites
            total_nonsyn_sites += nonsyn_sites
            total_syn_diffs += syn_d
            total_nonsyn_diffs += nonsyn_d

    pi_s = total_syn_diffs / total_syn_sites if total_syn_sites > 0 else float("nan")
    pi_n = total_nonsyn_diffs / total_nonsyn_sites if total_nonsyn_sites > 0 else float("nan")
    pi_ratio = pi_n / pi_s if pi_s and not math.isnan(pi_s) else float("nan")

    # MK-style counts
    fixed_syn = 0
    fixed_nonsyn = 0
    poly_syn = 0
    poly_nonsyn = 0

    og_seq = outgroup[1]

    for i in range(codon_count):
        og_codon = og_seq[i * 3 : i * 3 + 3]
        if not is_valid_codon(og_codon):
            continue
        og_aa = CODON_TABLE[og_codon]
        ingroup_codons = []
        for _, seq in ingroup:
            codon = seq[i * 3 : i * 3 + 3]
            if is_valid_codon(codon):
                ingroup_codons.append(codon)
        if not ingroup_codons:
            continue

        allele_set = set(ingroup_codons)
        if len(allele_set) == 1:
            ing_codon = next(iter(allele_set))
            if ing_codon == og_codon:
                continue
            ing_aa = CODON_TABLE[ing_codon]
            if ing_aa == og_aa:
                fixed_syn += 1
            else:
                fixed_nonsyn += 1
        else:
            # Polymorphic site
            ing_aas = {CODON_TABLE[c] for c in allele_set}
            if len(ing_aas) == 1 and og_aa in ing_aas:
                poly_syn += 1
            else:
                poly_nonsyn += 1

    mk_pvalue = None
    if fisher_exact is not None:
        table = [[fixed_nonsyn, fixed_syn], [poly_nonsyn, poly_syn]]
        _, mk_pvalue = fisher_exact(table, alternative="greater")

    return {
        "ingroup_n": len(ingroup),
        "codon_count": codon_count,
        "pi_n": pi_n,
        "pi_s": pi_s,
        "pi_n_over_pi_s": pi_ratio,
        "syn_sites": total_syn_sites,
        "nonsyn_sites": total_nonsyn_sites,
        "syn_diffs": total_syn_diffs,
        "nonsyn_diffs": total_nonsyn_diffs,
        "mk_fixed_syn": fixed_syn,
        "mk_fixed_nonsyn": fixed_nonsyn,
        "mk_poly_syn": poly_syn,
        "mk_poly_nonsyn": poly_nonsyn,
        "mk_pvalue": mk_pvalue,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--outgroup", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    records = read_fasta(args.fasta)
    stats = calc_pi(records, args.outgroup)

    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        for key in [
            "ingroup_n",
            "codon_count",
            "pi_n",
            "pi_s",
            "pi_n_over_pi_s",
            "syn_sites",
            "nonsyn_sites",
            "syn_diffs",
            "nonsyn_diffs",
            "mk_fixed_syn",
            "mk_fixed_nonsyn",
            "mk_poly_syn",
            "mk_poly_nonsyn",
            "mk_pvalue",
        ]:
            val = stats.get(key)
            if isinstance(val, float):
                if math.isnan(val):
                    val_str = "NA"
                else:
                    val_str = f"{val:.6g}"
            else:
                val_str = str(val)
            handle.write(f"{key}\t{val_str}\n")


if __name__ == "__main__":
    main()
