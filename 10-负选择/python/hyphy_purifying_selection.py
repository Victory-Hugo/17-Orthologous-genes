#!/usr/bin/env python3
import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
from collections import Counter
from datetime import datetime

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


class ConfigError(RuntimeError):
    pass


def load_config(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def ensure_file(path, desc):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise FileNotFoundError(f"{desc} 不存在或为空: {path}")


def parse_fasta(path):
    name = None
    seq = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if name is not None:
            yield name, "".join(seq)


def consensus_aa_by_codon(aln_path):
    seqs = [seq for _, seq in parse_fasta(aln_path)]
    if not seqs:
        return []
    length = len(seqs[0])
    for s in seqs:
        if len(s) != length:
            raise ValueError("对齐序列长度不一致")
    if length % 3 != 0:
        raise ValueError("对齐长度不是 3 的倍数")
    codon_count = length // 3
    consensus = []
    for i in range(codon_count):
        codons = []
        for s in seqs:
            codon = s[i * 3 : i * 3 + 3].upper()
            if "-" in codon or "N" in codon or len(codon) < 3:
                continue
            codons.append(codon)
        if not codons:
            consensus.append("")
            continue
        aas = [CODON_TABLE.get(c, "") for c in codons]
        aas = [a for a in aas if a and a != "*"]
        if not aas:
            consensus.append("")
            continue
        aa = Counter(aas).most_common(1)[0][0]
        consensus.append(aa)
    return consensus


def headers_to_index(headers, key_name):
    for i, h in enumerate(headers):
        if h[0] == key_name:
            return i
    return None


def fel_table(obj):
    headers = obj.get("MLE", {}).get("headers", [])
    rows = obj.get("MLE", {}).get("content", {}).get("0", [])
    idx_alpha = headers_to_index(headers, "alpha")
    idx_beta = headers_to_index(headers, "beta")
    idx_p = headers_to_index(headers, "p-value")
    if idx_alpha is None or idx_beta is None or idx_p is None:
        raise ValueError("FEL JSON 缺少关键字段（alpha/beta/p-value）")
    return rows, idx_alpha, idx_beta, idx_p


def fubar_table(obj):
    headers = obj.get("MLE", {}).get("headers", [])
    rows = obj.get("MLE", {}).get("content", {}).get("0", [])
    idx_prob = headers_to_index(headers, "Prob[alpha>beta]")
    if idx_prob is None:
        raise ValueError("FUBAR JSON 缺少 Prob[alpha>beta] 字段")
    return rows, idx_prob


def slac_table(obj):
    headers = obj.get("MLE", {}).get("headers", [])
    by_site = obj.get("MLE", {}).get("content", {}).get("0", {}).get("by-site", {})
    rows = by_site.get("AVERAGED", [])
    idx_dn = headers_to_index(headers, "dN")
    idx_ds = headers_to_index(headers, "dS")
    idx_p = headers_to_index(headers, "P [dN/dS < 1]")
    if idx_dn is None or idx_ds is None or idx_p is None:
        raise ValueError("SLAC JSON 缺少关键字段（dN/dS/P [dN/dS < 1]）")
    return rows, idx_dn, idx_ds, idx_p


def fmt_value(value):
    if value is None:
        return ""
    return value


def log_has_success(log_path):
    if not os.path.isfile(log_path):
        return False
    with open(log_path, "r", encoding="utf-8") as f:
        for line in f:
            if "STATUS=SUCCESS" in line:
                return True
    return False


def write_summary(ali_path, fel_path, fubar_path, slac_path, out_path):
    fel_obj = load_json(fel_path)
    fubar_obj = load_json(fubar_path)
    slac_obj = load_json(slac_path)

    fel_rows, fel_alpha_i, fel_beta_i, fel_p_i = fel_table(fel_obj)
    fubar_rows, fubar_pneg_i = fubar_table(fubar_obj)
    slac_rows, slac_dn_i, slac_ds_i, slac_p_i = slac_table(slac_obj)

    aa_list = consensus_aa_by_codon(ali_path)
    codon_count = len(aa_list)

    if not (len(fel_rows) == len(fubar_rows) == len(slac_rows) == codon_count):
        raise ValueError(
            "位点数不一致: FEL=%d, FUBAR=%d, SLAC=%d, AA=%d"
            % (len(fel_rows), len(fubar_rows), len(slac_rows), codon_count)
        )

    with open(out_path, "w", encoding="utf-8") as out:
        out.write(
            "codon_position\tamino_acid\tFEL_dN\tFEL_dS\tFEL_pvalue\t"
            "FUBAR_posterior_negative\tSLAC_dN\tSLAC_dS\tSLAC_pvalue\n"
        )
        for i in range(codon_count):
            pos = i + 1
            aa = aa_list[i]
            fel_dn = fel_rows[i][fel_beta_i]
            fel_ds = fel_rows[i][fel_alpha_i]
            fel_p = fel_rows[i][fel_p_i]
            fubar_pp = fubar_rows[i][fubar_pneg_i]
            slac_dn = slac_rows[i][slac_dn_i]
            slac_ds = slac_rows[i][slac_ds_i]
            slac_p = slac_rows[i][slac_p_i]
            out.write(
                f"{pos}\t{aa}\t{fmt_value(fel_dn)}\t{fmt_value(fel_ds)}\t"
                f"{fmt_value(fel_p)}\t{fmt_value(fubar_pp)}\t"
                f"{fmt_value(slac_dn)}\t{fmt_value(slac_ds)}\t{fmt_value(slac_p)}\n"
            )


def run_hyphy(cmd, log):
    log.info("CMD: %s", " ".join(cmd))
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    log.info(proc.stdout)
    if proc.returncode != 0:
        raise RuntimeError("HyPhy 运行失败")


def run(config_path, dataset_name, force=False, log_dir=None, tmp_dir=None):
    cfg = load_config(config_path)
    project_dir = cfg.get("project_dir")
    if not project_dir:
        raise ConfigError("配置缺少 project_dir")

    datasets = cfg.get("datasets", {})
    if dataset_name not in datasets:
        raise ConfigError(f"配置中未找到数据集: {dataset_name}")

    dataset = datasets[dataset_name]
    alignment = dataset.get("alignment")
    tree = dataset.get("tree")
    output_dir = dataset.get("output_dir")
    if not alignment or not tree or not output_dir:
        raise ConfigError(f"数据集 {dataset_name} 缺少 alignment/tree/output_dir")

    hyphy_bin = cfg.get("hyphy_bin", "hyphy")
    log_dir = log_dir or os.path.join(project_dir, "log")
    tmp_dir = tmp_dir or os.path.join(project_dir, "tmp")

    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)

    log_path = os.path.join(log_dir, f"hyphy_{dataset_name}.log")
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    log = logging.getLogger("hyphy")

    if not force and log_has_success(log_path):
        return 0

    ensure_file(alignment, "alignment")
    ensure_file(tree, "tree")

    if os.path.exists(output_dir) and not force:
        raise RuntimeError(f"输出目录已存在，请使用 --force 或清理: {output_dir}")

    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    tmp_run_dir = os.path.join(tmp_dir, f"hyphy_{dataset_name}_{run_id}")
    if os.path.exists(tmp_run_dir):
        shutil.rmtree(tmp_run_dir)
    os.makedirs(tmp_run_dir, exist_ok=True)

    fel_json = os.path.join(tmp_run_dir, "fel.json")
    fubar_json = os.path.join(tmp_run_dir, "fubar.json")
    slac_json = os.path.join(tmp_run_dir, "slac.json")
    tsv_tmp = os.path.join(tmp_run_dir, "summary.tsv")

    run_hyphy(
        [hyphy_bin, "FEL", "--alignment", alignment, "--tree", tree, "--branches", "All", "--output", fel_json],
        log,
    )
    run_hyphy(
        [hyphy_bin, "FUBAR", "--alignment", alignment, "--tree", tree, "--output", fubar_json],
        log,
    )
    run_hyphy(
        [hyphy_bin, "SLAC", "--alignment", alignment, "--tree", tree, "--branches", "All", "--output", slac_json],
        log,
    )

    ensure_file(fel_json, "FEL 输出")
    ensure_file(fubar_json, "FUBAR 输出")
    ensure_file(slac_json, "SLAC 输出")

    write_summary(alignment, fel_json, fubar_json, slac_json, tsv_tmp)
    ensure_file(tsv_tmp, "TSV 输出")

    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    shutil.move(fel_json, os.path.join(output_dir, "fel.json"))
    shutil.move(fubar_json, os.path.join(output_dir, "fubar.json"))
    shutil.move(slac_json, os.path.join(output_dir, "slac.json"))
    shutil.move(tsv_tmp, os.path.join(output_dir, "summary.tsv"))

    log.info("STATUS=SUCCESS dataset=%s", dataset_name)
    return 0


def list_datasets(config_path):
    cfg = load_config(config_path)
    datasets = cfg.get("datasets", {})
    for name in datasets.keys():
        print(name)


def main():
    parser = argparse.ArgumentParser(description="Run HyPhy FEL/FUBAR/SLAC for one dataset")
    parser.add_argument("--config", required=True, help="Path to JSON config")
    parser.add_argument("--dataset", help="Dataset name")
    parser.add_argument("--log-dir", help="Override log directory")
    parser.add_argument("--tmp-dir", help="Override tmp directory")
    parser.add_argument("--force", action="store_true", help="Overwrite existing outputs")
    parser.add_argument("--list-datasets", action="store_true", help="List dataset names")
    args = parser.parse_args()

    if args.list_datasets:
        list_datasets(args.config)
        return

    if not args.dataset:
        raise SystemExit("--dataset is required unless --list-datasets is set")

    run(
        config_path=args.config,
        dataset_name=args.dataset,
        force=args.force,
        log_dir=args.log_dir,
        tmp_dir=args.tmp_dir,
    )


if __name__ == "__main__":
    main()
