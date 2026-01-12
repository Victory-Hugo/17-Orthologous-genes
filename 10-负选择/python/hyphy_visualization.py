#!/usr/bin/env python3
import argparse
import csv
import json
import logging
import os
import shutil
from datetime import datetime

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


class ConfigError(RuntimeError):
    pass


def load_config(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def ensure_file(path, desc):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise FileNotFoundError(f"{desc} 不存在或为空: {path}")


def log_has_success(log_path):
    if not os.path.isfile(log_path):
        return False
    with open(log_path, "r", encoding="utf-8") as f:
        for line in f:
            if "STATUS=SUCCESS" in line:
                return True
    return False


def parse_float(value):
    if value is None:
        return None
    value = value.strip()
    if value == "" or value.lower() == "none":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def load_summary(summary_path):
    rows = []
    with open(summary_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def get_column(rows, key):
    values = []
    for row in rows:
        val = parse_float(row.get(key, ""))
        if val is not None:
            values.append(val)
    return values


def get_xy_pairs(rows, key_x, key_y):
    xs = []
    ys = []
    for row in rows:
        x_val = parse_float(row.get(key_x, ""))
        y_val = parse_float(row.get(key_y, ""))
        if x_val is None or y_val is None:
            continue
        xs.append(x_val)
        ys.append(y_val)
    return xs, ys


def plot_hist(values, title, xlabel, out_path):
    if not values:
        return False
    plt.figure(figsize=(6, 4))
    plt.hist(values, bins=50, color="#4C72B0", edgecolor="#1F2D3D")
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()
    return True


def plot_scatter(x, y, title, xlabel, ylabel, out_path):
    if not x or not y:
        return False
    plt.figure(figsize=(5, 5))
    plt.scatter(x, y, s=8, alpha=0.6, color="#2A9D8F")
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()
    return True


def write_plots(summary_path, output_dir, log):
    rows = load_summary(summary_path)

    fel_p = get_column(rows, "FEL_pvalue")
    fubar_pp = get_column(rows, "FUBAR_posterior_negative")
    slac_p = get_column(rows, "SLAC_pvalue")

    fel_dn, fel_ds = get_xy_pairs(rows, "FEL_dN", "FEL_dS")
    slac_dn, slac_ds = get_xy_pairs(rows, "SLAC_dN", "SLAC_dS")

    plots = [
        (plot_hist, fel_p, "FEL p-value", "FEL p-value", "fel_pvalue_hist.png"),
        (plot_hist, fubar_pp, "FUBAR posterior negative", "Posterior probability", "fubar_posterior_negative_hist.png"),
        (plot_hist, slac_p, "SLAC P[dN/dS < 1]", "SLAC P[dN/dS < 1]", "slac_pvalue_hist.png"),
    ]

    for func, values, title, xlabel, name in plots:
        out_path = os.path.join(output_dir, name)
        ok = func(values, title, xlabel, out_path)
        log.info("PLOT %s created=%s", name, ok)

    scatter_pairs = [
        (fel_dn, fel_ds, "FEL dN vs dS", "FEL dN", "FEL dS", "fel_dn_vs_ds.png"),
        (slac_dn, slac_ds, "SLAC dN vs dS", "SLAC dN", "SLAC dS", "slac_dn_vs_ds.png"),
    ]

    for x, y, title, xlabel, ylabel, name in scatter_pairs:
        out_path = os.path.join(output_dir, name)
        ok = plot_scatter(x, y, title, xlabel, ylabel, out_path)
        log.info("PLOT %s created=%s", name, ok)


def run(config_path, dataset_name, force=False, log_dir=None, tmp_dir=None):
    cfg = load_config(config_path)
    project_dir = cfg.get("project_dir")
    if not project_dir:
        raise ConfigError("配置缺少 project_dir")

    datasets = cfg.get("datasets", {})
    if dataset_name not in datasets:
        raise ConfigError(f"配置中未找到数据集: {dataset_name}")

    dataset = datasets[dataset_name]
    summary_tsv = dataset.get("summary_tsv")
    output_dir = dataset.get("output_dir")
    if not summary_tsv or not output_dir:
        raise ConfigError(f"数据集 {dataset_name} 缺少 summary_tsv/output_dir")

    log_dir = log_dir or os.path.join(project_dir, "log")
    tmp_dir = tmp_dir or os.path.join(project_dir, "tmp")

    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)

    log_path = os.path.join(log_dir, f"visualize_{dataset_name}.log")
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    log = logging.getLogger("visualize")

    if not force and log_has_success(log_path):
        return 0

    ensure_file(summary_tsv, "summary.tsv")

    if os.path.exists(output_dir) and not force:
        raise RuntimeError(f"输出目录已存在，请使用 --force 或清理: {output_dir}")

    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    tmp_run_dir = os.path.join(tmp_dir, f"visualize_{dataset_name}_{run_id}")
    if os.path.exists(tmp_run_dir):
        shutil.rmtree(tmp_run_dir)
    os.makedirs(tmp_run_dir, exist_ok=True)

    write_plots(summary_tsv, tmp_run_dir, log)

    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    for name in os.listdir(tmp_run_dir):
        src = os.path.join(tmp_run_dir, name)
        dst = os.path.join(output_dir, name)
        shutil.move(src, dst)

    log.info("STATUS=SUCCESS dataset=%s", dataset_name)
    return 0


def list_datasets(config_path):
    cfg = load_config(config_path)
    datasets = cfg.get("datasets", {})
    for name in datasets.keys():
        print(name)


def main():
    parser = argparse.ArgumentParser(description="Visualize HyPhy summary.tsv for one dataset")
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
