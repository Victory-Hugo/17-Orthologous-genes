#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LEVELS_ALL = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
LEVELS_HEATMAP = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


def sniff_sep(path: Path) -> str:
    first_line = path.read_text().splitlines()[0]
    if first_line.count("	") > first_line.count(","):
        return "	"
    return ","


def read_similarity_matrix(path: Path) -> pd.DataFrame:
    sep = sniff_sep(path)
    df = pd.read_csv(path, sep=sep, index_col=0)
    return df


def read_metadata(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    return df


def build_group_matrix(sim: pd.DataFrame, labels: pd.Series) -> pd.DataFrame:
    labels = labels.dropna()
    assemblies = labels.index.intersection(sim.index)
    labels = labels.loc[assemblies]
    sim = sim.loc[assemblies, assemblies]

    # Preserve appearance order of groups.
    seen = []
    for item in labels.tolist():
        if item not in seen:
            seen.append(item)

    index_map = {asm: i for i, asm in enumerate(sim.index)}
    group_indices = {}
    for group in seen:
        members = labels[labels == group].index
        group_indices[group] = np.array([index_map[a] for a in members], dtype=int)

    sim_vals = sim.values
    out = np.full((len(seen), len(seen)), np.nan, dtype=float)

    for i, gi in enumerate(seen):
        idx_i = group_indices[gi]
        for j, gj in enumerate(seen):
            idx_j = group_indices[gj]
            sub = sim_vals[np.ix_(idx_i, idx_j)]
            if i == j:
                if sub.size == 0:
                    val = np.nan
                else:
                    # Exclude diagonal self-comparisons for within-group.
                    mask = np.ones_like(sub, dtype=bool)
                    np.fill_diagonal(mask, False)
                    vals = sub[mask]
                    val = np.nan if vals.size == 0 else float(np.nanmean(vals))
            else:
                val = np.nan if sub.size == 0 else float(np.nanmean(sub))
            out[i, j] = val

    return pd.DataFrame(out, index=seen, columns=seen)


def save_heatmap(df: pd.DataFrame, out_path: Path, title: str) -> None:
    # Scale figure size based on number of groups.
    n = max(df.shape)
    base = 4
    size = max(base, min(18, 0.4 * n))
    fig, ax = plt.subplots(figsize=(size, size))
    im = ax.imshow(df.values, aspect="auto", cmap="viridis")
    ax.set_title(title)
    ax.set_xticks(range(df.shape[1]))
    ax.set_yticks(range(df.shape[0]))
    ax.set_xticklabels(df.columns, rotation=90, fontsize=8)
    ax.set_yticklabels(df.index, fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Mean similarity")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute taxonomic-level similarity relationships and heatmaps."
    )
    parser.add_argument(
        "--similarity",
        required=True,
        type=Path,
        help="Path to similarity_matrix.tsv",
    )
    parser.add_argument(
        "--metadata",
        required=True,
        type=Path,
        help="Path to metadata CSV with taxonomy columns",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory (default: similarity file directory)",
    )
    args = parser.parse_args()

    sim = read_similarity_matrix(args.similarity)
    meta = read_metadata(args.metadata)

    assembly_col = "Assembly" if "Assembly" in meta.columns else meta.columns[0]
    meta = meta.set_index(assembly_col, drop=False)

    # Normalize empty strings to NaN.
    for col in LEVELS_ALL:
        if col in meta.columns:
            meta[col] = meta[col].replace({"": np.nan, "NA": np.nan, "N/A": np.nan})

    common = sim.index.intersection(meta.index)
    if common.empty:
        raise SystemExit("No overlapping assemblies between similarity matrix and metadata.")

    out_dir = args.out_dir or args.similarity.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Similarity matrix assemblies: {len(sim)}")
    print(f"Metadata assemblies: {len(meta)}")
    print(f"Overlapping assemblies: {len(common)}")

    for level in LEVELS_ALL:
        if level not in meta.columns:
            print(f"Skip {level}: column missing in metadata")
            continue
        labels = meta.loc[common, level]
        labels = labels.dropna()
        if labels.empty:
            print(f"Skip {level}: no labels")
            continue

        group_matrix = build_group_matrix(sim.loc[common, common], labels)
        out_tsv = out_dir / f"taxonomy_relation_{level}.tsv"
        group_matrix.to_csv(out_tsv, sep="	", na_rep="NA")
        print(f"Wrote {out_tsv}")

        if level in LEVELS_HEATMAP:
            out_png = out_dir / f"taxonomy_heatmap_{level}.png"
            save_heatmap(group_matrix, out_png, f"{level} mean similarity")
            print(f"Wrote {out_png}")


if __name__ == "__main__":
    main()
