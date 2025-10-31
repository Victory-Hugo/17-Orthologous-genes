import pandas as pd
import re

# Load
file_path = "/home/luolintao/0_Github/17-Orthologous-genes/0-GTDB数据库/example/conf/gtdb_taxonomy.tsv"
out_path = "/home/luolintao/0_Github/17-Orthologous-genes/0-GTDB数据库/example/conf/整理_merge.tsv"


df = pd.read_csv(file_path, sep="\t", header=None, names=["accession", "taxonomy"])

# Split taxonomy into ranks (GTDB order: d, p, c, o, f, g, s)
split_cols = df["taxonomy"].str.split(";", n=6, expand=True)

# Ensure we always have 7 columns
while split_cols.shape[1] < 7:
    split_cols[split_cols.shape[1]] = None

split_cols.columns = ["d", "p", "c", "o", "f", "g", "s"]

# Clean: strip whitespace and drop the rank prefixes like "d__", "p__"
def strip_prefix(series):
    # remove leading whitespace, rank prefix, and surrounding spaces
    return (
        series.fillna("")
              .str.strip()
              .str.replace(r"^[a-z]__","", regex=True)
              .replace("", pd.NA)
    )

clean = split_cols.apply(strip_prefix)

# Assemble final tidy table
tidy = pd.concat([df["accession"], clean], axis=1)
tidy = tidy.rename(columns={
    "d": "domain",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species"
})

# Save
tidy.to_csv(out_path, sep="\t", index=False)
