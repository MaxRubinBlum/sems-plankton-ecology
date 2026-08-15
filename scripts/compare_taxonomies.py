#!/usr/bin/env python3
"""Compare original SILVA and GTDB r226 taxonomic assignments by ASV ID."""

import re
from pathlib import Path
import numpy as np
import pandas as pd

RAW = Path("data/raw/16S")
OUT = Path("results/taxonomy_comparison")
OUT.mkdir(parents=True, exist_ok=True)

silva = pd.read_csv(RAW / "taxonomy_SILVA.csv")
gtdb = pd.read_csv(RAW / "taxonomy_GTDB226.tsv", sep="\t")

ranks = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
prefixes = ["d__", "p__", "c__", "o__", "f__", "g__", "s__"]


def parse_gtdb(taxon):
    out = {rank: np.nan for rank in ranks}
    if pd.isna(taxon):
        return out
    for part in str(taxon).split(";"):
        for rank, prefix in zip(ranks, prefixes):
            if part.startswith(prefix):
                out[rank] = part
                break
    return out


def strip_prefix(value):
    if pd.isna(value):
        return np.nan
    return re.sub(r"^[dkpcofgs]__", "", str(value))


parsed = pd.DataFrame([parse_gtdb(x) for x in gtdb["Taxon"]])
gtdb2 = pd.concat(
    [gtdb.rename(columns={"Feature ID": "OTUID"}).reset_index(drop=True), parsed],
    axis=1,
)

merged = silva.merge(gtdb2, on="OTUID", suffixes=("_SILVA", "_GTDB"), how="outer", indicator=True)
merged.to_csv(OUT / "SILVA_GTDB226_joined_assignments.csv", index=False)

coverage = []
for rank, prefix in zip(ranks, prefixes):
    s = merged[f"{rank}_SILVA"]
    g = merged[f"{rank}_GTDB"]
    s_assigned = s.notna() & ~s.astype(str).isin([prefix, "Unassigned", "nan"])
    g_assigned = g.notna() & ~g.astype(str).isin([prefix, "Unassigned", "nan"])
    coverage.append(
        [rank, int(s_assigned.sum()), float(s_assigned.mean()), int(g_assigned.sum()), float(g_assigned.mean())]
    )

coverage = pd.DataFrame(
    coverage,
    columns=["rank", "SILVA_assigned", "SILVA_fraction", "GTDB_assigned", "GTDB_fraction"],
)
coverage.to_csv(OUT / "rank_coverage.csv", index=False)

agreement = []
for rank in ranks:
    s = merged[f"{rank}_SILVA"].map(strip_prefix)
    g = merged[f"{rank}_GTDB"].map(strip_prefix)
    both = s.notna() & g.notna() & (s != "") & (g != "")
    exact = s[both] == g[both]
    agreement.append(
        [rank, int(both.sum()), int(exact.sum()), float(exact.mean()) if both.sum() else np.nan]
    )

agreement = pd.DataFrame(
    agreement,
    columns=["rank", "both_assigned", "exact_same_name", "exact_agreement"],
)
agreement.to_csv(OUT / "exact_name_agreement.csv", index=False)

print("ASVs in SILVA:", len(silva))
print("ASVs in GTDB226:", len(gtdb))
print("Join status:")
print(merged["_merge"].value_counts())
print("\nRank coverage:\n", coverage.to_string(index=False))
print("\nExact name agreement:\n", agreement.to_string(index=False))
