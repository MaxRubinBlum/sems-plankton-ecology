#!/usr/bin/env python3
"""Rank microbial-eukaryote taxa for functional-trait annotation.

Inputs are the cleaned 18S taxonomy, feature table and metadata used by the
community analyses. The script collapses ASVs to the most reliable available
taxonomic label and ranks taxa for deep-sea trait curation.

No functional trait is inferred from depth distribution.
"""

from pathlib import Path
import argparse
import numpy as np
import pandas as pd

DEPTH_ORDER = [
    "A_surface", "B_nearsurface", "C_DCM", "D_below_DCM",
    "F_300-600", "G_below_600", "H_near_bottom",
]


def best_label(row):
    ranks = ["Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom"]
    suffixes = ("_X", "_XX", "_XXX", "_XXXX")
    for rank in ranks:
        value = str(row.get(rank, "")).strip()
        if value and value.lower() not in {"nan", "unclassified"} and not value.endswith(suffixes):
            return pd.Series([value, rank])
    for rank in ranks:
        value = str(row.get(rank, "")).strip()
        if value and value.lower() not in {"nan", "unclassified"}:
            return pd.Series([value, rank])
    return pd.Series(["Unclassified", "NA"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--taxonomy", required=True)
    ap.add_argument("--table", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--layer-field", default="ds2")
    args = ap.parse_args()

    tax = pd.read_csv(args.taxonomy).fillna("")
    otu = pd.read_csv(args.table)
    meta = pd.read_csv(args.metadata)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    sample_cols = [c for c in otu.columns if c in set(meta["sample-id"])]
    x = otu.set_index("OTUID")[sample_cols]
    tx = tax.set_index("OTUID").loc[x.index]

    microbial = ~tx["Kingdom"].astype(str).str.contains("Metazoa|Animalia", case=False, regex=True)
    x, tx = x.loc[microbial], tx.loc[microbial]

    labels = tx.apply(best_label, axis=1)
    tx[["best_taxon", "rank_used"]] = labels

    collapsed = x.copy()
    collapsed["best_taxon"] = tx["best_taxon"]
    collapsed = collapsed.groupby("best_taxon").sum()
    rel = collapsed.div(collapsed.sum(axis=0), axis=1) * 100

    smeta = meta.set_index("sample-id").loc[sample_cols]
    if args.layer_field not in smeta:
        raise ValueError(f"Missing metadata field: {args.layer_field}")

    layers = [d for d in DEPTH_ORDER if d in set(smeta[args.layer_field].astype(str))]
    layer_means = {}
    for layer in layers:
        cols = smeta.index[smeta[args.layer_field].astype(str) == layer]
        layer_means[layer] = rel[cols].mean(axis=1)
    layer_means = pd.DataFrame(layer_means)

    totals = x.sum(axis=1)
    rows = []
    for label, ids in tx.groupby("best_taxon").groups.items():
        rep = max(ids, key=lambda z: totals.get(z, 0))
        r = tx.loc[rep]
        rows.append({
            "taxon": label, "rank_used": r["rank_used"], "Kingdom": r["Kingdom"],
            "Phylum": r["Phylum"], "Class": r["Class"], "Order": r["Order"],
            "Family": r["Family"], "Genus": r["Genus"], "Species": r["Species"],
        })
    result = pd.DataFrame(rows).set_index("taxon").join(layer_means, how="left").fillna(0)

    upper = [d for d in ["A_surface", "B_nearsurface", "C_DCM"] if d in result]
    deep = [d for d in ["F_300-600", "G_below_600", "H_near_bottom"] if d in result]
    result["upper_mean_pct"] = result[upper].mean(axis=1)
    result["deep_mean_pct"] = result[deep].mean(axis=1)
    result["deep_enrichment_ratio"] = (result["deep_mean_pct"] + 1e-4) / (result["upper_mean_pct"] + 1e-4)
    result["dominant_layer"] = result[layers].idxmax(axis=1)

    result = result.sort_values(["deep_mean_pct", "deep_enrichment_ratio"], ascending=False)
    result.to_csv(out / "all_microeuk_taxa_depth_profiles.csv")

    targets = result[(result["deep_mean_pct"] >= 0.05) | (result["deep_enrichment_ratio"] >= 3)].copy()
    targets["recommended_priority"] = np.select(
        [targets["deep_mean_pct"] >= 1, targets["deep_mean_pct"] >= 0.2, targets["deep_mean_pct"] >= 0.05],
        ["A_critical", "B_high", "C_moderate"], default="D_low"
    )
    targets.to_csv(out / "deep_trait_targets.csv")


if __name__ == "__main__":
    main()
