#!/usr/bin/env python3
"""Summarize conservative evidence-backed protist traits across depth.

This script complements the full Ramond/Gower workflow with transparent trait
signals that are independently interpretable before clustering. It deliberately
keeps bona fide phototrophy separate from constitutive mixotrophy and does not
infer function from depth or co-occurrence.

Expected inputs (default paths can be edited or passed as arguments in later
workflow wrappers):
  data/processed/taxonomy_euk_clean.csv
  data/processed/otu_table_euk_clean.csv
  data/processed/metadata_euk_clean.csv

Outputs are written to results/protist_traits/.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests

BASE = Path("data/processed")
OUT = Path("results/protist_traits")
OUT.mkdir(parents=True, exist_ok=True)

TAX = BASE / "taxonomy_euk_clean.csv"
OTU = BASE / "otu_table_euk_clean.csv"
META = BASE / "metadata_euk_clean.csv"

DEPTH_ORDER = [
    "A_surface", "B_nearsurface", "C_DCM", "D_below_DCM",
    "F_300-600", "G_below_600", "H_near_bottom",
]


def trait_flags(row):
    text = " | ".join(str(row.get(c, "")) for c in
                      ["Phylum", "Class", "Order", "Family", "Genus", "Species"]).lower()
    flags = {k: False for k in [
        "bona_fide_phototrophy", "constitutive_mixotrophy", "parasitism",
        "phagotrophy_radiolaria", "diplonemid_heterotrophy",
    ]}

    # Conservative bona fide phototroph set: lineages in which photosynthesis
    # is the defining nutritional mode; ambiguous broad dinoflagellate and
    # haptophyte assignments are excluded here.
    photo_tokens = [
        "mamiellophyceae", "micromonas", "bathycoccus", "ostreococcus",
        "prasinococcales", "prasinoderma", "chloropicophyceae",
        "trebouxiophyceae", "pedinophyceae", "bacillariophyceae",
        "mediophyceae", "coscinodiscophyceae", "bolidophyceae",
        "pelagophyceae", "eustigmatophyceae", "pinguiophyceae",
    ]
    if any(t in text for t in photo_tokens):
        flags["bona_fide_phototrophy"] = True

    # Named lineages with strong evidence for constitutive phago-mixotrophy.
    mixo_tokens = [
        "chrysochromulina", "haptolina", "prymnesium", "pavlova",
        "dinophysis", "tripos", "ceratium", "karlodinium",
    ]
    if any(t in text for t in mixo_tokens):
        flags["constitutive_mixotrophy"] = True

    if any(t in text for t in [
        "dino-group-i", "dino-group-ii", "syndinial", "hematodinium", "thalassomyces"
    ]):
        flags["parasitism"] = True

    if any(t in text for t in [
        "radiolaria", "rad-a", "rad-b", "rad-c", "collodaria",
        "sphaerozo", "acantha", "astrosphaer", "polycyst",
    ]):
        flags["phagotrophy_radiolaria"] = True

    if any(t in text for t in ["eupelagon", "diplonem", "dspd"]):
        flags["diplonemid_heterotrophy"] = True

    return pd.Series(flags)


def main():
    tax = pd.read_csv(TAX).fillna("")
    otu = pd.read_csv(OTU)
    meta = pd.read_csv(META)
    sample_cols = [c for c in otu.columns if c in set(meta["sample-id"])]

    counts = otu.set_index("OTUID")[sample_cols]
    taxonomy = tax.set_index("OTUID").loc[counts.index]
    keep = ~taxonomy["Kingdom"].astype(str).str.contains(
        "Metazoa|Animalia", case=False, regex=True
    )
    counts, taxonomy = counts.loc[keep], taxonomy.loc[keep]
    flags = taxonomy.apply(trait_flags, axis=1)

    rel = counts.div(counts.sum(axis=0), axis=1) * 100
    profile = pd.DataFrame(index=sample_cols)
    for trait in flags.columns:
        profile[trait] = rel.loc[flags.index[flags[trait]]].sum(axis=0)
    profile = profile.join(meta.set_index("sample-id")[[
        "ds2", "depth", "station", "year", "season", "cruise"
    ]])

    summary_rows, tests = [], []
    for trait in flags.columns:
        groups = []
        for layer in DEPTH_ORDER:
            values = profile.loc[profile["ds2"] == layer, trait].dropna()
            groups.append(values.values)
            summary_rows.append({
                "trait": trait, "depth_layer": layer, "n": len(values),
                "mean_pct": values.mean(), "median_pct": values.median(),
                "q25": values.quantile(0.25), "q75": values.quantile(0.75),
            })
        groups = [g for g in groups if len(g)]
        stat, p = kruskal(*groups)
        tests.append({"trait": trait, "kruskal_H": stat, "p": p})

    summary = pd.DataFrame(summary_rows)
    tests = pd.DataFrame(tests)
    tests["FDR"] = multipletests(tests["p"], method="fdr_bh")[1]

    # Contribution of bona fide phototroph taxa by layer.
    photo_ids = flags.index[flags["bona_fide_phototrophy"]]
    photo_tax = taxonomy.loc[photo_ids].copy()
    def best_label(r):
        for col in ["Genus", "Family", "Order", "Class"]:
            value = str(r[col]).strip()
            if value and value.lower() != "nan":
                return value
        return "Unclassified"
    photo_tax["label"] = photo_tax.apply(best_label, axis=1)
    photo_rel = rel.loc[photo_ids].copy()
    photo_rel["label"] = photo_tax["label"]
    photo_rel = photo_rel.groupby("label").sum()
    photo_depth = pd.DataFrame({
        layer: photo_rel[profile.index[profile["ds2"] == layer]].mean(axis=1)
        for layer in DEPTH_ORDER
    })
    photo_depth["overall_mean"] = photo_depth.mean(axis=1)
    photo_depth = photo_depth.sort_values("overall_mean", ascending=False)

    profile.reset_index(names="sample-id").to_csv(
        OUT / "sample_level_core_functional_traits.csv", index=False)
    summary.to_csv(OUT / "core_functional_trait_summary_by_depth.csv", index=False)
    tests.to_csv(OUT / "core_functional_trait_depth_tests.csv", index=False)
    photo_depth.to_csv(OUT / "bona_fide_phototroph_taxa_by_depth.csv")


if __name__ == "__main__":
    main()
