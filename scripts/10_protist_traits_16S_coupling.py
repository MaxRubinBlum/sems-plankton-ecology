#!/usr/bin/env python3
"""Link sample-level protist functional traits to dominant 16S taxa.

The analysis is deliberately performed in two complementary ways:
1. residual correlations after removing categorical depth layer, season and year;
2. within-depth Spearman correlations combined across layers with Fisher-z weighting.

The second analysis is used to flag stable cross-domain associations that are
not driven solely by the vertical gradient. Associations are descriptive and
do not imply direct biological interactions.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, norm
from statsmodels.stats.multitest import multipletests

PROK = Path("data/processed")
TRAIT_FILE = Path("results/protist_traits/sample_level_core_functional_traits.csv")
OUT = Path("results/protist_traits/crossdomain")
OUT.mkdir(parents=True, exist_ok=True)

CORE_TRAITS = [
    "bona_fide_phototrophy", "constitutive_mixotrophy", "parasitism",
    "phagotrophy_radiolaria", "diplonemid_heterotrophy",
]
DEPTH_ORDER = [
    "A_surface", "B_nearsurface", "C_DCM", "D_below_DCM",
    "F_300-600", "G_below_600", "H_near_bottom",
]


def best_label(row):
    for rank in ["Genus", "Family", "Order", "Class", "Phylum"]:
        value = str(row[rank]).strip()
        if value and value.lower() not in {"nan", "unclassified"} and not value.endswith(("_X", "_XX", "_XXX", "_XXXX")):
            return f"{rank}:{value}"
    for rank in ["Genus", "Family", "Order", "Class", "Phylum"]:
        value = str(row[rank]).strip()
        if value and value.lower() not in {"nan", "unclassified"}:
            return f"{rank}:{value}"
    return "Unclassified"


def residualize(y, design):
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(y)
    out = np.full_like(y, np.nan, dtype=float)
    beta = np.linalg.lstsq(design[ok], y[ok], rcond=None)[0]
    out[ok] = y[ok] - design[ok] @ beta
    return out


def main():
    otu = pd.read_csv(PROK / "otu_table_prok_clean.csv")
    tax = pd.read_csv(PROK / "taxonomy_prok_clean.csv").fillna("")
    meta = pd.read_csv(PROK / "metadata_prok_clean.csv")
    traits = pd.read_csv(TRAIT_FILE)

    samples = sorted(set(meta["sample-id"]) & set(traits["sample-id"]))
    meta = meta.set_index("sample-id").loc[samples]
    traits = traits.set_index("sample-id").loc[samples]

    tax["taxon"] = tax.apply(best_label, axis=1)
    counts = otu.set_index("OTUID")[samples]
    counts["taxon"] = tax.set_index("OTUID").loc[counts.index, "taxon"]
    abundance = counts.groupby("taxon").sum()
    prop = abundance.div(abundance.sum(axis=0), axis=1)

    mean_ab = prop.mean(axis=1)
    prevalence = (prop > 0).mean(axis=1)
    keep = (mean_ab >= 0.001) & (prevalence >= 0.20)
    prop = prop.loc[keep]
    hellinger = np.sqrt(prop).T

    # Residual analysis.
    cov = traits[["ds2", "season", "year"]].astype(str)
    design = pd.get_dummies(cov, drop_first=True, dtype=float)
    design.insert(0, "Intercept", 1.0)
    M = design.values
    tax_res = {t: residualize(hellinger[t].values, M) for t in hellinger.columns}

    residual_rows = []
    for trait in CORE_TRAITS:
        y = traits[trait].astype(float).values
        yr = residualize(y, M)
        for tx in hellinger.columns:
            raw_r, raw_p = spearmanr(hellinger[tx].values, y, nan_policy="omit")
            res_r, res_p = spearmanr(tax_res[tx], yr, nan_policy="omit")
            residual_rows.append({
                "trait": trait, "prokaryote": tx,
                "rho_raw": raw_r, "p_raw": raw_p,
                "rho_residual": res_r, "p_residual": res_p,
                "mean_rel_abundance_pct": 100 * mean_ab[tx],
                "prevalence": prevalence[tx],
            })
    residual = pd.DataFrame(residual_rows)
    residual["FDR_residual"] = residual.groupby("trait")["p_residual"].transform(
        lambda p: multipletests(p, method="fdr_bh")[1]
    )
    residual.to_csv(OUT / "protist_traits_vs_16S_residual_correlations.csv", index=False)

    # Within-depth validation and fixed-effect Fisher-z combination.
    rows = []
    for trait in CORE_TRAITS:
        for tx in hellinger.columns:
            layer_stats, zs, weights = [], [], []
            for layer in DEPTH_ORDER:
                ids = [s for s in samples if traits.loc[s, "ds2"] == layer]
                if len(ids) < 8:
                    continue
                r, p = spearmanr(hellinger.loc[ids, tx], traits.loc[ids, trait])
                if not np.isfinite(r):
                    continue
                layer_stats.append((layer, len(ids), r, p))
                zs.append(np.arctanh(np.clip(r, -0.999999, 0.999999)))
                weights.append(max(len(ids) - 3, 1))
            if not zs:
                continue
            zbar = np.average(zs, weights=weights)
            pooled = np.tanh(zbar)
            se = 1 / np.sqrt(np.sum(weights))
            pmeta = 2 * norm.sf(abs(zbar / se))
            signs = [np.sign(x[2]) for x in layer_stats if abs(x[2]) >= 0.10]
            rows.append({
                "trait": trait, "prokaryote": tx,
                "pooled_within_depth_rho": pooled, "p_meta": pmeta,
                "n_layers": len(layer_stats),
                "positive_layers": sum(s > 0 for s in signs),
                "negative_layers": sum(s < 0 for s in signs),
                "median_layer_rho": np.median([x[2] for x in layer_stats]),
                "layer_details": "; ".join(
                    f"{d}:n={n},rho={r:.3f},p={p:.3g}" for d, n, r, p in layer_stats
                ),
            })
    stable = pd.DataFrame(rows)
    stable["FDR_meta"] = stable.groupby("trait")["p_meta"].transform(
        lambda p: multipletests(p, method="fdr_bh")[1]
    )
    stable["same_direction_layers"] = np.where(
        stable["pooled_within_depth_rho"] > 0,
        stable["positive_layers"], stable["negative_layers"]
    )
    stable["stability_fraction"] = stable["same_direction_layers"] / stable["n_layers"]
    stable.to_csv(OUT / "protist_traits_vs_16S_within_depth_meta.csv", index=False)

    strong = stable[
        (stable["FDR_meta"] < 0.01) &
        (stable["pooled_within_depth_rho"].abs() >= 0.20) &
        (stable["same_direction_layers"] >= 4) &
        (stable["stability_fraction"] >= 0.57)
    ].copy()
    strong.to_csv(OUT / "protist_traits_vs_16S_stable_links.csv", index=False)


if __name__ == "__main__":
    main()
