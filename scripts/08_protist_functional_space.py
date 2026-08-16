#!/usr/bin/env python3
"""Build functional trait space after evidence-based annotation is complete.

This is intentionally downstream of trait curation. It accepts a taxon x trait
matrix, computes Gower distances, performs PCoA, and writes coordinates for
cluster evaluation. Cluster labels are not biologically named here.
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def gower_distance(df):
    """Simple mixed-type Gower distance with pairwise missing-value handling."""
    n = len(df)
    out = np.zeros((n, n), dtype=float)
    numeric = {c: pd.api.types.is_numeric_dtype(df[c]) for c in df.columns}
    ranges = {}
    for c in df.columns:
        if numeric[c]:
            x = pd.to_numeric(df[c], errors="coerce")
            r = x.max() - x.min()
            ranges[c] = r if pd.notna(r) and r > 0 else 1.0

    for i in range(n):
        for j in range(i + 1, n):
            scores = []
            for c in df.columns:
                a, b = df.iloc[i][c], df.iloc[j][c]
                if pd.isna(a) or pd.isna(b) or str(a).strip() == "" or str(b).strip() == "":
                    continue
                if numeric[c]:
                    scores.append(abs(float(a) - float(b)) / ranges[c])
                else:
                    scores.append(0.0 if str(a) == str(b) else 1.0)
            d = float(np.mean(scores)) if scores else np.nan
            out[i, j] = out[j, i] = d
    return out


def pcoa(distance):
    n = distance.shape[0]
    d2 = distance ** 2
    h = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * h @ d2 @ h
    vals, vecs = np.linalg.eigh(b)
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    positive = vals > 0
    vals, vecs = vals[positive], vecs[:, positive]
    coords = vecs * np.sqrt(vals)
    return vals, coords


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--traits", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--taxon-col", default="taxon")
    args = ap.parse_args()

    x = pd.read_csv(args.traits)
    if args.taxon_col not in x:
        raise ValueError(f"Missing taxon column: {args.taxon_col}")
    x = x.set_index(args.taxon_col)

    # Remove provenance columns; retain only actual biological traits.
    provenance = [c for c in x.columns if c.startswith(("reference_", "evidence_", "confidence", "curator_"))]
    traits = x.drop(columns=provenance, errors="ignore")
    traits = traits.loc[traits.notna().sum(axis=1) > 0]

    d = gower_distance(traits)
    valid = ~np.isnan(d).any(axis=1)
    d = d[np.ix_(valid, valid)]
    taxa = traits.index[valid]
    vals, coords = pcoa(d)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(d, index=taxa, columns=taxa).to_csv(out / "protist_trait_gower_distance.csv")
    axes = [f"PCoA{i+1}" for i in range(coords.shape[1])]
    pd.DataFrame(coords, index=taxa, columns=axes).to_csv(out / "protist_trait_pcoa_coordinates.csv")
    pd.DataFrame({"axis": axes[:len(vals)], "eigenvalue": vals,
                  "variance_fraction": vals / vals.sum()}).to_csv(out / "protist_trait_pcoa_variance.csv", index=False)


if __name__ == "__main__":
    main()
