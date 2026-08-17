"""Reproduce the analytical data underlying manuscript Figure 3.

Run from repository root:
    python reproducibility/fig03/make_figure.py

The script writes deterministic derived tables to results/fig03/. Plot styling can
be extended here, but taxon selection/scaling must not be altered downstream.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
OUT = ROOT / "results" / "fig03"
OUT.mkdir(parents=True, exist_ok=True)

DEPTH_ORDER = ["A_surface","B_nearsurface","C_DCM","D_below_DCM","F_300-600","G_below_600","H_near_bottom"]
SEASONS = ["A_Winter", "Summer"]
LAYER_DEPTH = np.array([1, 50, 115, 190, 375, 1000, 1450], dtype=float)
TRAITS = ["bona_fide_phototrophy","constitutive_mixotrophy","parasitism","phagotrophy_radiolaria","diplonemid_heterotrophy"]


def clean_tax(x):
    x = str(x or "")
    if "__" in x:
        x = x.split("__", 1)[1]
    return (x.replace("_", " ").strip()
             .replace(" clade(Marine group B)", "")
             .replace(" X", " (unresolved)")
             .replace("c  ", ""))


def aggregate(otu, tax, meta):
    mm = meta[meta.ds2.isin(DEPTH_ORDER) & meta.season.isin(SEASONS)].dropna(subset=["depth"]).set_index("sample-id")
    ids = [s for s in mm.index if s in otu.columns]
    x = otu[ids].T.astype(float)
    x = x.div(x.sum(axis=1).replace(0, np.nan), axis=0) * 100
    labels = []
    for oid in otu.index:
        r = tax.loc[oid]
        vals = [clean_tax(r[c]) if c in r.index and str(r[c]).strip() else "" for c in ["Order", "Class", "Family"]]
        labels.append(vals[0] or vals[1] or vals[2] or "Unclassified")
    a = x.T.copy()
    a["taxon"] = pd.Series(labels, index=otu.index).reindex(a.index).values
    return a.groupby("taxon").sum().T, mm.loc[ids]


def select_and_split(g, mm, n=14, min_peak=0.20):
    pooled = g.groupby(mm.ds2).median().reindex(DEPTH_ORDER).T.fillna(0)
    eligible = pooled.loc[pooled.max(axis=1) >= min_peak]
    score = (eligible.max(axis=1) - eligible.min(axis=1)) * np.log1p(eligible.max(axis=1) * 5)
    p = pooled.loc[score.nlargest(n).index]
    centroid = (p.values * LAYER_DEPTH).sum(axis=1) / (p.values.sum(axis=1) + 1e-9)
    p = p.iloc[np.argsort(centroid)]
    seasonal = {}
    for season in SEASONS:
        ms = mm[mm.season == season]
        seasonal[season] = g.loc[ms.index].groupby(ms.ds2).median().reindex(DEPTH_ORDER).T.fillna(0).reindex(p.index)
    w, s = seasonal["A_Winter"], seasonal["Summer"]
    both = np.c_[w.values, s.values]
    mu = both.mean(axis=1, keepdims=True)
    sd = both.std(axis=1, ddof=1, keepdims=True)
    sd[sd == 0] = 1
    zw = pd.DataFrame((w.values-mu)/sd, index=p.index, columns=DEPTH_ORDER)
    zs = pd.DataFrame((s.values-mu)/sd, index=p.index, columns=DEPTH_ORDER)
    return p, w, s, zw, zs


def trait_summaries(df):
    grid = np.geomspace(8, 1700, 120)
    rows = []
    for trait in TRAITS:
        for season in SEASONS:
            z = df[(df.season == season)].dropna(subset=[trait, "depth"])
            dep, val = z.depth.to_numpy(float), z[trait].to_numpy(float)
            for focal in grid:
                v = val[np.abs(dep-focal) <= max(35, 0.22*focal)]
                if len(v) >= 6:
                    rows.append([trait, season, focal, len(v), np.median(v), np.quantile(v,.25), np.quantile(v,.75)])
    return pd.DataFrame(rows, columns=["trait","season","depth_m","n","median","q25","q75"])


def run(marker):
    otu = pd.read_csv(DATA / ("otu_table_prok_clean.csv" if marker == "16S" else "otu_table_euk_clean.csv")).set_index("OTUID")
    tax = pd.read_csv(DATA / ("taxonomy_prok_clean.csv" if marker == "16S" else "taxonomy_euk_clean.csv")).set_index("OTUID").fillna("")
    meta = pd.read_csv(DATA / f"{marker}_samples_CTD_chemistry.csv")
    g, mm = aggregate(otu, tax, meta)
    p,w,s,zw,zs = select_and_split(g, mm)
    p.to_csv(OUT / f"{marker}_selected_taxa_pooled.csv")
    pd.concat({"winter":w,"summer":s}, axis=1).to_csv(OUT / f"{marker}_winter_summer_medians.csv")
    pd.concat({"winter":zw,"summer":zs}, axis=1).to_csv(OUT / f"{marker}_winter_summer_zscores.csv")


if __name__ == "__main__":
    run("16S")
    run("18S")
    traits = pd.read_csv(DATA / "sample_level_core_functional_traits.csv")
    trait_summaries(traits).to_csv(OUT / "trait_depth_summaries.csv", index=False)
