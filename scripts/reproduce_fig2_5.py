#!/usr/bin/env python3
"""
Reproduce SEMS Fig. 2.5: depth-resolved 16S taxonomic heatmap.

The script deliberately preserves the historical SILVA-style marine labels
(SAR11 clades, SAR86, SAR202, SAR406, Marine Groups II/III, etc.) because
these are ecologically interpretable and closest to the previous figure.

Input:
  data/processed/16S/otu_table_prok_clean.csv
  data/processed/16S/taxonomy_prok_clean.csv
  data/processed/16S/metadata_prok_clean.csv

Outputs:
  results/figures/Fig2_5_16S_taxa_heatmap.{svg,pdf}
  results/tables/Fig2_5_16S_taxa_heatmap_values.csv

Cells are mean relative abundance (%) for each year × season × ds2 group.
Taxa are selected by overall mean abundance and can be fixed explicitly
with TAXA_KEEP below for exact figure-to-figure reproducibility.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "processed" / "16S"
FIG = ROOT / "results" / "figures"
TAB = ROOT / "results" / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

otu = pd.read_csv(DATA / "otu_table_prok_clean.csv").set_index("OTUID")
tax = pd.read_csv(DATA / "taxonomy_prok_clean.csv").set_index("OTUID")
meta = pd.read_csv(DATA / "metadata_prok_clean.csv").set_index("sample-id")

samples = [s for s in otu.columns if s in meta.index]
X = otu[samples].T.astype(float)
P = X.div(X.sum(axis=1), axis=0).fillna(0)

DEPTH_ORDER = [
    "A_surface", "B_nearsurface", "C_DCM", "D_below_DCM",
    "F_300-600", "G_below_600", "H_near_bottom"
]
DEPTH_LABELS = {
    "A_surface": "surface",
    "B_nearsurface": "near surface",
    "C_DCM": "DCM",
    "D_below_DCM": "below DCM",
    "F_300-600": "300–600 m",
    "G_below_600": ">600 m",
    "H_near_bottom": "near bottom",
}

TAXA_KEEP = [
    "Nitrosopumilaceae",
    "Candidatus_Nitrosopelagicus",
    "Clade_Ia",
    "Clade_Ib",
    "Clade_II",
    "Clade_IV",
    "SAR11_clade",
    "SAR86_clade",
    "AEGEAN-169_marine_group",
    "SAR202_clade",
    "Marinimicrobia_(SAR406_clade)",
    "SAR324_clade(Marine_group_B)",
    "Marine_Group_II",
    "Marine_Group_III",
    "SAR116_clade",
    "Synechococcus_CC9902",
    "Prochlorococcus_MIT9313",
    "NS5_marine_group",
    "UBA10353_marine_group",
    "Rhodobacteraceae",
    "Candidatus_Actinomarina",
    "KI89A_clade",
    "Sva0996_marine_group",
    "SUP05_cluster",
    "Nitrospira",
]
TOP_N = 25

def useful_label(row):
    for rank in ["Genus", "Family", "Order", "Class", "Phylum"]:
        if rank not in row.index:
            continue
        v = row[rank]
        if pd.isna(v):
            continue
        v = str(v)
        if v in {"Unassigned", ""}:
            continue
        for prefix in ("g__", "f__", "o__", "c__", "p__"):
            if v.startswith(prefix):
                v = v[len(prefix):]
        if v:
            return v
    return "Unclassified"

labels = pd.Series({oid: useful_label(tax.loc[oid]) for oid in P.columns})
A = P.copy()
A.columns = labels.reindex(P.columns).values
A = A.T.groupby(level=0).sum().T

if TAXA_KEEP:
    taxa = [t for t in TAXA_KEEP if t in A.columns]
else:
    taxa = list(A.mean(axis=0).sort_values(ascending=False).head(TOP_N).index)

season = meta.loc[samples, "season"].astype(str).replace({
    "A_Winter": "Winter",
    "B_Summer": "Summer",
    "Summer": "Summer",
    "Winter": "Winter",
})
md = meta.loc[samples].copy()
md["season_display"] = season

records = []
col_labels = []
for year in sorted(md["year"].dropna().unique()):
    for seas in ["Winter", "Summer"]:
        for dep in DEPTH_ORDER:
            ss = md.index[
                (md["year"] == year) &
                (md["season_display"] == seas) &
                (md["ds2"] == dep)
            ].intersection(A.index)
            if len(ss) == 0:
                continue
            vals = A.loc[ss, taxa].mean(axis=0) * 100
            records.append(vals)
            col_labels.append((int(year), seas, dep, len(ss)))

H = pd.DataFrame(records, index=pd.MultiIndex.from_tuples(
    col_labels, names=["year", "season", "ds2", "n"]
)).T

out = H.stack(["year", "season", "ds2", "n"], future_stack=True).reset_index()
out.columns = ["taxon", "year", "season", "ds2", "n", "mean_relative_abundance_pct"]
out.to_csv(TAB / "Fig2_5_16S_taxa_heatmap_values.csv", index=False)

depth_num = {d:i for i,d in enumerate(DEPTH_ORDER)}
centers = {}
for t in H.index:
    weights = []
    positions = []
    for col in H.columns:
        val = H.loc[t, col]
        if np.isfinite(val) and val > 0:
            weights.append(val)
            positions.append(depth_num[col[2]])
    centers[t] = np.average(positions, weights=weights) if weights else np.inf
H = H.loc[sorted(H.index, key=lambda x: centers[x])]

M = H.values
positive = M[M > 0]
vmin = max(np.quantile(positive, .05), 0.005) if len(positive) else 0.005
vmax = np.quantile(positive, .995) if len(positive) else 10
display = np.where(M > 0, M, vmin / 2)

fig_w = max(12, 0.22 * H.shape[1])
fig_h = max(8, 0.32 * H.shape[0])
fig, ax = plt.subplots(figsize=(fig_w, fig_h))

im = ax.imshow(
    display,
    aspect="auto",
    interpolation="nearest",
    norm=LogNorm(vmin=vmin/2, vmax=vmax),
)

ax.set_yticks(np.arange(H.shape[0]))
ax.set_yticklabels(H.index, fontsize=8)

depth_short = {
    "A_surface":"S", "B_nearsurface":"NS", "C_DCM":"D",
    "D_below_DCM":"BD", "F_300-600":"M", "G_below_600":"B",
    "H_near_bottom":"NB",
}
ax.set_xticks(np.arange(H.shape[1]))
ax.set_xticklabels([depth_short[c[2]] for c in H.columns],
                   rotation=90, fontsize=6)
ax.set_xlabel("Depth category within year and season")
ax.set_ylabel("16S taxon")

last_year = None
last_season = None
for i, col in enumerate(H.columns):
    year, seas, dep, n = col
    if i and (year != last_year or seas != last_season):
        ax.axvline(i - .5, linewidth=.45)
    last_year, last_season = year, seas

groups = []
start = 0
cols = list(H.columns)
for i in range(1, len(cols)+1):
    if i == len(cols) or cols[i][:2] != cols[start][:2]:
        groups.append((start, i-1, cols[start][0], cols[start][1]))
        start = i
for a,b,year,seas in groups:
    ax.text((a+b)/2, -1.15, f"{year}\n{seas}",
            ha="center", va="bottom", fontsize=7)

cbar = fig.colorbar(im, ax=ax, fraction=.025, pad=.015)
cbar.set_label("Mean relative abundance (%)")

ax.set_title("16S vertical taxonomic succession", loc="left", fontweight="bold")
fig.tight_layout()

fig.savefig(FIG / "Fig2_5_16S_taxa_heatmap.svg", format="svg", bbox_inches="tight")
fig.savefig(FIG / "Fig2_5_16S_taxa_heatmap.pdf", format="pdf", bbox_inches="tight")
plt.close(fig)

print(f"Samples: {len(samples)}")
print(f"Taxa shown: {len(H)}")
print(f"Year-season-depth combinations: {H.shape[1]}")
