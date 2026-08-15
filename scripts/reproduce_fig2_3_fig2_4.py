#!/usr/bin/env python3
"""
Reproduce SEMS 16S Fig. 2.3 (Shannon + archaeal proportion)
and Fig. 2.4 (Bray-Curtis PCoA).

Inputs expected:
  data/processed/16S/otu_table_prok_clean.csv
  data/processed/16S/taxonomy_prok_clean.csv
  data/processed/16S/metadata_prok_clean.csv

Outputs:
  results/figures/Fig2_3_16S_diversity_archaea.{svg,pdf}
  results/figures/Fig2_4_16S_pcoa.{svg,pdf}
  results/tables/Fig2_3_values.csv
  results/tables/Fig2_4_PCoA_coordinates.csv
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

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

depth_order = [
    "surface", "nearsurface", "DCM", "below DCM",
    "300–600 m", "below 600 m", "near bottom"
]

def pretty_depth(x):
    s = str(x).lower().replace("_", " ").strip()
    if "near bottom" in s: return "near bottom"
    if "below 600" in s: return "below 600 m"
    if "300" in s and "600" in s: return "300–600 m"
    if "below dcm" in s: return "below DCM"
    if s == "dcm" or s.endswith(" dcm"): return "DCM"
    if "nearsurface" in s or "near surface" in s: return "nearsurface"
    if "surface" in s: return "surface"
    return str(x)

depth_col = "ds2" if "ds2" in meta.columns else "ds3"
depth = meta.loc[samples, depth_col].map(pretty_depth)

season_raw = meta.loc[samples, "season"].astype(str)
season = season_raw.replace({
    "A_Winter": "Winter",
    "B_Summer": "Summer",
    "Winter": "Winter",
    "Summer": "Summer",
})

# Relative abundance
P = X.div(X.sum(axis=1), axis=0).fillna(0)

# Shannon diversity
shannon = -(P.where(P > 0) * np.log(P.where(P > 0))).sum(axis=1)

# Archaeal fraction
arch_ids = tax.index[tax["Kingdom"].astype(str).str.contains("Archaea", case=False, na=False)]
arch_prop = X[arch_ids].sum(axis=1) / X.sum(axis=1)

summary = pd.DataFrame({
    "sample": samples,
    "depth": depth.values,
    "season": season.values,
    "Shannon_H": shannon.values,
    "Archaea_proportion": arch_prop.values,
})
summary.to_csv(TAB / "Fig2_3_values.csv", index=False)

# ---- Fig 2.3 ----
rng = np.random.default_rng(12)
fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.8), sharey=True)

for ax, var, title, xlab in [
    (axes[0], "Shannon_H", "A", "Shannon diversity (H′)"),
    (axes[1], "Archaea_proportion", "B", "Archaeal proportion"),
]:
    for yi, d in enumerate(depth_order):
        vals = summary.loc[summary.depth == d, var].dropna().values
        if not len(vals):
            continue
        jitter = rng.normal(0, 0.055, len(vals))
        ax.scatter(vals, yi + jitter, s=12, alpha=.35)
        q1, med, q3 = np.quantile(vals, [.25, .5, .75])
        lo, hi = np.quantile(vals, [.05, .95])
        ax.hlines(yi, lo, hi, linewidth=.8)
        ax.hlines(yi, q1, q3, linewidth=3.0)
        ax.scatter(med, yi, s=38, facecolors="white",
                   edgecolors="black", linewidth=.9, zorder=5)
        ax.text(hi, yi - .20, f"n={len(vals)}", fontsize=8.5, ha="left")

    ax.set_xlabel(xlab)
    ax.set_title(title, loc="left", fontweight="bold", fontsize=12)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)

axes[0].set_yticks(range(len(depth_order)))
axes[0].set_yticklabels(depth_order)
axes[0].invert_yaxis()
axes[1].tick_params(labelleft=False)

fig.tight_layout(w_pad=2.2)
fig.savefig(FIG / "Fig2_3_16S_diversity_archaea.svg", format="svg", bbox_inches="tight")
fig.savefig(FIG / "Fig2_3_16S_diversity_archaea.pdf", format="pdf", bbox_inches="tight")
plt.close(fig)

# ---- Fig 2.4 Bray-Curtis PCoA ----
D = squareform(pdist(P.values, metric="braycurtis"))
n = len(D)
J = np.eye(n) - np.ones((n, n)) / n
B = -0.5 * J @ (D ** 2) @ J
eigval, eigvec = np.linalg.eigh(B)
ix = np.argsort(eigval)[::-1]
eigval = eigval[ix]
eigvec = eigvec[:, ix]
pos = eigval > 0
coords = eigvec[:, pos] * np.sqrt(eigval[pos])
expl = eigval[pos] / eigval[pos].sum() * 100

pcoa = pd.DataFrame({
    "sample": samples,
    "PCo1": coords[:, 0],
    "PCo2": coords[:, 1],
    "depth": depth.values,
    "season": season.values,
})
pcoa.to_csv(TAB / "Fig2_4_PCoA_coordinates.csv", index=False)

markers = {
    "surface": "o",
    "nearsurface": "^",
    "DCM": "s",
    "below DCM": "P",
    "300–600 m": "v",
    "below 600 m": "X",
    "near bottom": "D",
}

seasons = ["Winter", "Summer"]
cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
season_color = {s: cycle[i] for i, s in enumerate(seasons)}

fig, ax = plt.subplots(figsize=(8.2, 6.2))
for d in depth_order:
    for seas in seasons:
        z = pcoa[(pcoa.depth == d) & (pcoa.season == seas)]
        if len(z):
            ax.scatter(
                z.PCo1, z.PCo2,
                marker=markers[d],
                s=34,
                alpha=.72,
                color=season_color[seas],
                edgecolors="none",
            )

ax.set_xlabel(f"PCo1 [{expl[0]:.1f}%]")
ax.set_ylabel(f"PCo2 [{expl[1]:.1f}%]")
ax.spines[["top", "right"]].set_visible(False)

depth_handles = [
    Line2D([], [], marker=markers[d], linestyle="None",
           markerfacecolor="0.35", markeredgecolor="0.35",
           markersize=7, label=d)
    for d in depth_order
]
leg1 = ax.legend(
    handles=depth_handles,
    title="Depth",
    frameon=False,
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
)
ax.add_artist(leg1)

season_handles = [
    Line2D([], [], marker="o", linestyle="None",
           markerfacecolor=season_color[s], markeredgecolor="none",
           markersize=7, label=s)
    for s in seasons
]
ax.legend(
    handles=season_handles,
    title="Season",
    frameon=False,
    bbox_to_anchor=(1.02, .42),
    loc="upper left",
)

fig.tight_layout()
fig.savefig(FIG / "Fig2_4_16S_pcoa.svg", format="svg", bbox_inches="tight")
fig.savefig(FIG / "Fig2_4_16S_pcoa.pdf", format="pdf", bbox_inches="tight")
plt.close(fig)

print(f"PCo1: {expl[0]:.1f}%")
print(f"PCo2: {expl[1]:.1f}%")
