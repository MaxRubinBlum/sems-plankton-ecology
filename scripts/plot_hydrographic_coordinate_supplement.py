from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results" / "hydrography"
OUT = ROOT / "figures" / "supplementary"
OUT.mkdir(parents=True, exist_ok=True)

# This plotting script expects the mapped bootstrap table produced by the hydrographic-coordinate analysis.
boot = pd.read_csv(RES / "boundary_bootstraps_mapped_to_hydrography.csv")
align = pd.read_csv(RES / "hydrographic_coordinate_seasonal_alignment.csv")
boot["season_display"] = boot["season"].replace({"A_Winter": "Winter", "Summer": "Summer"})

labels = {"16S": "16S composition", "18S": "18S composition", "18S_function": "18S trophic structure"}
order = ["16S", "18S", "18S_function"]

fig = plt.figure(figsize=(8.5, 12))
gs = fig.add_gridspec(3, 1, hspace=0.50)

def box_panel(ax, col, ylabel, letter, invert=False):
    positions = [1, 2, 4, 5, 7, 8]
    data = []
    for ds in order:
        for season in ["Winter", "Summer"]:
            data.append(boot[(boot.dataset == ds) & (boot.season_display == season)][col].dropna().values)
    ax.boxplot(data, positions=positions, widths=0.62, showfliers=False,
               medianprops={"linewidth": 1.7}, boxprops={"linewidth": 1.1},
               whiskerprops={"linewidth": 1.0}, capprops={"linewidth": 1.0})
    ax.set_xticks(positions)
    ax.set_xticklabels(["Winter", "Summer"] * 3)
    ax.set_ylabel(ylabel)
    ax.axvline(3, lw=0.6, alpha=0.25)
    ax.axvline(6, lw=0.6, alpha=0.25)
    for center, ds in zip([1.5, 4.5, 7.5], order):
        ax.text(center, -0.17, labels[ds], ha="center", va="top", transform=ax.get_xaxis_transform(), fontsize=10)
    ax.text(-0.07, 1.02, letter, transform=ax.transAxes, fontsize=13, fontweight="bold")
    if invert:
        ax.invert_yaxis()

box_panel(fig.add_subplot(gs[0]), "depth_m", "Bootstrap breakpoint depth (m)", "a", True)
box_panel(fig.add_subplot(gs[1]), "salinity", "Salinity at bootstrap breakpoint", "b")

ax = fig.add_subplot(gs[2])
coords = ["depth_m", "salinity", "sigma_theta", "oxygen_umol_kg", "hydro_PC1"]
coord_labels = ["Depth", "Salinity", "σθ", "O₂", "Hydro-PC1"]
x = np.arange(len(coords))
for offset, ds in zip([-0.04, 0, 0.04], order):
    g = align[align.dataset == ds].set_index("coordinate")
    vals = np.array([g.loc[c, "standardized_seasonal_separation"] for c in coords], float)
    line, = ax.plot(x + offset, vals, marker="o", lw=1.6, ms=4)
    ax.text(4.18, vals[-1], labels[ds], va="center", ha="left", fontsize=9, color=line.get_color())
ax.set_xticks(x)
ax.set_xticklabels(coord_labels)
ax.set_ylabel("Standardized winter–summer separation\n(smaller = stronger seasonal alignment)")
ax.set_xlim(-0.25, 5.35)
ax.set_ylim(-0.08, 3.0)
ax.text(-0.07, 1.02, "c", transform=ax.transAxes, fontsize=13, fontweight="bold")

fig.subplots_adjust(left=0.14, right=0.88, top=0.98, bottom=0.06)
for ext in ["svg", "png", "pdf"]:
    fig.savefig(OUT / f"Supplementary_Figure_hydrographic_coordinate_boundaries.{ext}", dpi=450, bbox_inches="tight")
plt.close(fig)
