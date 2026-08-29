#!/usr/bin/env python3
"""Generate the final Figure 2 from the matched transition analysis outputs."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import gaussian_kde


def add_dcm_density(ax, dcm, seasons, xlim, colors):
    grid = np.linspace(xlim[0], xlim[1], 360)
    for season in seasons:
        vals = dcm.loc[dcm["season"] == season, "dcm_peak_depth_m"].dropna().to_numpy(float)
        if len(vals) < 3:
            continue
        dens = gaussian_kde(vals)(grid)
        dens /= dens.max()
        for x0, x1, a in zip(grid[:-1], grid[1:], dens[:-1]):
            if a > 0.015:
                ax.axvspan(x0, x1, facecolor=colors[season], alpha=float(0.16 * a), linewidth=0, zorder=0)


def clean_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def draw_scan(ax, state, matched, long16, dcm, letter, xlim, seasons, layer_colors, dcm_colors, show_legend=False):
    labels = {"16S": "16S composition", "18S": "18S composition", "traits": "18S trophic representation"}
    order = ["16S", "18S", "traits"]
    ax.set_xlim(*xlim)
    ax.set_ylim(0, 0.58)
    add_dcm_density(ax, dcm, seasons, xlim, dcm_colors)
    handles = []
    for layer in order:
        z = matched[(matched["state"] == state) & (matched["layer"] == layer)].sort_values("threshold_m")
        imax = int(np.argmax(z["partial_R2"].to_numpy()))
        h, = ax.plot(z["threshold_m"], z["partial_R2"], linewidth=2.4, marker="o", markevery=[imax], markersize=6, color=layer_colors[layer], label=labels[layer], zorder=3)
        handles.append(h)
    zl = long16[long16["state"] == state].sort_values("threshold_m")
    imax = int(np.argmax(zl["partial_R2"].to_numpy()))
    hlong, = ax.plot(zl["threshold_m"], zl["partial_R2"], linestyle="--", linewidth=1.6, alpha=0.7, color=layer_colors["16S"], label="16S long record (2018–2026)", zorder=2)
    ax.plot(zl["threshold_m"].iloc[imax], zl["partial_R2"].iloc[imax], marker="o", fillstyle="none", markersize=5.5, color=layer_colors["16S"], zorder=4)
    handles.append(hlong)
    ax.set_xlabel("Candidate transition depth (m)")
    ax.set_ylabel("Incremental variance explained (partial R²)")
    ax.text(0.01, 0.99, letter, transform=ax.transAxes, ha="left", va="top", fontweight="bold", fontsize=13)
    clean_axes(ax)
    if show_legend:
        leg1 = ax.legend(handles=handles, frameon=False, fontsize=8.2, loc="upper right")
        ax.add_artist(leg1)
        ax.legend(handles=[Patch(facecolor=dcm_colors["Winter"], alpha=0.16, label="Winter DCM"), Patch(facecolor=dcm_colors["Summer"], alpha=0.16, label="Summer DCM")], frameon=False, fontsize=8.0, loc="lower right")


def main(root):
    root = Path(root)
    out = root / "results" / "fig02"
    matched = pd.read_csv(out / "matched_boundary_scans.csv")
    long16 = pd.read_csv(out / "long_16S_2018_2026_scans.csv")
    bsum = pd.read_csv(out / "matched_bootstrap_summary.csv")
    dcm = pd.read_csv(out / "dcm_peaks_from_ctd.csv")
    cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    layer_colors = {"16S": cycle[0], "18S": cycle[1], "traits": cycle[2]}
    dcm_colors = {"Winter": cycle[0], "Summer": cycle[1]}
    fig = plt.figure(figsize=(13.2, 9.0))
    gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.34)
    axa = fig.add_subplot(gs[0, 0]); axb = fig.add_subplot(gs[0, 1]); axc = fig.add_subplot(gs[1, 0]); axd = fig.add_subplot(gs[1, 1])
    draw_scan(axa, "pooled", matched, long16, dcm, "a", (0, 420), ["Winter", "Summer"], layer_colors, dcm_colors, True)
    draw_scan(axc, "A_Winter", matched, long16, dcm, "c", (0, 370), ["Winter"], layer_colors, dcm_colors)
    draw_scan(axd, "Summer", matched, long16, dcm, "d", (0, 420), ["Summer"], layer_colors, dcm_colors)
    labels = {"16S": "16S composition", "18S": "18S composition", "traits": "18S trophic representation"}
    ypos = {"16S": 2, "18S": 1, "traits": 0}
    for layer in ["16S", "18S", "traits"]:
        r = bsum[bsum["layer"] == layer].iloc[0]
        med, lo, hi = float(r["median_optimum_m"]), float(r["ci95_low_m"]), float(r["ci95_high_m"])
        axb.errorbar(med, ypos[layer], xerr=np.array([[med - lo], [hi - med]]), fmt="o", markersize=7, linewidth=2.2, capsize=4, color=layer_colors[layer])
    axb.set_yticks([2, 1, 0], [labels["16S"], labels["18S"], labels["traits"]])
    axb.set_xlim(0, 275); axb.set_ylim(-0.75, 2.65)
    axb.set_xlabel("Bootstrapped transition depth (m)")
    axb.text(0.01, 0.99, "b", transform=axb.transAxes, ha="left", va="top", fontweight="bold", fontsize=13)
    clean_axes(axb)
    fig.subplots_adjust(left=0.10, right=0.98, top=0.98, bottom=0.08)
    fig.savefig(out / "Figure2_final_clean.svg", bbox_inches="tight")
    fig.savefig(out / "Figure2_final_clean.pdf", bbox_inches="tight")
    fig.savefig(out / "Figure2_final_clean.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=Path(__file__).resolve().parents[2])
    args = ap.parse_args()
    main(args.root)
