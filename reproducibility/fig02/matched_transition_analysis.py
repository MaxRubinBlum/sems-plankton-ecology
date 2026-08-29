#!/usr/bin/env python3
"""Authoritative Figure 2 transition analysis.

The primary comparison uses the exact sample universe shared by the 16S and
18S 99%-OTU tables and the five-component 18S trophic representation. Pooled
scans control station and monitoring state; state-specific scans control
station. Candidate-depth scans require at least 12 independent sampling profiles
spanning both sides of the split. The joint bootstrap retains the original,
more conservative 15-profile support criterion; sensitivity to minimum scan
support of 10, 12 and 15 profiles is exported explicitly.

Complete sampling profiles are resampled jointly across all three biological
layers, allowing direct bootstrap inference on pairwise transition-depth
differences. Scan-wide significance is assessed by profile-restricted depth
permutations. The longer 2018-2026 16S record is analysed separately with the
same nuisance structure and is used only as contextual validation in Figure 2.
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd

THRESHOLDS = np.arange(10, 421, 10, dtype=float)
MAX_DEPTH = 650.0
MIN_SPANNING = 12
BOOT_MIN_SPANNING = 15
SUPPORT_SENSITIVITY = (10, 12, 15)
N_BOOT = 500
N_SEASON_BOOT = 300
N_PERM = 499
SEED = 20260829
TRAITS = [
    "bona_fide_phototrophy",
    "constitutive_mixotrophy",
    "parasitism",
    "phagotrophy_radiolaria",
    "diplonemid_heterotrophy",
]


def pick(*paths):
    for p in paths:
        if p.exists():
            return p
    raise FileNotFoundError("None of these inputs exists: " + ", ".join(map(str, paths)))


def profile_id(meta):
    cruise = meta["cruise"].astype("string")
    fallback = (
        meta["year"].astype("Int64").astype(str) + "|" +
        meta["month"].astype(str) + "|" +
        meta["station"].astype(str)
    )
    return np.where(
        cruise.notna() & (cruise.str.len().fillna(0) > 0),
        cruise.astype(str) + "|" + meta["station"].astype(str),
        fallback,
    )


def hellinger_gram(otu, ids):
    x = otu[ids].T.to_numpy(float)
    rs = x.sum(axis=1, keepdims=True)
    x = np.divide(x, rs, out=np.zeros_like(x), where=rs > 0)
    x = np.sqrt(x)
    return x @ x.T


def trait_gram(traits, ids):
    x = traits.loc[ids, TRAITS].fillna(0).to_numpy(float)
    rs = x.sum(axis=1, keepdims=True)
    x = np.divide(x, rs, out=np.zeros_like(x), where=rs > 0)
    x = np.sqrt(x)
    return x @ x.T


def design(meta, mode):
    parts = [np.ones((len(meta), 1))]
    if mode in ("station", "station_season"):
        parts.append(pd.get_dummies(meta["station"].astype(str), drop_first=True, dtype=float).to_numpy())
        if mode == "station_season":
            parts.append(pd.get_dummies(meta["season"].astype(str), drop_first=True, dtype=float).to_numpy())
    elif mode == "profile":
        parts.append(pd.get_dummies(meta["profile"].astype(str), drop_first=True, dtype=float).to_numpy())
    else:
        raise ValueError(mode)
    return np.column_stack(parts)


def residual_maker(meta, mode):
    z = design(meta, mode)
    return np.eye(len(meta)) - z @ np.linalg.pinv(z)


def valid_thresholds(meta, min_spanning=MIN_SPANNING, depths=None, profile_labels=None):
    if depths is None:
        depths = meta["depth"].to_numpy(float)
    if profile_labels is None:
        profile_labels = meta["profile"].to_numpy()
    out = []
    for t in THRESHOLDS:
        shallow = depths <= t
        span = 0
        for p in np.unique(profile_labels):
            z = shallow[profile_labels == p]
            span += int(z.any() and (~z).any())
        if span >= min_spanning and shallow.sum() >= 2 and (~shallow).sum() >= 2:
            out.append((t, span))
    return out


def scan(K, meta, mode, min_spanning=MIN_SPANNING, depths=None, profile_labels=None):
    if depths is None:
        depths = meta["depth"].to_numpy(float)
    if profile_labels is None:
        profile_labels = meta["profile"].to_numpy()
    M = residual_maker(meta, mode)
    ss0 = float(np.trace(M @ K))
    rows = []
    for t, span in valid_thresholds(meta, min_spanning, depths, profile_labels):
        shallow = depths <= t
        g = (~shallow).astype(float)
        gr = M @ g
        den = float(gr @ gr)
        if den <= 1e-12:
            continue
        ss_inc = float(gr @ K @ gr / den)
        rows.append([t, ss_inc / ss0, int(shallow.sum()), int((~shallow).sum()), span])
    return pd.DataFrame(
        rows,
        columns=["threshold_m", "partial_R2", "n_shallow", "n_deep", "n_profiles_spanning"],
    )


def joint_bootstrap(Ks, meta, nboot, seed, mode, min_spanning=MIN_SPANNING):
    rng = np.random.default_rng(seed)
    profiles = meta["profile"].unique()
    parr = meta["profile"].to_numpy()
    row_indices = {p: np.where(parr == p)[0] for p in profiles}
    out = []
    for b in range(nboot):
        chosen = rng.choice(profiles, size=len(profiles), replace=True)
        idx, bp = [], []
        for j, p in enumerate(chosen):
            ii = row_indices[p]
            idx.extend(ii.tolist())
            bp.extend([f"{j}|{p}"] * len(ii))
        idx = np.asarray(idx, dtype=int)
        bp = np.asarray(bp)
        mb = meta.iloc[idx].copy()
        depths = mb["depth"].to_numpy(float)
        for name, K in Ks.items():
            kb = K[np.ix_(idx, idx)]
            s = scan(kb, mb, mode, min_spanning=min_spanning, depths=depths, profile_labels=bp)
            if len(s):
                r = s.loc[s["partial_R2"].idxmax()]
                out.append([b, name, r["threshold_m"], r["partial_R2"]])
    return pd.DataFrame(out, columns=["bootstrap", "layer", "optimum_m", "max_partial_R2"])


def bootstrap_summary(boot):
    rows = []
    for layer, z in boot.groupby("layer"):
        rows.append([
            layer,
            z["optimum_m"].median(),
            z["optimum_m"].quantile(0.025),
            z["optimum_m"].quantile(0.975),
            len(z),
        ])
    return pd.DataFrame(rows, columns=["layer", "median_optimum_m", "ci95_low_m", "ci95_high_m", "n_bootstrap"])


def profile_shuffle_maxstat(Ks, meta, mode, nperm, seed, min_spanning=MIN_SPANNING):
    M = residual_maker(meta, mode)
    d0 = meta["depth"].to_numpy(float)
    prof = meta["profile"].to_numpy()
    valid = np.array([t for t, _ in valid_thresholds(meta, min_spanning=min_spanning)], dtype=float)
    G0 = (d0[:, None] > valid[None, :]).astype(float)
    MG0 = M @ G0
    den0 = np.sum(G0 * MG0, axis=0)

    prep, observed = {}, {}
    for name, K in Ks.items():
        Kr = M @ K @ M
        ss0 = float(np.trace(M @ K))
        num = np.sum(G0 * (Kr @ G0), axis=0)
        r = np.divide(num, den0, out=np.full_like(num, np.nan), where=den0 > 1e-12) / ss0
        j = int(np.nanargmax(r))
        observed[name] = (valid[j], float(r[j]))
        prep[name] = (Kr, ss0)

    rng = np.random.default_rng(seed)
    pidx = {p: np.where(prof == p)[0] for p in np.unique(prof)}
    null = {name: np.empty(nperm) for name in Ks}
    for b in range(nperm):
        dp = d0.copy()
        for _, ii in pidx.items():
            dp[ii] = rng.permutation(dp[ii])
        G = (dp[:, None] > valid[None, :]).astype(float)
        MG = M @ G
        den = np.sum(G * MG, axis=0)
        for name, (Kr, ss0) in prep.items():
            num = np.sum(G * (Kr @ G), axis=0)
            r = np.divide(num, den, out=np.full_like(num, np.nan), where=den > 1e-12) / ss0
            null[name][b] = np.nanmax(r)

    rows = []
    for name, arr in null.items():
        t, r = observed[name]
        p = (1 + np.sum(arr >= r)) / (nperm + 1)
        rows.append([name, t, r, p, np.median(arr), np.quantile(arr, 0.95)])
    return pd.DataFrame(
        rows,
        columns=[
            "layer", "observed_optimum_m", "observed_max_partial_R2",
            "profile_shuffle_maxstat_p", "null_median", "null_95"
        ],
    )


def difference_summary(boot):
    wide = boot.pivot(index="bootstrap", columns="layer", values="optimum_m").dropna()
    rows = []
    for a, b in [("traits", "16S"), ("traits", "18S"), ("16S", "18S")]:
        d = wide[b] - wide[a]
        rows.append([
            a, b, np.median(d), np.quantile(d, 0.025), np.quantile(d, 0.975),
            np.mean(wide[a] < wide[b]), np.mean(wide[a] == wide[b]), np.mean(wide[a] > wide[b])
        ])
    return pd.DataFrame(rows, columns=[
        "layer_a", "layer_b", "median_b_minus_a_m", "ci95_low_m", "ci95_high_m",
        "prob_a_shallower_b", "prob_equal", "prob_a_deeper_b"
    ])


def load_inputs(root):
    o16 = pick(
        root / "data/raw/16S/otu_table_prok_clean.csv.gz",
        root / "data/raw/16S/otu_table_prok_clean.csv",
    )
    o18 = pick(
        root / "data/raw/18S/otu_table_euk_clean.csv.gz",
        root / "data/raw/18S/otu_table_euk_clean.csv",
    )
    m16p = root / "data/processed/metadata/16S_samples_CTD_chemistry.csv"
    m18p = root / "data/processed/metadata/18S_samples_CTD_chemistry.csv"
    trp = root / "data/processed/traits/sample_level_core_functional_traits.csv"
    return (
        pd.read_csv(o16).set_index("OTUID"),
        pd.read_csv(o18).set_index("OTUID"),
        pd.read_csv(m16p).set_index("sample-id"),
        pd.read_csv(m18p).set_index("sample-id"),
        pd.read_csv(trp).set_index("sample-id"),
    )


def main(root):
    root = Path(root)
    out = root / "results" / "fig02"
    out.mkdir(parents=True, exist_ok=True)
    otu16, otu18, m16, m18, traits = load_inputs(root)

    common = sorted(set(otu16.columns) & set(otu18.columns) & set(m16.index) & set(m18.index) & set(traits.index))
    meta = m18.loc[common].copy()
    meta = meta[meta["depth"] <= MAX_DEPTH].copy()
    meta["profile"] = profile_id(meta)
    ids = meta.index.tolist()

    audit = []
    for sid in ids:
        audit.append([
            sid,
            float(m16.loc[sid, "depth"]), float(m18.loc[sid, "depth"]),
            str(m16.loc[sid, "station"]), str(m18.loc[sid, "station"]),
            str(m16.loc[sid, "season"]), str(m18.loc[sid, "season"]),
        ])
    audit = pd.DataFrame(audit, columns=[
        "sample_id", "depth_16S", "depth_18S", "station_16S", "station_18S",
        "season_16S", "season_18S"
    ])
    audit["depth_abs_diff_m"] = (audit["depth_16S"] - audit["depth_18S"]).abs()
    audit["station_match"] = audit["station_16S"] == audit["station_18S"]
    audit["season_match"] = audit["season_16S"] == audit["season_18S"]
    audit.to_csv(out / "paired_metadata_audit.csv", index=False)

    Ks = {
        "16S": hellinger_gram(otu16, ids),
        "18S": hellinger_gram(otu18, ids),
        "traits": trait_gram(traits, ids),
    }

    scan_rows = []
    for name, K in Ks.items():
        s = scan(K, meta, "station_season")
        s["layer"] = name
        s["state"] = "pooled"
        scan_rows.append(s)
    for state in ["A_Winter", "Summer"]:
        mask = (meta["season"] == state).to_numpy()
        mm = meta.loc[mask].copy()
        idx = np.where(mask)[0]
        for name, K in Ks.items():
            s = scan(K[np.ix_(idx, idx)], mm, "station")
            s["layer"] = name
            s["state"] = state
            scan_rows.append(s)
    scans = pd.concat(scan_rows, ignore_index=True)
    scans.to_csv(out / "matched_boundary_scans.csv", index=False)
    best = scans.loc[scans.groupby(["state", "layer"])["partial_R2"].idxmax()].copy()
    best = best.sort_values(["state", "layer"]).reset_index(drop=True)
    best.to_csv(out / "matched_boundary_optima.csv", index=False)

    boot = joint_bootstrap(Ks, meta, N_BOOT, SEED, "station_season", min_spanning=BOOT_MIN_SPANNING)
    boot.to_csv(out / "matched_joint_profile_bootstrap.csv", index=False)
    bootstrap_summary(boot).to_csv(out / "matched_bootstrap_summary.csv", index=False)
    difference_summary(boot).to_csv(out / "matched_breakpoint_difference_summary.csv", index=False)

    state_boots = []
    for j, state in enumerate(["A_Winter", "Summer"], start=1):
        mask = (meta["season"] == state).to_numpy()
        mm = meta.loc[mask].copy()
        idx = np.where(mask)[0]
        kb = {name: K[np.ix_(idx, idx)] for name, K in Ks.items()}
        sb = joint_bootstrap(kb, mm, N_SEASON_BOOT, SEED + j, "station", min_spanning=BOOT_MIN_SPANNING)
        sb["state"] = state
        state_boots.append(sb)
    pd.concat(state_boots, ignore_index=True).to_csv(out / "matched_state_profile_bootstrap.csv", index=False)

    profile_shuffle_maxstat(Ks, meta, "station_season", N_PERM, SEED + 100).to_csv(
        out / "matched_profile_shuffle_maxstat.csv", index=False
    )

    nuisance_rows = []
    for mode in ["station", "station_season", "profile"]:
        for name, K in Ks.items():
            s = scan(K, meta, mode)
            r = s.loc[s["partial_R2"].idxmax()]
            nuisance_rows.append([mode, name, r["threshold_m"], r["partial_R2"]])
    pd.DataFrame(nuisance_rows, columns=["nuisance_model", "layer", "optimum_m", "partial_R2"]).to_csv(
        out / "matched_nuisance_sensitivity.csv", index=False
    )

    support_rows = []
    for minspan in SUPPORT_SENSITIVITY:
        for state in ["pooled", "A_Winter", "Summer"]:
            if state == "pooled":
                idx = np.arange(len(meta)); mm = meta.copy(); mode = "station_season"
            else:
                idx = np.where((meta["season"] == state).to_numpy())[0]
                mm = meta.iloc[idx].copy(); mode = "station"
            for name, K in Ks.items():
                s = scan(K[np.ix_(idx, idx)], mm, mode, min_spanning=minspan)
                r = s.loc[s["partial_R2"].idxmax()]
                support_rows.append([
                    minspan, state, name, r["threshold_m"], r["partial_R2"],
                    r["n_profiles_spanning"], s["threshold_m"].max()
                ])
    pd.DataFrame(support_rows, columns=[
        "min_spanning_profiles", "state", "layer", "optimum_m", "max_partial_R2",
        "profiles_spanning_at_optimum", "deepest_supported_threshold_m"
    ]).to_csv(out / "minimum_profile_sensitivity.csv", index=False)

    # Long 16S record: contextual validation only, never mixed into the matched bootstrap.
    m16full = m16[m16["depth"].notna() & (m16["depth"] <= MAX_DEPTH) & m16.index.isin(otu16.columns)].copy()
    m16full["profile"] = profile_id(m16full)
    ids16 = m16full.index.tolist()
    K16full = hellinger_gram(otu16, ids16)
    long_rows = []
    for state in ["pooled", "A_Winter", "Summer"]:
        if state == "pooled":
            idx = np.arange(len(m16full)); mm = m16full.copy(); mode = "station_season"
        else:
            idx = np.where((m16full["season"] == state).to_numpy())[0]
            mm = m16full.iloc[idx].copy(); mode = "station"
        s = scan(K16full[np.ix_(idx, idx)], mm, mode)
        s["state"] = state
        long_rows.append(s)
    pd.concat(long_rows, ignore_index=True).to_csv(out / "long_16S_2018_2026_scans.csv", index=False)

    pd.DataFrame([{
        "n_paired_all_depths": len(common),
        "n_paired_le650m": len(meta),
        "n_profiles_le650m": meta["profile"].nunique(),
        "n_winter_le650m": int((meta["season"] == "A_Winter").sum()),
        "n_summer_le650m": int((meta["season"] == "Summer").sum()),
        "minimum_spanning_profiles_primary": MIN_SPANNING,
        "minimum_spanning_profiles_bootstrap": BOOT_MIN_SPANNING,
        "n_metadata_station_mismatches_le650m": int((~audit["station_match"]).sum()),
        "max_depth_metadata_difference_m": float(audit["depth_abs_diff_m"].max()),
    }]).to_csv(out / "matched_analysis_sample_summary.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=Path(__file__).resolve().parents[2])
    args = ap.parse_args()
    main(args.root)
