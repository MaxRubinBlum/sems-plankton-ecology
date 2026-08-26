#!/usr/bin/env python3
"""Phase 3: test Syndiniales partner-set turnover across the 220-m boundary.

Input is the Phase-2 `cross_regime_edge_retesting.csv` table produced by
`reproducibility/syndiniales/phase2/analyze_phase2.py`.

The strict null preserves, for each parasite lineage and regime, observed degree,
partner-guild composition, regime-specific candidate availability, and the
parasite-specific candidate-partner universe.
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--cross-regime", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--seed", type=int, default=20260826)
    p.add_argument("--permutations", type=int, default=20000)
    return p.parse_args()


def guild_counts(ids: set[str], group_map: dict[str, str]) -> Counter:
    return Counter(group_map.get(x, "Other protists") for x in ids)


def main() -> None:
    args = parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    x = pd.read_csv(args.cross_regime)
    required = {
        "parasite_id", "parasite", "partner_id", "partner", "partner_group",
        "upper_supported", "deep_supported", "upper_both_features", "deep_both_features",
        "upper_resid_rho", "deep_resid_rho",
    }
    missing = required - set(x.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    for c in ["upper_supported", "deep_supported", "upper_both_features", "deep_both_features"]:
        if x[c].dtype == object:
            x[c] = x[c].astype(str).str.lower().map({"true": True, "false": False})
        x[c] = x[c].fillna(False).astype(bool)

    partner_group = dict(zip(x.partner_id, x.partner_group))
    rng = np.random.default_rng(args.seed)
    B = args.permutations

    lineage_rows = []
    shared_rows = []
    null_rows = []
    global_null = np.zeros(B, dtype=int)
    observed_global = 0

    for pid, sub in x.groupby("parasite_id", sort=True):
        name = sub.parasite.iloc[0]
        upper = set(sub.loc[sub.upper_supported, "partner_id"])
        deep = set(sub.loc[sub.deep_supported, "partner_id"])
        shared = upper & deep
        union = upper | deep

        same_sign = 0
        reversed_sign = 0
        for q in sorted(shared):
            r = sub[sub.partner_id == q].iloc[0]
            behavior = "same sign" if np.sign(r.upper_resid_rho) == np.sign(r.deep_resid_rho) else "sign reversal"
            same_sign += behavior == "same sign"
            reversed_sign += behavior == "sign reversal"
            shared_rows.append([
                pid, name, q, r.partner, r.partner_group,
                r.upper_resid_rho, r.deep_resid_rho, behavior,
            ])

        lineage_rows.append([
            pid, name, len(upper), len(deep), len(shared), same_sign, reversed_sign,
            len(shared) / len(union) if union else np.nan,
        ])

        if not upper or not deep:
            continue

        candidate_union = set(sub.partner_id)
        available_upper = set(sub.loc[sub.upper_both_features, "partner_id"]) & candidate_union
        available_deep = set(sub.loc[sub.deep_both_features, "partner_id"]) & candidate_union
        c_upper = guild_counts(upper, partner_group)
        c_deep = guild_counts(deep, partner_group)

        overlap_null = np.zeros(B, dtype=int)
        for guild in set(c_upper) | set(c_deep):
            ku, kd = c_upper.get(guild, 0), c_deep.get(guild, 0)
            if ku == 0 or kd == 0:
                continue
            gu = {q for q in available_upper if partner_group.get(q, "Other protists") == guild}
            gd = {q for q in available_deep if partner_group.get(q, "Other protists") == guild}
            Nu, Nd, Nc = len(gu), len(gd), len(gu & gd)
            if Nu == 0 or Nd == 0 or Nc == 0:
                continue
            ku = min(ku, Nu)
            kd = min(kd, Nd)

            common_selected_upper = rng.hypergeometric(Nc, Nu - Nc, ku, size=B)
            guild_overlap = rng.hypergeometric(common_selected_upper, Nd - common_selected_upper, kd)
            overlap_null += guild_overlap

        observed_shared = len(shared)
        observed_jaccard = observed_shared / len(union)
        null_jaccard = overlap_null / (len(upper) + len(deep) - overlap_null)
        p_less = (1 + np.sum(null_jaccard <= observed_jaccard)) / (B + 1)
        p_more = (1 + np.sum(null_jaccard >= observed_jaccard)) / (B + 1)

        null_rows.append([
            pid, name, len(upper), len(deep), len(candidate_union),
            len(available_upper), len(available_deep), observed_shared, observed_jaccard,
            overlap_null.mean(), null_jaccard.mean(), np.quantile(null_jaccard, 0.025),
            np.quantile(null_jaccard, 0.975), p_less, p_more,
        ])
        global_null += overlap_null
        observed_global += observed_shared

    lineage = pd.DataFrame(lineage_rows, columns=[
        "parasite_id", "parasite", "upper_degree", "deep_degree", "shared_partners",
        "shared_same_sign", "shared_sign_reversal", "partner_jaccard",
    ])
    shared_df = pd.DataFrame(shared_rows, columns=[
        "parasite_id", "parasite", "partner_id", "partner", "partner_group",
        "upper_resid_rho", "deep_resid_rho", "sign_behavior",
    ])
    null = pd.DataFrame(null_rows, columns=[
        "parasite_id", "parasite", "upper_degree", "deep_degree", "candidate_union_size",
        "candidate_available_upper", "candidate_available_deep", "observed_shared",
        "observed_jaccard", "null_mean_shared", "null_mean_jaccard", "null_CI2.5",
        "null_CI97.5", "p_less_overlap_than_null", "p_more_overlap_than_null",
    ])
    if len(null):
        null["q_less_overlap_than_null"] = multipletests(null.p_less_overlap_than_null, method="fdr_bh")[1]
        null["q_more_overlap_than_null"] = multipletests(null.p_more_overlap_than_null, method="fdr_bh")[1]

    global_summary = pd.DataFrame([{
        "eligible_parasite_lineages": len(null),
        "observed_shared_partner_edges": observed_global,
        "null_mean_shared_partner_edges": global_null.mean(),
        "null_CI2.5": np.quantile(global_null, 0.025),
        "null_CI97.5": np.quantile(global_null, 0.975),
        "p_less_overlap_than_candidate_universe_null": (1 + np.sum(global_null <= observed_global)) / (B + 1),
        "p_more_overlap_than_candidate_universe_null": (1 + np.sum(global_null >= observed_global)) / (B + 1),
        "permutations": B,
        "seed": args.seed,
    }])

    lineage.to_csv(out / "syndiniales_lineage_partner_turnover.csv", index=False)
    shared_df.to_csv(out / "shared_partner_edges_and_sign_behavior.csv", index=False)
    null.to_csv(out / "partner_turnover_candidate_universe_null_tests.csv", index=False)
    global_summary.to_csv(out / "global_partner_turnover_candidate_universe_null_test.csv", index=False)

    summary = pd.DataFrame([
        ["upper_residual_supported_edges", int(x.upper_supported.sum())],
        ["deep_residual_supported_edges", int(x.deep_supported.sum())],
        ["shared_residual_supported_edges", int((x.upper_supported & x.deep_supported).sum())],
        ["shared_edges_same_sign", int((shared_df.sign_behavior == "same sign").sum())],
        ["shared_edges_sign_reversal", int((shared_df.sign_behavior == "sign reversal").sum())],
        ["comparable_parasite_lineages", len(null)],
        ["observed_shared_partner_edges", observed_global],
        ["null_mean_shared_partner_edges", float(global_null.mean())],
        ["global_p_less_overlap", float(global_summary.iloc[0].p_less_overlap_than_candidate_universe_null)],
        ["global_p_more_overlap", float(global_summary.iloc[0].p_more_overlap_than_candidate_universe_null)],
        ["lineages_FDR_excess_turnover", int((null.q_less_overlap_than_null < 0.05).sum())],
        ["lineages_FDR_excess_conservation", int((null.q_more_overlap_than_null < 0.05).sum())],
    ], columns=["metric", "value"])
    summary.to_csv(out / "summary.csv", index=False)


if __name__ == "__main__":
    main()
