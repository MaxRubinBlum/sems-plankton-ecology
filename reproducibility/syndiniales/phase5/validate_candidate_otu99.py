#!/usr/bin/env python3
"""Targeted validation of aggregate parasite–protist associations at 99%-OTU level.

The current 18S feature table was produced by DADA2 and then clustered de novo at
99% identity with q2-vsearch. Consequently, feature IDs are 99%-OTU centroids,
not unclustered ASVs. This script preserves the archived Phase-5 calculations
while naming the analytical unit correctly.
"""

import argparse
import numpy as np
import pandas as pd
from scipy.stats import rankdata, t
from statsmodels.stats.multitest import multipletests

RANKS = ["Species", "Genus", "Family", "Order", "Class", "Phylum"]


def norm(x):
    x = str(x).strip()
    if x.lower() in {"", "nan", "unclassified", "uncultured", "unknown", "none"}:
        return ""
    if "__" in x:
        x = x.split("__", 1)[1]
    return x.replace("_", " ").strip()


def label(row):
    for rank in RANKS:
        value = norm(row.get(rank, ""))
        if value:
            return value
    return "Unclassified Eukaryota"


def ecological_group(row):
    lineage = ";".join(norm(row.get(c, "")) for c in ["Phylum", "Class", "Order", "Family", "Genus", "Species"]).lower()
    if "metazoa" in lineage:
        return "Metazoa"
    if "syndin" in lineage or "dino-group-i" in lineage or "dino-group-ii" in lineage:
        return "Syndiniales"
    if "hematodinium" in lineage:
        return "Hematodinium"
    return "Other eukaryote"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--otu", required=True)
    ap.add_argument("--taxonomy", required=True)
    ap.add_argument("--mapping", required=True)
    ap.add_argument("--features-parent", required=True)
    ap.add_argument("--metadata-parent", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    otu = pd.read_csv(args.otu).set_index("OTUID")
    tax = pd.read_csv(args.taxonomy).set_index("OTUID").fillna("")
    mapping = pd.read_csv(args.mapping, sep="\t").set_index("feature_id")
    fagg = pd.read_csv(args.features_parent, sep="\t", index_col=0)
    meta = pd.read_csv(args.metadata_parent, sep="\t", index_col=0).loc[fagg.index]
    candidates = pd.read_csv(args.candidates)

    ann = pd.DataFrame(index=tax.index)
    ann["taxon_label"] = tax.apply(label, axis=1)
    ann["ecological_group"] = tax.apply(ecological_group, axis=1)
    ann["agg_key"] = ann.ecological_group + " | " + ann.taxon_label
    feature_key = {fid: f"{r.ecological_group} | {r.taxon_label}" for fid, r in mapping.iterrows()}

    raw = otu[fagg.index].astype(float)
    prevalence = (raw > 0).mean(axis=1)
    total_reads = raw.sum(axis=1)
    eligible = (prevalence >= 0.02) & (total_reads >= 50)

    def trim(ids):
        if not ids:
            return []
        s = total_reads.loc[ids].sort_values(ascending=False)
        cumulative = s.cumsum() / s.sum()
        keep = list(s.index[cumulative <= 0.90]) or [s.index[0]]
        if len(keep) < len(s):
            keep.append(s.index[len(keep)])
        return keep[:20]

    specs = []
    all_ids = set()
    for _, row in candidates.iterrows():
        parasite_feature, partner_feature = row.edge_key.split("|")
        parasite_ids = trim(ann.index[(ann.agg_key == feature_key[parasite_feature]) & eligible].tolist())
        partner_ids = trim(ann.index[(ann.agg_key == feature_key[partner_feature]) & eligible].tolist())
        if parasite_ids and partner_ids:
            specs.append((row.edge_key, row.parasite, row.partner, row.tier, parasite_ids, partner_ids))
            all_ids.update(parasite_ids)
            all_ids.update(partner_ids)

    # Residualize log1p OTU abundances against the exact FlashWeave metadata matrix.
    cov = meta.astype(float)
    cov = cov.loc[:, cov.nunique() > 1]
    design = np.column_stack([np.ones(len(cov)), cov.values])
    ids = list(all_ids)
    y = np.log1p(raw.loc[ids, fagg.index].T.values)
    residual = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    pos = {x: i for i, x in enumerate(ids)}

    ranked = np.column_stack([rankdata(residual[:, j]) for j in range(residual.shape[1])])
    ranked = (ranked - ranked.mean(0)) / ranked.std(0, ddof=1)
    n = len(cov)

    rows = []
    for edge, parasite_name, partner_name, tier, parasite_ids, partner_ids in specs:
        for parasite_otu in parasite_ids:
            for partner_otu in partner_ids:
                copresent = int(((raw.loc[parasite_otu, fagg.index] > 0) & (raw.loc[partner_otu, fagg.index] > 0)).sum())
                if copresent < 20:
                    continue
                rho = float(np.dot(ranked[:, pos[parasite_otu]], ranked[:, pos[partner_otu]]) / (n - 1))
                statistic = rho * np.sqrt((n - 2) / max(1e-12, 1 - rho * rho))
                p = float(2 * t.sf(abs(statistic), df=n - 2))
                rows.append([edge, parasite_name, partner_name, tier, parasite_otu, partner_otu, copresent, rho, p])

    result = pd.DataFrame(rows, columns=[
        "aggregate_edge", "parasite_taxon", "partner_taxon", "aggregate_tier",
        "parasite_OTU99", "partner_OTU99", "n_copresent", "residual_spearman_rho", "p"
    ])
    result["q_global"] = multipletests(result.p, method="fdr_bh")[1]
    result["supported"] = (result.q_global < 0.05) & (result.residual_spearman_rho >= 0.20)
    result.to_csv(args.out_prefix + "_pair_tests.csv", index=False)

    summary = result.groupby("aggregate_edge").agg(
        OTU99_pairs_tested=("supported", "size"),
        OTU99_pairs_FDR_supported=("supported", "sum")
    ).reset_index()
    summary.to_csv(args.out_prefix + "_summary.csv", index=False)


if __name__ == "__main__":
    main()
