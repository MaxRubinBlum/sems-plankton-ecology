#!/usr/bin/env python3
"""Prepare priority 99%-OTU centroid IDs for Phase-6 sequence validation."""

import argparse
import pandas as pd


def norm(x):
    x = str(x).strip()
    if x.lower() in {"", "nan", "unclassified", "uncultured", "unknown", "none"}:
        return ""
    if "__" in x:
        x = x.split("__", 1)[1]
    return x.replace("_", " ").strip()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--otu", required=True)
    ap.add_argument("--taxonomy", required=True)
    ap.add_argument("--mapping", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--out-table", required=True)
    ap.add_argument("--out-ids", required=True)
    args = ap.parse_args()

    otu = pd.read_csv(args.otu).set_index("OTUID")
    tax = pd.read_csv(args.taxonomy).set_index("OTUID").fillna("")
    mapping = pd.read_csv(args.mapping, sep="\t").set_index("feature_id")
    candidates = pd.read_csv(args.candidates)

    def label(row):
        for rank in ["Species", "Genus", "Family", "Order", "Class", "Phylum"]:
            value = norm(row.get(rank, ""))
            if value:
                return value
        return "Unclassified Eukaryota"

    def group(row):
        lineage = ";".join(norm(row.get(c, "")) for c in ["Phylum", "Class", "Order", "Family", "Genus", "Species"]).lower()
        if "metazoa" in lineage:
            return "Metazoa"
        if "syndin" in lineage or "dino-group-i" in lineage or "dino-group-ii" in lineage:
            return "Syndiniales"
        if "hematodinium" in lineage:
            return "Hematodinium"
        return "Other eukaryote"

    ann = pd.DataFrame(index=tax.index)
    ann["taxon_label"] = tax.apply(label, axis=1)
    ann["ecological_group"] = tax.apply(group, axis=1)
    ann["agg_key"] = ann.ecological_group + " | " + ann.taxon_label
    feature_key = {fid: f"{r.ecological_group} | {r.taxon_label}" for fid, r in mapping.iterrows()}
    total_reads = otu.sum(axis=1)

    rows = []
    requested_ids = []
    priority = candidates[candidates.tier.str.startswith(("Tier 1", "Tier 2"))].sort_values("evidence_score", ascending=False).head(20)

    for _, row in priority.iterrows():
        parasite_feature, partner_feature = row.edge_key.split("|")
        record = [row.edge_key, row.parasite, row.partner, row.tier, row.literature_class]
        for feature_id in [parasite_feature, partner_feature]:
            members = ann.index[ann.agg_key == feature_key[feature_id]].tolist()
            ordered = total_reads.loc[members].sort_values(ascending=False) if members else pd.Series(dtype=float)
            top = ordered.head(5)
            share = (top.iloc[0] / ordered.sum()) if len(top) and ordered.sum() > 0 else float("nan")
            record.extend([
                feature_id,
                len(members),
                top.index[0] if len(top) else "",
                share,
                ";".join(top.index),
            ])
            requested_ids.extend(top.index[:3])
        rows.append(record)

    columns = [
        "edge_key", "parasite", "partner", "tier", "literature_class",
        "parasite_feature_id", "parasite_constituent_OTU99s", "parasite_top_OTU99",
        "parasite_top_OTU99_read_share", "parasite_top5_OTU99s",
        "partner_feature_id", "partner_constituent_OTU99s", "partner_top_OTU99",
        "partner_top_OTU99_read_share", "partner_top5_OTU99s",
    ]
    pd.DataFrame(rows, columns=columns).to_csv(args.out_table, index=False)

    with open(args.out_ids, "w") as handle:
        for otuid in sorted(set(requested_ids)):
            handle.write(otuid + "\n")


if __name__ == "__main__":
    main()
