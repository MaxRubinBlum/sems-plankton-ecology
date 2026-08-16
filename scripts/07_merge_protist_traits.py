#!/usr/bin/env python3
"""Merge Ramond-derived traits with literature-curated deep-sea extensions.

The curated extension must contain one row per taxon-trait assertion with
reference and confidence. Missing traits remain NA. This script does not
infer traits from depth, abundance or 16S-18S associations.
"""

import argparse
from pathlib import Path
import pandas as pd

REQUIRED_EXTENSION = {
    "taxon", "trait_name", "trait_state", "evidence_taxon", "evidence_rank",
    "reference_doi", "reference_citation", "confidence"
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ramond", required=True, help="Ramond trait table matched to study taxa")
    ap.add_argument("--extension", required=True, help="Curated deep-sea trait assertions")
    ap.add_argument("--profiles", required=True, help="Taxon abundance/depth profiles from script 06")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    ramond = pd.read_csv(args.ramond)
    ext = pd.read_csv(args.extension)
    profiles = pd.read_csv(args.profiles)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    missing = REQUIRED_EXTENSION - set(ext.columns)
    if missing:
        raise ValueError(f"Extension table missing required columns: {sorted(missing)}")

    # Reject unsupported assertions before they enter the analysis matrix.
    assertions = ext[
        ext["trait_name"].notna()
        & ext["trait_state"].notna()
        & (ext["trait_name"].astype(str).str.strip() != "")
        & (ext["trait_state"].astype(str).str.strip() != "")
    ].copy()

    duplicate = assertions.duplicated(["taxon", "trait_name"], keep=False)
    if duplicate.any():
        d = assertions.loc[duplicate, ["taxon", "trait_name"]]
        raise ValueError("Multiple extension states for the same taxon/trait:\n" + d.to_string(index=False))

    # Convert extension assertions from long to wide.
    ext_wide = assertions.pivot(index="taxon", columns="trait_name", values="trait_state")

    if "taxon" not in ramond.columns:
        raise ValueError("Ramond matched table must contain a 'taxon' column")
    ramond = ramond.set_index("taxon")

    # Ramond is primary; extensions fill only missing cells unless explicitly
    # curated outside this script after biological review.
    merged = ramond.copy()
    for trait in ext_wide.columns:
        if trait not in merged:
            merged[trait] = pd.NA
        merged[trait] = merged[trait].combine_first(ext_wide[trait])

    merged.reset_index().to_csv(out / "protist_traits_merged.csv", index=False)
    assertions.to_csv(out / "deep_trait_assertions_used.csv", index=False)

    # Coverage is abundance-weighted and therefore reported separately from
    # the trait matrix itself.
    if "taxon" in profiles.columns:
        p = profiles.set_index("taxon")
    else:
        p = profiles.set_index(profiles.columns[0])
    annotated = merged.notna().any(axis=1)
    coverage_rows = []
    for col in [c for c in p.columns if c.startswith(("A_", "B_", "C_", "D_", "F_", "G_", "H_"))]:
        common = p.index.intersection(annotated.index)
        total = p.loc[common, col].sum()
        covered = p.loc[common[annotated.reindex(common).fillna(False)], col].sum()
        coverage_rows.append({"layer": col, "annotated_abundance_pct": covered, "total_abundance_pct": total,
                              "coverage_fraction": covered / total if total else pd.NA})
    pd.DataFrame(coverage_rows).to_csv(out / "trait_coverage_by_depth.csv", index=False)


if __name__ == "__main__":
    main()
