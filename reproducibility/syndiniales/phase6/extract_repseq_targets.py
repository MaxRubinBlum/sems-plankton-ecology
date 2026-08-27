#!/usr/bin/env python3
"""Validate rep-seqs-99.qza and extract priority 99%-OTU centroid sequences."""

import argparse
import hashlib
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd


def read_qza_fasta(qza_path):
    with zipfile.ZipFile(qza_path) as z:
        members = z.namelist()
        fasta_member = next(x for x in members if x.endswith("/data/dna-sequences.fasta"))
        fasta_text = z.read(fasta_member).decode()
        top_action = next(x for x in members if x.endswith("/provenance/action/action.yaml"))
        top_action_text = z.read(top_action).decode()
        provenance_actions = [z.read(x).decode() for x in members if x.endswith("/action/action.yaml")]

    seqs = {}
    current = None
    for line in fasta_text.splitlines():
        if line.startswith(">"):
            current = line[1:].split()[0]
            seqs[current] = ""
        else:
            seqs[current] += line.strip()
    return seqs, top_action_text, provenance_actions


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qza", required=True)
    ap.add_argument("--otu", required=True)
    ap.add_argument("--taxonomy", required=True)
    ap.add_argument("--target-ids", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    qza = Path(args.qza)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    seqs, top_action, provenance_actions = read_qza_fasta(qza)
    otu = pd.read_csv(args.otu).set_index("OTUID")
    tax = pd.read_csv(args.taxonomy).set_index("OTUID").fillna("")
    targets = [x.strip() for x in Path(args.target_ids).read_text().splitlines() if x.strip()]

    if set(otu.index) != set(seqs):
        raise ValueError("QZA sequence IDs and OTU-table feature IDs are not identical")

    missing = [x for x in targets if x not in seqs]
    if missing:
        raise ValueError(f"Missing requested OTUIDs: {missing}")

    provenance_blob = "\n".join(provenance_actions)
    clustered_99 = (
        "cluster_features_de_novo" in top_action
        and re.search(r"perc_identity:\s*0\.99", top_action) is not None
    )
    dada2_source = "denoise_paired" in provenance_blob
    hashed_source = "hashed_feature_ids: true" in provenance_blob

    fasta_out = out / "priority_99pct_OTU_centroid_sequences.fasta"
    with fasta_out.open("w") as handle:
        for otuid in targets:
            tr = tax.loc[otuid]
            lineage = ";".join(str(tr.get(c, "")) for c in ["Phylum", "Class", "Order", "Family", "Genus", "Species"])
            handle.write(f">{otuid} taxonomy={lineage}\n{seqs[otuid]}\n")

    rows = []
    for otuid in targets:
        tr = tax.loc[otuid]
        rows.append({
            "OTUID": otuid,
            "sequence_length": len(seqs[otuid]),
            "sequence": seqs[otuid],
            "total_reads": float(otu.loc[otuid].sum()),
            "sample_prevalence": float((otu.loc[otuid] > 0).mean()),
            **{c: tr[c] for c in ["Phylum", "Class", "Order", "Family", "Genus", "Species"]},
        })
    pd.DataFrame(rows).sort_values("total_reads", ascending=False).to_csv(
        out / "priority_99pct_OTU_centroid_sequence_metadata.csv", index=False
    )

    lengths = [len(seqs[x]) for x in targets]
    summary = pd.DataFrame([{
        "qza_file": qza.name,
        "qza_sha256": hashlib.sha256(qza.read_bytes()).hexdigest(),
        "representative_features": len(seqs),
        "OTU_table_features": len(otu),
        "feature_id_sets_identical": True,
        "target_sequences_requested": len(targets),
        "target_sequences_found": len(targets),
        "target_sequences_unique": len({seqs[x] for x in targets}),
        "target_length_min": min(lengths),
        "target_length_median": float(np.median(lengths)),
        "target_length_max": max(lengths),
        "q2_vsearch_cluster_features_de_novo_99pct": clustered_99,
        "source_q2_dada2_denoise_paired": dada2_source,
        "source_hashed_feature_ids": hashed_source,
    }])
    summary.to_csv(out / "rep_seqs_99_validation_summary.csv", index=False)


if __name__ == "__main__":
    main()
