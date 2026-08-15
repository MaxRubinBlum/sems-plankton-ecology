# Hybrid SILVA–GTDB226 prokaryotic annotation

This project retains both the original SILVA-style 16S taxonomy and a second classification against GTDB release 226. The canonical ecological annotation is a **hybrid crosswalk**, not a wholesale replacement of SILVA by GTDB.

## GTDB classification provenance

```bash
qiime feature-classifier classify-sklearn \
  --i-classifier /media/bioinf/Data/GTDB_qiime2/last/gtdb-226-v4-classifier.qza \
  --i-reads rep-seqs-99.qza \
  --o-classification taxonomyGTDB226.qza \
  --p-n-jobs 48 \
  --verbose

qiime tools export \
  --input-path taxonomyGTDB226.qza \
  --output-path ./gtdb226
```

The SILVA and GTDB tables contain exactly the same 24,909 ASV IDs.

## Annotation strategy

The crosswalk keeps the full taxonomy from both databases, GTDB classifier confidence, a `preferred_label`, an `ecological_guild`, and `annotation_source` describing how the preferred annotation was chosen.

The guiding rule is to preserve established marine ecological names when they carry information that is lost by the GTDB amplicon classifier, while using GTDB names where genome-based nomenclature improves resolution.

### Ecologically curated groups

- **SAR11:** retain the SILVA clade designation where available (for example Clade Ia, Ib, II, IV) and append a GTDB genus when informative (`Pelagibacter`, `Pelagibacter_A`, `AG-414-E02`, etc.).
- **SAR86:** retain SAR86 as the ecological umbrella and append GTDB subdivisions such as `D2472`, `TMED112`, `AEGEAN-183`, and `AG-339-G14`.
- **AEGEAN-169:** retain the ecological name and append GTDB HIMB59-radiation labels where resolved.
- **Nitrosopumilaceae/AOA:** prefer GTDB genus-level assignments. Most abundant ASVs resolve as `Nitrosopelagicus`.
- **SAR202:** retain SAR202 as the ecological umbrella and append GTDB genera where available (`Lucifugimonas`, `UBA3495`, `TMED-70`, etc.).
- **SAR406/Marinimicrobia:** retain SAR406/Marinimicrobia and append GTDB genera such as `Marinisoma`, `TCS55`, and `UBA2126`.
- **SAR324:** retain the SAR324 ecological label and append GTDB resolution where available.
- **Marine Group II / Poseidoniales:** retain the legacy ecological name while recording GTDB Poseidoniales/genus assignments.
- **Marine Group III:** retain the legacy ecological label and append GTDB resolution.
- **SUP05/sulfur oxidizers:** use GTDB names such as `Pseudothioglobus` where resolved but retain sulfur-oxidizer guild assignment.
- **Methylococcales:** assign the aerobic-methanotroph guild conservatively at order/family/genus level.
- **IheB2-23:** explicitly curated as **UBA1147 / IheB2-23**, a genomically supported aerobic methanotroph. This identity is not attributed to the GTDB226 sklearn output; the dominant ASV is classified by GTDB226 only to `Methylococcales`.

## Why hybrid taxonomy is necessary

GTDB226 is more conservative at lower ranks in this amplicon classification. In the complete ASV set, GTDB assigns 71.1% to order and 50.8% to genus, versus 100% and 80.3% respectively in the original SILVA table. Abundance-weighted coverage is higher (88.2% to order and 78.9% to genus), showing that dominant planktonic ASVs are generally better resolved.

However, some highly abundant and ecologically familiar marine lineages become less informative under GTDB classification. For example, the most abundant SILVA SAR11 Clade Ia ASV (~7.1% of selected reads) is returned only as `d__Bacteria` by the GTDB226 sklearn classifier. The hybrid crosswalk therefore preserves both systems and records the provenance of each ecological label.

## Output

The full generated table is named:

`prokaryote_taxonomy_crosswalk_SILVA_GTDB226_curated.csv`

A smaller key-lineage subset contains only ASVs assigned to one of the explicitly curated ecological guilds.
