# SILVA vs GTDB release 226: initial comparison

Both taxonomy files contain the same 24,909 ASV identifiers and join one-to-one.

## Assignment depth

GTDB226 is substantially more conservative than the original SILVA taxonomy at lower ranks. ASV-level assignment fractions are:

| Rank | SILVA | GTDB226 |
|---|---:|---:|
| Phylum | 100.0% | 82.6% |
| Class | 100.0% | 81.0% |
| Order | 100.0% | 71.1% |
| Family | 85.1% | 65.1% |
| Genus | 80.3% | 50.8% |
| Species | 52.4% | 32.8% |

When weighted by read abundance in samples with `select == a`, GTDB coverage is higher than these ASV counts suggest because abundant ASVs are generally better classified: 88.2% of reads receive an order and 78.9% a genus in GTDB226.

## Nomenclature agreement

Exact taxon-name agreement is low at several ranks because GTDB and SILVA use different nomenclature and delimit ecological clades differently. This should not be treated as classification error. Examples include Proteobacteria/Pseudomonadota and the restructuring of SAR11, SAR202, SAR406 and marine archaeal groups.

## Ecologically important observations

- SILVA `Nitrosopumilaceae` is resolved predominantly as GTDB `Nitrosopelagicus`; this is a clear gain in ecological resolution.
- SILVA SAR11 clades split among `Pelagibacter`, `Pelagibacter_A`, `AG-414-E02`, `AAA240-E13`, `AG-422-B15` and related GTDB lineages. However, the single most abundant SILVA Clade Ia ASV (~7.1% of selected reads) is classified only as `d__Bacteria` by the GTDB sklearn classifier, so GTDB is not uniformly more informative.
- SILVA `AEGEAN-169 marine group` maps mainly to GTDB HIMB59-family lineages such as `HIMB59` and `AG-337-I02`.
- SILVA SAR86 splits into multiple GTDB genera including `D2472`, `TMED112`, `AEGEAN-183`, and `AG-339-G14`.
- SILVA SAR202 and SAR406/Marinimicrobia are decomposed into many GTDB genome-defined genera, which can be useful for lineage-level analyses but requires a crosswalk to retain ecological continuity with the literature.
- Methylococcales abundance is almost identical between systems at the order level (~0.397% of selected reads). The dominant SILVA `IheB2-23` ASV (~0.281% of reads) is classified by GTDB226 only to `o__Methylococcales` (confidence ~0.775), not to UBA1147. Therefore the UBA1147 = IheB2-23 biological identification should be retained as a curated literature/genome-supported annotation rather than claimed as an output of the GTDB amplicon classifier.

## Recommendation

Do not replace SILVA with GTDB. Retain both taxonomies and build a curated ecological crosswalk. Use GTDB for genome-consistent nomenclature and finer resolution where supported, while preserving SILVA ecological clade labels when they are more informative for marine-plankton interpretation.
