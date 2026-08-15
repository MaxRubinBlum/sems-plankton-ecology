# Taxonomy provenance

## Prokaryotic GTDB taxonomy

GTDB release 226 assignments were generated in QIIME 2 using the sklearn classifier:

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

The exported taxonomy table is retained alongside the original SILVA-based 16S taxonomy. Neither taxonomy table should be overwritten; downstream comparisons should join assignments by feature/ASV identifier.

## 18S taxonomy

The original 18S taxonomy is retained, including Metazoa. The main microbial-eukaryote dataset is produced by filtering Metazoa and other excluded groups reproducibly from the original table rather than altering the raw file.
