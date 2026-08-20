"""Reproduce Figure 1 analytical panels and export plotted data.

Figure 1 contains the sampling design, Bray-Curtis PCoA for 16S and 18S,
and vertical Shannon-diversity profiles. The CTD-only hydrographic PCA was
moved to the supplementary workflow.

Expected repository inputs (see inputs.tsv): cleaned 16S/18S feature tables,
integrated sample metadata, and original sample metadata containing water depth.
Run from any directory; paths are repo-relative.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from statsmodels.nonparametric.smoothers_lowess import lowess

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'data'
OUT = ROOT / 'results' / 'fig01'
OUT.mkdir(parents=True, exist_ok=True)

DEPTH_ORDER = [
    'A_surface','B_nearsurface','C_DCM','D_below_DCM',
    'F_300-600','G_below_600','H_near_bottom'
]


def pcoa_bray(otu, sample_ids):
    """Bray-Curtis PCoA on sample-wise relative abundances."""
    ids = [x for x in sample_ids if x in otu.columns]
    x = otu[ids].T.astype(float).to_numpy()
    x = x / np.maximum(x.sum(1, keepdims=True), 1)
    d = squareform(pdist(x, metric='braycurtis'))
    n = len(ids)
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ (d * d) @ j
    eigval, eigvec = np.linalg.eigh(b)
    ix = np.argsort(eigval)[::-1]
    eigval = eigval[ix]
    eigvec = eigvec[:, ix]
    pos = eigval > 0
    coords = eigvec[:, pos] * np.sqrt(eigval[pos])
    pct = 100 * eigval[pos] / eigval[pos].sum()
    return pd.DataFrame(coords[:, :2], index=ids, columns=['PCo1', 'PCo2']), pct[:2]


def shannon_from_table(otu, sample_ids):
    """Natural-log Shannon diversity for each sample from ASV counts."""
    ids = [x for x in sample_ids if x in otu.columns]
    x = otu[ids].T.astype(float)
    p = x.div(x.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    arr = p.to_numpy()
    h = -(np.where(arr > 0, arr * np.log(arr), 0)).sum(axis=1)
    return pd.Series(h, index=ids, name='Shannon')


def lowess_curve(df, frac):
    z = df.dropna(subset=['depth', 'Shannon']).sort_values('depth')
    sm = lowess(z['Shannon'].to_numpy(), z['depth'].to_numpy(),
                frac=frac, it=2, return_sorted=True)
    return pd.DataFrame({'depth': sm[:, 0], 'Shannon_smooth': sm[:, 1]})


def main():
    m16 = pd.read_csv(DATA / '16S_samples_CTD_chemistry.csv')
    m18 = pd.read_csv(DATA / '18S_samples_CTD_chemistry.csv')
    raw18 = pd.read_csv(DATA / 'metadata_euk_clean.csv')

    all_shannon = []
    for marker, meta in [('16S', m16), ('18S', m18)]:
        otu = pd.read_csv(DATA / ('otu_table_prok_clean.csv' if marker == '16S' else 'otu_table_euk_clean.csv')).set_index('OTUID')
        ids = meta['sample-id'].dropna().tolist()

        co, pct = pcoa_bray(otu, ids)
        dat = meta.set_index('sample-id').reindex(co.index)[['depth','ds2','season','station','cruise']].join(co)
        dat.to_csv(OUT / f'{marker}_bray_pcoa_scores.csv')
        pd.DataFrame({'axis':['PCo1','PCo2'], 'variance_percent':pct}).to_csv(
            OUT / f'{marker}_bray_pcoa_variance.csv', index=False)

        h = shannon_from_table(otu, ids)
        hd = meta.set_index('sample-id').reindex(h.index)[['depth','ds2','season','station','cruise']].copy()
        hd['Shannon'] = h
        hd['domain'] = marker
        all_shannon.append(hd.reset_index())

    sh = pd.concat(all_shannon, ignore_index=True)
    sh.to_csv(OUT / 'Shannon_by_sample.csv', index=False)

    smooth_rows = []
    for marker in ['16S','18S']:
        d = sh[sh['domain'] == marker].copy()
        c = lowess_curve(d, frac=0.28)
        c['domain'] = marker
        c['season'] = 'Combined'
        smooth_rows.append(c)
        for season in ['A_Winter','Summer']:
            q = d[d['season'] == season]
            if len(q) > 10:
                s = lowess_curve(q, frac=0.34)
                s['domain'] = marker
                s['season'] = 'Winter' if season == 'A_Winter' else 'Summer'
                smooth_rows.append(s)
    pd.concat(smooth_rows, ignore_index=True).to_csv(OUT / 'Shannon_lowess_curves.csv', index=False)

    # Sampling-design table used for the left panel.
    keep = ['sample-id','station','station_norm','depth','ds2','season','cruise']
    sampling = m18[keep].merge(
        raw18[['sample-id','waterdepth']].drop_duplicates('sample-id'),
        on='sample-id', how='left'
    )
    sampling.to_csv(OUT / 'sampling_design_data.csv', index=False)

    print('16S PCoA variance:', pd.read_csv(OUT/'16S_bray_pcoa_variance.csv')['variance_percent'].tolist())
    print('18S PCoA variance:', pd.read_csv(OUT/'18S_bray_pcoa_variance.csv')['variance_percent'].tolist())


if __name__ == '__main__':
    main()
