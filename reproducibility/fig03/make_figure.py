#!/usr/bin/env python3
"""Reproduce manuscript Figure 3 from repository inputs.

Figure 3 links taxonomic succession to the matched transition analysis in
Figure 2. Panels a-b use the full marker-specific records to display abundant,
depth-responsive taxa. Panel c uses the same five-component trophic
normalization as Figure 2 and shows pooled transition references at 110 m
(trophic representation) and 220 m (shared matched 16S/18S composition).
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / 'data' / 'raw'
META = ROOT / 'data' / 'processed' / 'metadata'
TRAIT_FILE = ROOT / 'data' / 'processed' / 'traits' / 'sample_level_core_functional_traits.csv'
OUT = ROOT / 'results' / 'fig03'
OUT.mkdir(parents=True, exist_ok=True)

DEPTH_ORDER = ['A_surface','B_nearsurface','C_DCM','D_below_DCM','F_300-600','G_below_600','H_near_bottom']
DEPTH_NAMES = ['Surface','Near-surf.','DCM','Below DCM','300–600','>600','Near-bottom']
SEASONS = ['A_Winter','Summer']
LAYER_DEPTH = np.array([1,50,115,190,375,1000,1450], dtype=float)
TRAITS = [
    ('bona_fide_phototrophy','Bona fide\nphototrophy'),
    ('constitutive_mixotrophy','Constitutive\nmixotrophy'),
    ('parasitism','Parasitism'),
    ('phagotrophy_radiolaria','Radiolarian\nphagotrophy'),
    ('diplonemid_heterotrophy','Diplonemid\nheterotrophy'),
]
POOLED_TROPHIC_OPTIMUM_M = 110
POOLED_TAXONOMIC_OPTIMUM_M = 220


def strip_prefix(x):
    x = str(x or '').strip()
    if '__' in x:
        x = x.split('__', 1)[1]
    return x.strip()


def unresolved_placeholder(x):
    raw = strip_prefix(x)
    if not raw:
        return True
    if re.search(r'(_X+|_XX+|_XXX+)$', raw):
        return True
    return raw.lower() in {'uncultured','unclassified','unknown','nan'}


def display_tax(x):
    x = strip_prefix(x).replace('_', ' ').strip()
    replacements = {
        'Marinimicrobia (SAR406 clade)': 'Marinimicrobia (SAR406)',
        'SAR86 clade': 'SAR86',
        'SAR11 clade': 'SAR11',
        'SAR202 clade': 'SAR202',
        'SAR324 clade(Marine group B)': 'SAR324',
        'Acantharea 1': 'Acantharea',
        'Dino-Group-I': 'Syndiniales Group I',
        'Dino-Group-II': 'Syndiniales Group II',
        'Dino-Group-IV': 'Syndiniales Group IV',
        'RAD-B': 'RAD-B',
    }
    return replacements.get(x, x)


def load_marker(marker):
    if marker == '16S':
        otu_path = RAW / '16S' / 'otu_table_prok_clean.csv.gz'
        tax_path = RAW / '16S' / 'taxonomy_prok_clean.csv'
    else:
        otu_path = RAW / '18S' / 'otu_table_euk_clean.csv.gz'
        tax_path = RAW / '18S' / 'taxonomy_euk_clean.csv'
    otu = pd.read_csv(otu_path).set_index('OTUID')
    tax = pd.read_csv(tax_path).set_index('OTUID').fillna('')
    meta = pd.read_csv(META / f'{marker}_samples_CTD_chemistry.csv')
    return otu, tax, meta


def aggregate(marker):
    otu, tax, meta = load_marker(marker)
    mm = meta[
        meta.ds2.isin(DEPTH_ORDER) & meta.season.isin(SEASONS)
    ].dropna(subset=['depth']).set_index('sample-id')
    ids = [s for s in mm.index if s in otu.columns]
    x = otu[ids].T.astype(float)
    x = x.div(x.sum(axis=1).replace(0, np.nan), axis=0) * 100

    labels = []
    for oid in otu.index:
        r = tax.loc[oid]
        chosen = ''
        for rank in ['Order','Class','Family']:
            if rank in r.index and not unresolved_placeholder(r[rank]):
                chosen = strip_prefix(r[rank])
                break
        if not chosen:
            chosen = 'Unclassified'
        labels.append(display_tax(chosen))

    a = x.T.copy()
    a['taxon'] = pd.Series(labels, index=otu.index).reindex(a.index).values
    return a.groupby('taxon').sum().T, mm.loc[ids]


def select_and_split(g, mm, n=14, min_peak=.20):
    pooled = g.groupby(mm.ds2).median().reindex(DEPTH_ORDER).T.fillna(0)
    eligible = pooled.loc[pooled.max(axis=1) >= min_peak]
    score = (eligible.max(axis=1)-eligible.min(axis=1)) * np.log1p(eligible.max(axis=1)*5)
    p = pooled.loc[score.nlargest(n).index]
    centroid = (p.values * LAYER_DEPTH).sum(axis=1) / (p.values.sum(axis=1)+1e-9)
    p = p.iloc[np.argsort(centroid)]

    seasonal = {}
    for season in SEASONS:
        ms = mm[mm.season == season]
        seasonal[season] = (
            g.loc[ms.index].groupby(ms.ds2).median().reindex(DEPTH_ORDER).T.fillna(0).reindex(p.index)
        )

    w, s = seasonal['A_Winter'], seasonal['Summer']
    both = np.c_[w.values, s.values]
    mu = both.mean(axis=1, keepdims=True)
    sd = both.std(axis=1, ddof=1, keepdims=True)
    sd[sd == 0] = 1
    zw = pd.DataFrame((w.values-mu)/sd, index=p.index, columns=DEPTH_ORDER)
    zs = pd.DataFrame((s.values-mu)/sd, index=p.index, columns=DEPTH_ORDER)
    return p, w, s, zw, zs


def trait_summaries(df):
    trait_cols = [x[0] for x in TRAITS]
    vals = df[trait_cols].fillna(0).to_numpy(float)
    rs = vals.sum(axis=1, keepdims=True)
    vals = np.divide(vals, rs, out=np.zeros_like(vals), where=rs > 0) * 100
    for j, trait in enumerate(trait_cols):
        df[trait + '_selected_pct'] = vals[:, j]

    grid = np.geomspace(8, 1700, 120)
    rows = []
    for trait in trait_cols:
        vcol = trait + '_selected_pct'
        for season in SEASONS:
            z = df[df.season == season].dropna(subset=[vcol, 'depth'])
            dep = z.depth.to_numpy(float)
            val = z[vcol].to_numpy(float)
            for focal in grid:
                v = val[np.abs(dep-focal) <= max(35, .22*focal)]
                if len(v) >= 6:
                    rows.append([
                        trait, season, focal, len(v),
                        np.median(v), np.quantile(v, .25), np.quantile(v, .75)
                    ])
    return pd.DataFrame(rows, columns=[
        'trait','season','depth_m','n','median_selected_pct','q25_selected_pct','q75_selected_pct'
    ])


def bubble_size(v):
    return 10 + 20 * np.sqrt(np.maximum(v, 0))


def plot_heatmap(ax, zw, zs, mw, ms, panel_letter, norm):
    data = np.concatenate([zw.values, zs.values], axis=1)
    maxab = pd.concat([mw, ms], axis=1).max(axis=1).reindex(zw.index)
    im = ax.imshow(data, aspect='auto', norm=norm)
    ax.axvline(6.5, linewidth=.8)
    ax.scatter(np.full(len(maxab), 14.45), np.arange(len(maxab)),
               s=bubble_size(maxab.values), alpha=.7)
    ax.set_xlim(-.5, 15.05)
    ax.set_ylim(len(zw.index)-.5, -1.35)
    ax.set_yticks(np.arange(len(zw.index)))
    ax.set_yticklabels(zw.index, fontsize=7.7)
    ax.set_xticks(np.arange(14))
    ax.set_xticklabels(DEPTH_NAMES + DEPTH_NAMES, rotation=48, ha='right', fontsize=6.8)
    ax.text(3.0, -1.02, 'Winter/mixed', ha='center', va='bottom', fontsize=8.5)
    ax.text(10.0, -1.02, 'Summer/stratified', ha='center', va='bottom', fontsize=8.5)
    ax.text(14.45, -1.02, 'Max.\nmedian', ha='center', va='bottom', fontsize=7.4, linespacing=.9)
    ax.text(-.035, 1.045, panel_letter, transform=ax.transAxes,
            ha='left', va='bottom', fontweight='bold', fontsize=13)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    return im


def main():
    g16, m16 = aggregate('16S')
    g18, m18 = aggregate('18S')
    p16, w16, s16, z16w, z16s = select_and_split(g16, m16)
    p18, w18, s18, z18w, z18s = select_and_split(g18, m18)

    p16.to_csv(OUT/'16S_selected_taxa_pooled.csv')
    p18.to_csv(OUT/'18S_selected_taxa_pooled.csv')
    pd.concat({'winter':w16,'summer':s16}, axis=1).to_csv(OUT/'16S_winter_summer_medians.csv')
    pd.concat({'winter':w18,'summer':s18}, axis=1).to_csv(OUT/'18S_winter_summer_medians.csv')
    pd.concat({'winter':z16w,'summer':z16s}, axis=1).to_csv(OUT/'16S_winter_summer_zscores.csv')
    pd.concat({'winter':z18w,'summer':z18s}, axis=1).to_csv(OUT/'18S_winter_summer_zscores.csv')

    traits = pd.read_csv(TRAIT_FILE)
    ts = trait_summaries(traits)
    ts.to_csv(OUT/'trait_depth_summaries_selected_five_normalized.csv', index=False)

    zall = np.concatenate([z16w.values.ravel(), z16s.values.ravel(),
                           z18w.values.ravel(), z18s.values.ravel()])
    zlim = max(2.0, float(np.nanquantile(np.abs(zall), .985)))
    norm = TwoSlopeNorm(vmin=-zlim, vcenter=0, vmax=zlim)

    fig = plt.figure(figsize=(13.2, 8.4))
    outer = GridSpec(2, 2, figure=fig, height_ratios=[1.06,1.0], hspace=.40, wspace=.34)

    axa = fig.add_subplot(outer[0,0])
    im = plot_heatmap(axa, z16w, z16s, w16, s16, 'a', norm)
    axb = fig.add_subplot(outer[0,1])
    plot_heatmap(axb, z18w, z18s, w18, s18, 'b', norm)

    cax = fig.add_axes([.925,.635,.014,.19])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label('Within-taxon z-score', fontsize=7.5)
    cb.ax.tick_params(labelsize=7)

    lax = fig.add_axes([.902,.495,.082,.105])
    lax.set_xlim(0,1); lax.set_ylim(0,1); lax.axis('off')
    lax.text(.5,.95,'Max. median\nrelative abundance',ha='center',va='top',fontsize=7.4)
    xs = np.array([.18,.50,.82]); vals = np.array([1,10,30],dtype=float)
    lax.scatter(xs, np.repeat(.42,3), s=bubble_size(vals), alpha=.7)
    for x,v in zip(xs,vals):
        lax.text(x,.12,f'{int(v)}%',ha='center',va='center',fontsize=6.8)

    sub = GridSpecFromSubplotSpec(1, 5, subplot_spec=outer[1,:], wspace=.28)
    for i,(trait,title) in enumerate(TRAITS):
        ax = fig.add_subplot(sub[0,i])
        for season in SEASONS:
            z = ts[(ts.trait==trait)&(ts.season==season)].sort_values('depth_m')
            ax.plot(z.median_selected_pct, z.depth_m, linewidth=1.6,
                    linestyle='-' if season=='A_Winter' else '--',
                    label='Winter/mixed' if season=='A_Winter' else 'Summer/stratified')
            ax.fill_betweenx(z.depth_m, z.q25_selected_pct, z.q75_selected_pct, alpha=.10)
        ax.axhline(POOLED_TROPHIC_OPTIMUM_M, linestyle=':', linewidth=.9)
        ax.axhline(POOLED_TAXONOMIC_OPTIMUM_M, linestyle=':', linewidth=.9)
        ax.set_yscale('log'); ax.set_ylim(1700,8)
        ax.set_title(title, fontsize=8.8, pad=5)
        ax.tick_params(labelsize=7)
        if i == 0:
            ax.set_ylabel('Depth (m)', fontsize=8)
            ax.text(-.08,1.04,'c',transform=ax.transAxes,ha='left',va='bottom',fontweight='bold',fontsize=13)
            ax.legend(frameon=False,fontsize=6.8,loc='lower right')
        else:
            ax.tick_params(labelleft=False)
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

    fig.text(.50,.025,'Share of selected trophic representation (%)',ha='center',va='center',fontsize=8)
    fig.subplots_adjust(left=.12,right=.89,top=.96,bottom=.075)
    fig.savefig(OUT/'Figure3_horizontal_publication.pdf', bbox_inches='tight')
    fig.savefig(OUT/'Figure3_horizontal_publication.png', dpi=600, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()
