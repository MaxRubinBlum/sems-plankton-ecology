"""Reproduce Figure 1 analytical panels and export plotted data.

Expected repository inputs (see inputs.tsv): cleaned 16S/18S feature tables and
integrated sample metadata. Run from any directory; paths are repo-relative.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'data'; OUT=ROOT/'results'/'fig01'; OUT.mkdir(parents=True,exist_ok=True)
ENV=['ctd_temperature','ctd_salinity','ctd_oxygen','ctd_fluorescence']


def pcoa_bray(otu, sample_ids):
    ids=[x for x in sample_ids if x in otu.columns]
    x=otu[ids].T.astype(float).to_numpy()
    x=x/np.maximum(x.sum(1,keepdims=True),1)
    d=squareform(pdist(x,metric='braycurtis'))
    n=len(ids); j=np.eye(n)-np.ones((n,n))/n
    b=-0.5*j@(d*d)@j
    eigval,eigvec=np.linalg.eigh(b); ix=np.argsort(eigval)[::-1]
    eigval=eigval[ix]; eigvec=eigvec[:,ix]
    pos=eigval>0
    coords=eigvec[:,pos]*np.sqrt(eigval[pos])
    pct=100*eigval[pos]/eigval[pos].sum()
    return pd.DataFrame(coords[:,:2],index=ids,columns=['PCo1','PCo2']),pct[:2]


def main():
    m16=pd.read_csv(DATA/'16S_samples_CTD_chemistry.csv')
    m18=pd.read_csv(DATA/'18S_samples_CTD_chemistry.csv')
    # Hydrographic PCA is fitted to the 18S-associated environmental observations,
    # matching the final Figure 1 checkpoint; duplicate environmental rows are retained
    # because the panel describes conditions associated with microbial samples.
    env=m18.dropna(subset=ENV).copy()
    z=StandardScaler().fit_transform(env[ENV])
    p=PCA(n_components=4).fit(z); scores=p.transform(z)
    sc=env[['sample-id','depth','ds2','season','station','cruise']].copy()
    sc['PC1']=scores[:,0]; sc['PC2']=scores[:,1]
    sc.to_csv(OUT/'hydrographic_pca_scores.csv',index=False)
    pd.DataFrame(p.components_[:2].T,index=ENV,columns=['PC1','PC2']).to_csv(OUT/'hydrographic_pca_loadings.csv')
    pd.DataFrame({'axis':['PC1','PC2'],'variance_percent':100*p.explained_variance_ratio_[:2]}).to_csv(OUT/'hydrographic_pca_variance.csv',index=False)

    for marker,meta in [('16S',m16),('18S',m18)]:
        otu=pd.read_csv(DATA/('otu_table_prok_clean.csv' if marker=='16S' else 'otu_table_euk_clean.csv')).set_index('OTUID')
        ids=meta['sample-id'].dropna().tolist()
        co,pct=pcoa_bray(otu,ids)
        dat=meta.set_index('sample-id').reindex(co.index)[['depth','ds2','season','station','cruise']].join(co)
        dat.to_csv(OUT/f'{marker}_bray_pcoa_scores.csv')
        pd.DataFrame({'axis':['PCo1','PCo2'],'variance_percent':pct}).to_csv(OUT/f'{marker}_bray_pcoa_variance.csv',index=False)

    # Exact environmental plotting table for the depth/context panel.
    keep=['sample-id','station','station_norm','depth','ds2','season','cruise']+ENV+['dcm_pressure_db']
    pd.concat([m16.assign(marker='16S'),m18.assign(marker='18S')],ignore_index=True)[['marker']+keep].to_csv(OUT/'environmental_plot_data.csv',index=False)

    print('Hydrographic PCA PC1/PC2:',100*p.explained_variance_ratio_[:2])

if __name__=='__main__': main()
