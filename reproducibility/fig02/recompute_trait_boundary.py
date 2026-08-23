"""Recompute the 18S trophic-trait component of Figure 2 from repository data."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
ROOT=Path(__file__).resolve().parents[2]; TRAIT_FILE=ROOT/'data'/'processed'/'traits'/'sample_level_core_functional_traits.csv'; OUT=ROOT/'results'/'fig02'; OUT.mkdir(parents=True,exist_ok=True)
TRAITS=['bona_fide_phototrophy','constitutive_mixotrophy','parasitism','phagotrophy_radiolaria','diplonemid_heterotrophy']; THRESHOLDS=np.arange(10,421,10,dtype=float); MAX_DEPTH=650.; MIN_SPANNING=15; N_BOOT=250; SEED=20260817
def gower(x):
    d=squareform(pdist(x,'euclidean')); n=len(d); j=np.eye(n)-np.ones((n,n))/n; return -.5*j@(d*d)@j
def dummy(v): return pd.get_dummies(pd.Series(v).astype(str),drop_first=True,dtype=float).to_numpy()
def pr2(g,grp,z):
    n=len(grp); x0=np.c_[np.ones(n),z]; h0=x0@np.linalg.pinv(x0); ss0=np.trace((np.eye(n)-h0)@g); x1=np.c_[x0,np.asarray(grp,float)]; h1=x1@np.linalg.pinv(x1); ss1=np.trace((np.eye(n)-h1)@g); return float((ss0-ss1)/ss0)
def profile(m): return m.cruise.astype(str)+'|'+m.station.astype(str)
def prepare():
    m=pd.read_csv(TRAIT_FILE).dropna(subset=['depth','station','cruise']).copy(); m=m[m.depth<=MAX_DEPTH].reset_index(drop=True); x=m[TRAITS].fillna(0).to_numpy(float); x=x/np.maximum(x.sum(1,keepdims=True),1e-12); return np.sqrt(x),m
def scan(x,m,season='pooled'):
    if season!='pooled': keep=(m.season==season).to_numpy(); x=x[keep]; m=m.loc[keep].reset_index(drop=True)
    gg=gower(x); z=dummy(m.station); pid=profile(m); rows=[]
    for t in THRESHOLDS:
        sh=(m.depth.to_numpy()<=t); q=pd.DataFrame({'p':pid,'s':sh}); span=int(q.groupby('p').s.agg(lambda a:a.any() and (~a).any()).sum())
        if span<MIN_SPANNING: continue
        rows.append([t,pr2(gg,(~sh).astype(int),z),int(sh.sum()),int((~sh).sum()),span,'18S_function',season])
    return pd.DataFrame(rows,columns=['threshold_m','partial_R2','n_shallow','n_deep','n_profiles_spanning','dataset','season'])
def boot(x,m,season='pooled'):
    if season!='pooled': keep=(m.season==season).to_numpy(); x=x[keep]; m=m.loc[keep].reset_index(drop=True)
    rng=np.random.default_rng(SEED); pid=profile(m).to_numpy(); u=np.unique(pid); rows=[]
    for b in range(N_BOOT):
        chosen=rng.choice(u,len(u),replace=True); xx=[]; mm=[]
        for k,p in enumerate(chosen):
            ix=np.where(pid==p)[0]; xx.append(x[ix]); q=m.iloc[ix].copy(); q['_bp']=str(k); mm.append(q)
        xb=np.vstack(xx); mb=pd.concat(mm,ignore_index=True); gg=gower(xb); z=dummy(mb.station); cand=[]
        for t in THRESHOLDS:
            sh=(mb.depth.to_numpy()<=t); q=pd.DataFrame({'p':mb._bp,'s':sh})
            if int(q.groupby('p').s.agg(lambda a:a.any() and (~a).any()).sum())<MIN_SPANNING: continue
            cand.append((t,pr2(gg,(~sh).astype(int),z)))
        if cand:
            t,r=max(cand,key=lambda q:q[1]); rows.append([b,t,r,'18S_function',season])
    return pd.DataFrame(rows,columns=['bootstrap','threshold_m','max_partial_R2','dataset','season'])
if __name__=='__main__':
    x,m=prepare(); scans=[]; boots=[]
    for s in ['pooled','A_Winter','Summer']: scans.append(scan(x,m,s)); boots.append(boot(x,m,s))
    pd.concat(scans,ignore_index=True).to_csv(OUT/'trait_boundary_scans_recomputed.csv',index=False); pd.concat(boots,ignore_index=True).to_csv(OUT/'trait_boundary_bootstrap_recomputed.csv',index=False); print('Manuscript validation targets: pooled=120 m, winter=120 m, summer=120 m.')
