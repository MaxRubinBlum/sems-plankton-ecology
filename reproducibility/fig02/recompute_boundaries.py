"""Recompute the Figure 2 objective upper-boundary analysis from repository data."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

ROOT=Path(__file__).resolve().parents[2]
RAW=ROOT/'data'/'raw'; META=ROOT/'data'/'processed'/'metadata'; OUT=ROOT/'results'/'fig02'; OUT.mkdir(parents=True,exist_ok=True)
SEED=20260817; THRESHOLDS=np.arange(10,421,10,dtype=float); MAX_DEPTH=650.0; MIN_SPANNING_PROFILES=15; N_BOOT=250; N_PERM=499

def hellinger(otu,ids):
    x=otu[ids].T.astype(float).to_numpy(); x=x/np.maximum(x.sum(1,keepdims=True),1); return np.sqrt(x)
def gower_from_euclidean(x):
    d=squareform(pdist(x,metric='euclidean')); n=len(d); j=np.eye(n)-np.ones((n,n))/n; return -.5*j@(d*d)@j
def dummy(vals): return pd.get_dummies(pd.Series(vals).astype(str),drop_first=True,dtype=float).to_numpy()
def partial_r2(gower,group,nuisance):
    n=len(group); x0=np.c_[np.ones(n),nuisance] if nuisance.size else np.ones((n,1)); h0=x0@np.linalg.pinv(x0); ss0=np.trace((np.eye(n)-h0)@gower)
    x1=np.c_[x0,np.asarray(group,float)]; h1=x1@np.linalg.pinv(x1); ss1=np.trace((np.eye(n)-h1)@gower); return float((ss0-ss1)/ss0)
def prepare(marker):
    otu_path=RAW/marker/('otu_table_prok_clean.csv.gz' if marker=='16S' else 'otu_table_euk_clean.csv.gz')
    otu=pd.read_csv(otu_path).set_index('OTUID'); meta=pd.read_csv(META/f'{marker}_samples_CTD_chemistry.csv')
    meta=meta[meta.depth.notna() & (meta.depth<=MAX_DEPTH)].copy(); meta=meta[meta['sample-id'].isin(otu.columns)].reset_index(drop=True)
    ids=meta['sample-id'].tolist(); return hellinger(otu,ids),meta
def profile_id(meta): return meta.cruise.astype(str)+'|'+meta.station.astype(str)
def scan_matrix(x,meta,dataset,season=None):
    if season is not None:
        keep=(meta.season==season).to_numpy(); x=x[keep]; meta=meta.loc[keep].reset_index(drop=True)
    gower=gower_from_euclidean(x); nuisance=dummy(meta.station); pid=profile_id(meta); rows=[]
    for t in THRESHOLDS:
        shallow=(meta.depth.to_numpy()<=t); tmp=pd.DataFrame({'pid':pid,'shallow':shallow}); span=int(tmp.groupby('pid').shallow.agg(lambda z:z.any() and (~z).any()).sum())
        if span<MIN_SPANNING_PROFILES or shallow.sum()<2 or (~shallow).sum()<2: continue
        rows.append([t,partial_r2(gower,(~shallow).astype(int),nuisance),int(shallow.sum()),int((~shallow).sum()),span,dataset,season or 'pooled'])
    return pd.DataFrame(rows,columns=['threshold_m','partial_R2','n_shallow','n_deep','n_profiles_spanning','dataset','season'])
def bootstrap_profiles(x,meta,dataset,season=None,rng=None):
    rng=np.random.default_rng(SEED if rng is None else rng)
    if season is not None:
        keep=(meta.season==season).to_numpy(); x=x[keep]; meta=meta.loc[keep].reset_index(drop=True)
    pid=profile_id(meta).to_numpy(); unique=np.unique(pid); out=[]
    for b in range(N_BOOT):
        chosen=rng.choice(unique,size=len(unique),replace=True); xx=[]; mm=[]
        for k,p in enumerate(chosen):
            ix=np.where(pid==p)[0]; xx.append(x[ix]); z=meta.iloc[ix].copy(); z['_boot_profile']=f'{k}:{p}'; mm.append(z)
        xb=np.vstack(xx); mb=pd.concat(mm,ignore_index=True); gb=gower_from_euclidean(xb); nuisance=dummy(mb.station); rows=[]
        for t in THRESHOLDS:
            shallow=(mb.depth.to_numpy()<=t); tmp=pd.DataFrame({'pid':mb['_boot_profile'],'shallow':shallow}); span=int(tmp.groupby('pid').shallow.agg(lambda z:z.any() and (~z).any()).sum())
            if span<MIN_SPANNING_PROFILES: continue
            rows.append((t,partial_r2(gb,(~shallow).astype(int),nuisance)))
        if rows:
            t,r=max(rows,key=lambda q:q[1]); out.append([b,t,r,dataset,season or 'pooled'])
    return pd.DataFrame(out,columns=['bootstrap','threshold_m','max_partial_R2','dataset','season'])
def maxstat_permutation(x,meta,observed,dataset,rng=None):
    rng=np.random.default_rng(SEED+1 if rng is None else rng); gower=gower_from_euclidean(x); nuisance=dummy(meta.station); null=[]; valid=observed.threshold_m.to_numpy(); depth=meta.depth.to_numpy()
    for _ in range(N_PERM):
        mx=-np.inf
        for t in valid:
            grp=(depth>t).astype(int).copy()
            for st in meta.station.astype(str).unique():
                ix=np.where(meta.station.astype(str).to_numpy()==st)[0]; grp[ix]=rng.permutation(grp[ix])
            mx=max(mx,partial_r2(gower,grp,nuisance))
        null.append(mx)
    obs=float(observed.partial_R2.max()); p=(1+np.sum(np.asarray(null)>=obs))/(N_PERM+1)
    return pd.DataFrame({'dataset':[dataset],'observed_best_threshold_m':[float(observed.loc[observed.partial_R2.idxmax(),'threshold_m'])],'observed_max_partial_R2':[obs],'permutation_p_maxstat':[p],'null_median_max_R2':[np.median(null)],'null_95pct_max_R2':[np.quantile(null,.95)]}),pd.Series(null,name='null_max_partial_R2')
def main():
    scans=[]; boots=[]; tests=[]; nulls=[]
    for marker in ['16S','18S']:
        x,m=prepare(marker); pooled=scan_matrix(x,m,marker); scans.append(pooled); boots.append(bootstrap_profiles(x,m,marker)); test,null=maxstat_permutation(x,m,pooled,marker); tests.append(test); nulls.append(pd.DataFrame({'dataset':marker,'null_max_partial_R2':null}))
        for season in ['A_Winter','Summer']:
            scans.append(scan_matrix(x,m,marker,season)); boots.append(bootstrap_profiles(x,m,marker,season))
    pd.concat(scans,ignore_index=True).to_csv(OUT/'boundary_scans_recomputed.csv',index=False); pd.concat(boots,ignore_index=True).to_csv(OUT/'boundary_bootstrap_recomputed.csv',index=False); pd.concat(tests,ignore_index=True).to_csv(OUT/'boundary_maxstat_recomputed.csv',index=False); pd.concat(nulls,ignore_index=True).to_csv(OUT/'boundary_maxstat_nulls_recomputed.csv',index=False)
    print('Expected manuscript peaks for validation: pooled 16S=170 m, pooled 18S=220 m; winter 16S=220 m, 18S=220 m; summer 16S=140 m, 18S=160 m.')
if __name__=='__main__': main()
