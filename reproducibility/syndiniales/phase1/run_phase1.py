#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import f_oneway, spearmanr
from numpy.linalg import pinv, eigh

SEED=20260826
THRESHOLDS=np.arange(50,401,10)

def design(df, season=True, station=True):
    xs=[np.ones((len(df),1))]
    if station:
        xs.append(pd.get_dummies(df['station'],drop_first=True,dtype=float).values)
    if season and df['season'].nunique()>1:
        xs.append(pd.get_dummies(df['season'],drop_first=True,dtype=float).values)
    return np.column_stack(xs)

def projection(G,df,season=True,station=True):
    X=design(df,season,station)
    Q=X@pinv(X.T@X)@X.T
    return Q,np.trace(G)-np.trace(Q@G)

def partial_r2(G,z,Q,denom):
    z=z.astype(float); zr=z-Q@z; zz=float(zr@zr)
    if zz<1e-12:return np.nan
    return max(0,float(zr@(G@zr))/zz/denom)

def clade_label(r):
    for c in ['Family','Genus','Order']:
        v=str(r.get(c,'')).strip()
        if v and v.lower() not in {'nan','unclassified'}: return v.replace('_',' ')
    return 'Unresolved Syndiniales'

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--otu',required=True); ap.add_argument('--taxonomy',required=True); ap.add_argument('--metadata',required=True); ap.add_argument('--outdir',required=True)
    a=ap.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True); rng=np.random.default_rng(SEED)
    otu=pd.read_csv(a.otu).set_index('OTUID'); tax=pd.read_csv(a.taxonomy).set_index('OTUID').fillna(''); meta=pd.read_csv(a.metadata).set_index('sample-id')
    meta=meta[meta['select'].astype(str).eq('a')].copy(); meta['depth']=pd.to_numeric(meta['depth'],errors='coerce'); samples=[s for s in meta.index if s in otu.columns and pd.notna(meta.loc[s,'depth'])]; meta=meta.loc[samples]
    syn_ids=tax.index[tax['Class'].astype(str).str.contains('Syndiniales',case=False,na=False)]; syn_ids=[x for x in syn_ids if x in otu.index]; syn=otu.loc[syn_ids,samples].astype(float); st=syn.sum(axis=0); samples=st[st>0].index.tolist(); meta=meta.loc[samples]; syn=syn[samples]; st=st[samples]; total=otu[samples].sum(axis=0).astype(float)
    P=(syn/st).T; P=P.loc[:,P.sum(axis=0)>0]; D=squareform(pdist(P.values,metric='braycurtis')); n=len(P); J=np.eye(n)-np.ones((n,n))/n; G=-0.5*J@(D**2)@J; Q,den=projection(G,meta,True,True)
    scan=[]
    for t in THRESHOLDS:
        z=meta.depth.values>t
        if min((~z).sum(),z.sum())>=25:
            prof=meta.assign(side=np.where(z,'deep','upper')).groupby(['station','season'])['side'].nunique(); scan.append([t,partial_r2(G,z,Q,den),(~z).sum(),z.sum(),int((prof==2).sum())])
    scan=pd.DataFrame(scan,columns=['threshold_m','partial_R2','n_shallow','n_deep','n_profiles_spanning']); bp=float(scan.loc[scan.partial_R2.idxmax(),'threshold_m']); obs=float(scan.partial_R2.max()); scan.to_csv(out/'syndiniales_objective_boundary_scan.csv',index=False)
    prof_indices=[np.array(v,dtype=int) for v in meta.groupby(['station','season']).indices.values()]; depth=meta.depth.values.copy(); null=[]
    for _ in range(199):
        dp=depth.copy()
        for ii in prof_indices: dp[ii]=rng.permutation(dp[ii])
        null.append(max(partial_r2(G,dp>t,Q,den) for t in scan.threshold_m))
    pmax=(1+np.sum(np.array(null)>=obs))/(len(null)+1); keys=list(meta.groupby(['station','season']).indices.values()); boot=[]
    for _ in range(150):
        picks=rng.choice(len(keys),len(keys),replace=True); idx=[]; frames=[]
        for j,k in enumerate(picks):
            ii=np.array(keys[k],dtype=int); idx.extend(ii.tolist()); dd=meta.iloc[ii].copy(); dd['station']=dd['station'].astype(str)+f'#{j}'; frames.append(dd)
        bm=pd.concat(frames); BG=G[np.ix_(idx,idx)]; BQ,Bden=projection(BG,bm,True,True); vals=[]
        for t in scan.threshold_m:
            z=bm.depth.values>t
            if min((~z).sum(),z.sum())>=20: vals.append((t,partial_r2(BG,z,BQ,Bden)))
        boot.append(max(vals,key=lambda x:x[1])[0])
    boot=np.asarray(boot); pd.DataFrame([{'n_samples':len(meta),'n_syndiniales_ASVs':len(syn_ids),'observed_best_threshold_m':bp,'observed_max_partial_R2':obs,'maxstat_permutation_p':pmax,'n_permutations':len(null),'bootstrap_n':len(boot),'bootstrap_median_m':np.median(boot),'bootstrap_CI2.5_m':np.percentile(boot,2.5),'bootstrap_CI97.5_m':np.percentile(boot,97.5),'bootstrap_mode_m':pd.Series(boot).mode().iloc[0]}]).to_csv(out/'syndiniales_boundary_summary.csv',index=False)
    ssall=[]; ssbest=[]
    for sea,ii in meta.groupby('season').indices.items():
        ii=np.array(ii,dtype=int); sm=meta.iloc[ii]; SG=G[np.ix_(ii,ii)]; SQ,Sden=projection(SG,sm,False,True); rows=[]
        for t in THRESHOLDS:
            z=sm.depth.values>t
            if min((~z).sum(),z.sum())>=15: rows.append([sea,t,partial_r2(SG,z,SQ,Sden),(~z).sum(),z.sum()])
        sdf=pd.DataFrame(rows,columns=['season','threshold_m','partial_R2','n_shallow','n_deep']); ssall.append(sdf); ssbest.append(sdf.loc[sdf.partial_R2.idxmax()].to_dict())
    pd.concat(ssall).to_csv(out/'syndiniales_season_boundary_scans.csv',index=False); pd.DataFrame(ssbest).to_csv(out/'syndiniales_season_best_boundaries.csv',index=False)
    ev,ec=eigh(G); order=np.argsort(ev)[::-1]; ev=ev[order]; ec=ec[:,order]; pos=ev>1e-10; ev=ev[pos]; ec=ec[:,pos]; coords=ec*np.sqrt(ev); var=ev/ev.sum(); pd.DataFrame({'sample_id':meta.index,'PCo1':coords[:,0],'PCo2':coords[:,1],'depth':meta.depth.values,'season':meta.season.values,'station':meta.station.values}).to_csv(out/'syndiniales_PCoA_scores.csv',index=False); pd.DataFrame({'axis':np.arange(1,len(var)+1),'variance_fraction':var}).to_csv(out/'syndiniales_PCoA_variance.csv',index=False)
    zone=np.where(meta.depth.values<=bp,'upper','deep'); dc=np.zeros(n)
    for g in np.unique(zone):
        ii=np.where(zone==g)[0]; cen=coords[ii].mean(axis=0); dc[ii]=np.sqrt(((coords[ii]-cen)**2).sum(axis=1))
    F=f_oneway(dc[zone=='upper'],dc[zone=='deep']).statistic; pf=[]
    for _ in range(499):
        zp=zone.copy()
        for ii in prof_indices: zp[ii]=rng.permutation(zp[ii])
        dd=np.zeros(n)
        for g in np.unique(zp):
            jj=np.where(zp==g)[0]; cen=coords[jj].mean(axis=0); dd[jj]=np.sqrt(((coords[jj]-cen)**2).sum(axis=1))
        pf.append(f_oneway(dd[zp=='upper'],dd[zp=='deep']).statistic)
    pp=(1+np.sum(np.array(pf)>=F))/(len(pf)+1); pd.DataFrame([{'breakpoint_m':bp,'F':F,'permutation_p':pp,'mean_distance_upper':dc[zone=='upper'].mean(),'mean_distance_deep':dc[zone=='deep'].mean()}]).to_csv(out/'syndiniales_PERMDISP.csv',index=False)
    eps=1e-15; sh=-(P*np.log(P+eps)).sum(axis=1); rich=(syn.T>0).sum(axis=1); m=meta[['station','season','depth']].copy(); m['syndiniales_fraction_18S']=st/total; m['syndiniales_richness']=rich; m['syndiniales_shannon']=sh; m.to_csv(out/'syndiniales_sample_metrics.csv')
    labs=tax.loc[syn_ids].apply(clade_label,axis=1); C=syn.copy(); C['clade']=labs.reindex(C.index).values; C=C.groupby('clade').sum().T; CP=C.div(C.sum(axis=1),axis=0); rows=[]
    for c in CP:
        r,p=spearmanr(meta.depth,CP[c]); u=CP.loc[meta.depth<=bp,c]; d=CP.loc[meta.depth>bp,c]; rows.append([c,(CP[c]>0).mean(),CP[c].mean(),u.median(),d.median(),r,p])
    pd.DataFrame(rows,columns=['clade','prevalence','mean_fraction_within_Syndiniales','median_upper','median_deep','spearman_depth_rho','spearman_depth_p']).sort_values('mean_fraction_within_Syndiniales',ascending=False).to_csv(out/'syndiniales_clade_vertical_statistics.csv',index=False)

if __name__=='__main__': main()
