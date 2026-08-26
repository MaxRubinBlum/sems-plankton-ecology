#!/usr/bin/env python3
import argparse,numpy as np,pandas as pd
from scipy.spatial.distance import pdist,squareform
from numpy.linalg import eigh,pinv,svd

def pcoa(P,k=10):
 D=squareform(pdist(P.values,metric='braycurtis'));n=len(P);J=np.eye(n)-np.ones((n,n))/n;G=-.5*J@(D**2)@J;v,U=eigh(G);o=np.argsort(v)[::-1];v=v[o];U=U[:,o];ok=v>1e-10;v=v[ok];U=U[:,ok];return (U*np.sqrt(v))[:,:min(k,len(v))]
def proc(X,Y):
 k=min(X.shape[1],Y.shape[1]);X=X[:,:k]-X[:,:k].mean(0);Y=Y[:,:k]-Y[:,:k].mean(0);X/=np.linalg.norm(X);Y/=np.linalg.norm(Y);return float(svd(X.T@Y,compute_uv=False).sum())
def residualize(X,m,env):
 C=m[['depth']+env].astype(float).copy()
 for c in C:
  sd=C[c].std(ddof=0);C[c]=(C[c]-C[c].mean())/sd if sd>0 else 0
 dum=pd.concat([pd.get_dummies(m.season,drop_first=True,dtype=float),pd.get_dummies(m.station,drop_first=True,dtype=float)],axis=1);A=np.column_stack([np.ones(len(m)),C.values,dum.values]);H=A@pinv(A.T@A)@A.T;return X-H@X
def boundary(ids,otu,meta,samples):
 C=otu.loc[ids,samples].astype(float);C=C.loc[C.sum(axis=1)>0];P=(C/C.sum(axis=0)).T;D=squareform(pdist(P.values,metric='braycurtis'));n=len(P);J=np.eye(n)-np.ones((n,n))/n;G=-.5*J@(D**2)@J;X=np.column_stack([np.ones(n),pd.get_dummies(meta.station,drop_first=True,dtype=float).values,pd.get_dummies(meta.season,drop_first=True,dtype=float).values]);H=X@pinv(X.T@X)@X.T;den=np.trace(G)-np.trace(H@G);rows=[]
 for t in np.arange(50,401,10):
  z=(meta.depth.values>t).astype(float)
  if min((z==0).sum(),(z==1).sum())<25:continue
  zr=z-H@z;rows.append([t,float(zr@(G@zr)/(zr@zr))/den])
 s=pd.DataFrame(rows,columns=['threshold_m','partial_R2']);return s,s.loc[s.partial_R2.idxmax()]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--otu',required=True);ap.add_argument('--taxonomy',required=True);ap.add_argument('--metadata',required=True);ap.add_argument('--ctd',required=True);ap.add_argument('--outdir',required=True);ap.add_argument('--seed',type=int,default=20260826);ap.add_argument('--permutations',type=int,default=999);a=ap.parse_args();from pathlib import Path;out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
 otu=pd.read_csv(a.otu).set_index('OTUID');tax=pd.read_csv(a.taxonomy).set_index('OTUID').fillna('');m0=pd.read_csv(a.metadata).set_index('sample-id');m0=m0[m0['select']=='a'];ctd=pd.read_csv(a.ctd).set_index('sample-id');env=['ctd_temperature','ctd_salinity','ctd_oxygen','ctd_fluorescence'];m=m0.join(ctd[env],how='left');m['depth']=pd.to_numeric(m.depth,errors='coerce');samples=[s for s in m.index if s in otu.columns and pd.notna(m.loc[s,'depth'])];m=m.loc[samples]
 lin=tax[['Phylum','Class','Order','Family','Genus','Species']].astype(str).agg(';'.join,axis=1).str.lower();syn=set(tax.index[tax.Class.str.contains('Syndiniales',case=False,na=False)]);met=set(tax.index[lin.str.contains('metazoa')]);syn=[x for x in syn if x in otu.index];non=[x for x in otu.index if x not in set(syn) and x not in met]
 def comp(ids):
  C=otu.loc[ids,samples].astype(float);C=C.loc[C.sum(axis=1)>0];return (C/C.sum(axis=0)).T
 Ps=comp(syn);Pn=comp(non);common=Ps.index.intersection(Pn.index);Ps=Ps.loc[common];Pn=Pn.loc[common];m=m.loc[common];Xs=pcoa(Ps);Xn=pcoa(Pn);rng=np.random.default_rng(a.seed);raw=proc(Xs,Xn);groups=[np.array(v,dtype=int) for v in m.groupby(['station','season']).indices.values()]
 def rp():
  p=np.arange(len(m))
  for g in groups:p[g]=rng.permutation(g)
  return p
 nul=np.array([proc(Xs,Xn[rp()]) for _ in range(a.permutations)]);rows=[['global raw',len(m),raw,(1+(nul>=raw).sum())/(a.permutations+1)]];ok=m[env].notna().all(axis=1);mc=m.loc[ok];Xsc=Xs[ok.values];Xnc=Xn[ok.values];rs=residualize(Xsc,mc,env);rn=residualize(Xnc,mc,env);rr=proc(rs,rn);nul=np.array([proc(rs,rn[rng.permutation(len(mc))]) for _ in range(a.permutations)]);rows.append(['global environment-adjusted',len(mc),rr,(1+(nul>=rr).sum())/(a.permutations+1)])
 for label,mask in [('<=220 m',mc.depth<=220),('>220 m',mc.depth>220)]:
  idx=np.where(mask.values)[0];mm=mc.loc[mask];u=residualize(Xsc[idx],mm,env);v=residualize(Xnc[idx],mm,env);q=proc(u,v);nul=np.array([proc(u,v[rng.permutation(len(mm))]) for _ in range(a.permutations)]);rows.append([label+' environment-adjusted',len(mm),q,(1+(nul>=q).sum())/(a.permutations+1)])
 pd.DataFrame(rows,columns=['analysis','n','procrustes_r','permutation_p']).to_csv(out/'community_concordance.csv',index=False);hs=[]
 for hab in ['A_surface','B_nearsurface','C_DCM','D_below_DCM','F_300-600','G_below_600','H_near_bottom']:
  mask=(mc.ds3==hab).values
  if mask.sum()<15:continue
  mm=mc.loc[mask];u=residualize(Xsc[mask],mm,env);v=residualize(Xnc[mask],mm,env);q=proc(u,v);nul=np.array([proc(u,v[rng.permutation(len(mm))]) for _ in range(a.permutations)]);hs.append([hab,len(mm),q,(1+(nul>=q).sum())/(a.permutations+1)])
 pd.DataFrame(hs,columns=['habitat','n','adjusted_procrustes_r','permutation_p']).to_csv(out/'within_habitat_concordance.csv',index=False);micro=[x for x in otu.index if x not in met];nonsyn=[x for x in micro if x not in set(syn)];s1,b1=boundary(micro,otu,m,samples);s2,b2=boundary(nonsyn,otu,m,samples);s1.to_csv(out/'boundary_scan_with_syndiniales.csv',index=False);s2.to_csv(out/'boundary_scan_without_syndiniales.csv',index=False);pd.DataFrame([['with Syndiniales',b1.threshold_m,b1.partial_R2],['without Syndiniales',b2.threshold_m,b2.partial_R2]],columns=['dataset','best_threshold_m','partial_R2']).to_csv(out/'boundary_independence_summary.csv',index=False)
if __name__=='__main__':main()
