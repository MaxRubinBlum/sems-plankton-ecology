#!/usr/bin/env python3
import argparse,numpy as np,pandas as pd
from scipy.stats import rankdata,t
from statsmodels.stats.multitest import multipletests

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--otu',required=True);ap.add_argument('--taxonomy',required=True);ap.add_argument('--mapping',required=True);ap.add_argument('--features-parent',required=True);ap.add_argument('--metadata-parent',required=True);ap.add_argument('--candidates',required=True);ap.add_argument('--out-prefix',required=True);a=ap.parse_args()
 otu=pd.read_csv(a.otu).set_index('OTUID');tax=pd.read_csv(a.taxonomy).set_index('OTUID').fillna('');mp=pd.read_csv(a.mapping,sep='\t').set_index('feature_id');Fagg=pd.read_csv(a.features_parent,sep='\t',index_col=0);M=pd.read_csv(a.metadata_parent,sep='\t',index_col=0).loc[Fagg.index];cand=pd.read_csv(a.candidates)
 ranks=['Species','Genus','Family','Order','Class','Phylum']
 def norm(x):
  x=str(x).strip()
  if x.lower() in {'','nan','unclassified','uncultured','unknown','none'}:return ''
  if '__' in x:x=x.split('__',1)[1]
  return x.replace('_',' ').strip()
 def label(r):
  for c in ranks:
   v=norm(r.get(c,''))
   if v:return v
  return 'Unclassified Eukaryota'
 def group(r):
  lin=';'.join(norm(r.get(c,'')) for c in ['Phylum','Class','Order','Family','Genus','Species']).lower();syn=('syndin' in lin or 'dino-group-i' in lin or 'dino-group-ii' in lin)
  if 'metazoa' in lin:return 'Metazoa'
  if syn:return 'Syndiniales'
  if 'hematodinium' in lin:return 'Hematodinium'
  return 'Other eukaryote'
 ann=pd.DataFrame(index=tax.index);ann['taxon_label']=tax.apply(label,axis=1);ann['ecological_group']=tax.apply(group,axis=1);ann['agg_key']=ann.ecological_group+' | '+ann.taxon_label;fidkey={fid:f'{r.ecological_group} | {r.taxon_label}' for fid,r in mp.iterrows()}
 raw=otu[Fagg.index].astype(float);prev=(raw>0).mean(axis=1);tot=raw.sum(axis=1);eligible=(prev>=.02)&(tot>=50)
 def trim(ids):
  if not ids:return []
  s=tot.loc[ids].sort_values(ascending=False);cs=s.cumsum()/s.sum();k=list(s.index[cs<=.90]);k=k or [s.index[0]]
  if len(k)<len(s):k.append(s.index[len(k)])
  return k[:20]
 specs=[];allids=set()
 for _,r in cand.iterrows():
  p,q=r.edge_key.split('|');pids=trim(ann.index[(ann.agg_key==fidkey[p])&eligible].tolist());qids=trim(ann.index[(ann.agg_key==fidkey[q])&eligible].tolist())
  if pids and qids:specs.append((r.edge_key,r.parasite,r.partner,r.tier,pids,qids));allids.update(pids);allids.update(qids)
 A=M.astype(float);A=A.loc[:,A.nunique()>1];X=np.column_stack([np.ones(len(A)),A.values]);ids=list(allids);Y=np.log1p(raw.loc[ids,Fagg.index].T.values);R=Y-X@np.linalg.lstsq(X,Y,rcond=None)[0];pos={x:i for i,x in enumerate(ids)};Z=np.column_stack([rankdata(R[:,j]) for j in range(R.shape[1])]);Z=(Z-Z.mean(0))/Z.std(0,ddof=1);n=len(A);rows=[]
 for edge,pn,qn,tier,pids,qids in specs:
  for p in pids:
   for q in qids:
    co=int(((raw.loc[p,Fagg.index]>0)&(raw.loc[q,Fagg.index]>0)).sum())
    if co<20:continue
    rho=float(np.dot(Z[:,pos[p]],Z[:,pos[q]])/(n-1));tt=rho*np.sqrt((n-2)/max(1e-12,1-rho*rho));pv=float(2*t.sf(abs(tt),df=n-2));rows.append([edge,pn,qn,tier,p,q,co,rho,pv])
 R=pd.DataFrame(rows,columns=['aggregate_edge','parasite_taxon','partner_taxon','aggregate_tier','parasite_ASV','partner_ASV','n_copresent','residual_spearman_rho','p']);R['q_global']=multipletests(R.p,method='fdr_bh')[1];R['supported']=(R.q_global<.05)&(R.residual_spearman_rho>=.20);R.to_csv(a.out_prefix+'_pair_tests.csv',index=False)
 sm=R.groupby('aggregate_edge').agg(ASV_pairs_tested=('supported','size'),ASV_pairs_FDR_supported=('supported','sum')).reset_index();sm.to_csv(a.out_prefix+'_summary.csv',index=False)
if __name__=='__main__':main()
