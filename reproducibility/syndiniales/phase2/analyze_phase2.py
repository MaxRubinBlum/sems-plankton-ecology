#!/usr/bin/env python3
import argparse,re
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr,fisher_exact
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm

def read_edges(path):
 rows=[]
 for line in Path(path).read_text(errors='ignore').splitlines():
  if not line.strip() or line.startswith('#'): continue
  p=line.split('\t')
  if len(p)<3: p=re.split(r'\s+',line.strip())
  if len(p)>=3:
   try: rows.append((p[0],p[1],float(p[2])))
   except: pass
 return pd.DataFrame(rows,columns=['a','b','weight'])

def classify(x):
 z=str(x).lower(); rules=[('Dinoflagellates',['dinoflag','gymnodinium','prorocentrum','tripos','lepidodinium','dinophy']),('Radiolaria / Acantharea',['rad-a','rad-b','rad-c','radiolaria','acanth','polycyst','collodaria','nassell','spumell']),('MAST',['mast-','mast ']),('Diplonemids',['diplonem','eupelagonem']),('Telonemids',['telonem']),('Ciliates',['cilioph','spirotrich','strombid','tintinn']),('Chlorarachniophytes',['chlorarach']),('Haptophytes',['haptophy','phaeocyst','prymnes']),('Picozoa',['picozo','picomonas']),('Choanoflagellates',['choanoflag']),('Pelagophytes',['pelagomon','pelagophy']),('Diatoms',['diatom','bacillario']),('Chlorophytes',['chlorophy','micromonas','bathycoccus','ostreococcus']),('Fungi',['fung'])]
 for g,ks in rules:
  if any(k in z for k in ks): return g
 return 'Other protists'

def extract(edgefile,mapping):
 E=read_edges(edgefile); mp=mapping.set_index('feature_id'); rows=[]
 for _,e in E.iterrows():
  if e.a not in mp.index or e.b not in mp.index: continue
  A,B=mp.loc[e.a],mp.loc[e.b]
  if bool(A.parasite_flag) and B.ecological_group=='Other eukaryote': p,q=e.a,e.b
  elif bool(B.parasite_flag) and A.ecological_group=='Other eukaryote': p,q=e.b,e.a
  else: continue
  rows.append([p,q,e.weight,mp.loc[p,'taxon_label'],mp.loc[q,'taxon_label']])
 R=pd.DataFrame(rows,columns=['parasite_id','partner_id','weight','parasite','partner']); R['partner_group']=R.partner.map(classify); R['edge_key']=R.parasite_id+'|'+R.partner_id
 return E,R

def residual_matrix(F,M,needed):
 ids=[x for x in needed if x in F.columns]; A=M.select_dtypes(include=[np.number]).astype(float); A=A.loc[:,A.nunique()>1]; A=(A-A.mean())/A.std(ddof=0); X=np.column_stack([np.ones(len(A)),A.values]); Y=np.log1p(F[ids].values.astype(float)); B=np.linalg.lstsq(X,Y,rcond=None)[0]; return pd.DataFrame(Y-X@B,index=F.index,columns=ids)

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--input-dir',required=True); ap.add_argument('--upper-edges',required=True); ap.add_argument('--deep-edges',required=True); ap.add_argument('--outdir',required=True); a=ap.parse_args(); inp=Path(a.input_dir); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
 FU=pd.read_csv(inp/'features_upper_le_220m.tsv',sep='\t',index_col=0); FD=pd.read_csv(inp/'features_deep_gt_220m.tsv',sep='\t',index_col=0); MU=pd.read_csv(inp/'metadata_upper_le_220m.tsv',sep='\t',index_col=0); MD=pd.read_csv(inp/'metadata_deep_gt_220m.tsv',sep='\t',index_col=0); mapU=pd.read_csv(inp/'mapping_upper_le_220m.tsv',sep='\t'); mapD=pd.read_csv(inp/'mapping_deep_gt_220m.tsv',sep='\t')
 _,U=extract(a.upper_edges,mapU); _,D=extract(a.deep_edges,mapD); U.to_csv(out/'upper_parasite_protist_edges.csv',index=False); D.to_csv(out/'deep_parasite_protist_edges.csv',index=False)
 gr=[]
 for reg,R,mp in [('upper',U,mapU),('deep',D,mapD)]:
  pool=mp[(~mp.parasite_flag.astype(bool))&(mp.ecological_group=='Other eukaryote')].copy(); pool['partner_group']=pool.taxon_label.map(classify); np_=R.parasite_id.nunique()
  for g in sorted(set(pool.partner_group)|set(R.partner_group)):
   e=int((R.partner_group==g).sum()); av=int((pool.partner_group==g).sum()); gr.append([reg,g,e,av,np_,e/(np_*av) if np_*av else np.nan])
 pd.DataFrame(gr,columns=['regime','group','edges','available','nparas','density']).to_csv(out/'guild_normalized.csv',index=False)
 union=pd.concat([U,D]).drop_duplicates('edge_key')[['edge_key','parasite_id','partner_id','parasite','partner','partner_group']]; needed=set(union.parasite_id)|set(union.partner_id); RU=residual_matrix(FU,MU,needed); RD=residual_matrix(FD,MD,needed); rows=[]
 for _,e in union.iterrows():
  row=e.to_dict()
  for reg,F,Rm in [('upper',FU,RU),('deep',FD,RD)]:
   p,q=e.parasite_id,e.partner_id; ok=p in F.columns and q in F.columns; row[f'{reg}_both_features']=ok
   if ok:
    pres=(F[p]>0)&(F[q]>0); rr,pp=spearmanr(Rm[p],Rm[q]); row[f'{reg}_n_copresent']=int(pres.sum()); row[f'{reg}_resid_rho']=rr; row[f'{reg}_resid_p']=pp
    if pres.sum()>=10: cr,cp=spearmanr(np.log1p(F.loc[pres,p]),np.log1p(F.loc[pres,q]))
    else: cr=cp=np.nan
    row[f'{reg}_copresent_rho']=cr; row[f'{reg}_copresent_p']=cp
   else:
    for x in ['n_copresent','resid_rho','resid_p','copresent_rho','copresent_p']: row[f'{reg}_{x}']=np.nan
  rows.append(row)
 X=pd.DataFrame(rows).merge(U[['edge_key','weight']].rename(columns={'weight':'upper_FW'}),on='edge_key',how='left').merge(D[['edge_key','weight']].rename(columns={'weight':'deep_FW'}),on='edge_key',how='left')
 for reg in ['upper','deep']:
  mask=X[f'{reg}_resid_p'].notna(); X[f'{reg}_resid_q']=np.nan; X.loc[mask,f'{reg}_resid_q']=multipletests(X.loc[mask,f'{reg}_resid_p'],method='fdr_bh')[1]; X[f'{reg}_supported']=(X[f'{reg}_both_features']==True)&(X[f'{reg}_resid_q']<.05)&(X[f'{reg}_resid_rho'].abs()>=.2)
 X.to_csv(out/'cross_regime_edge_retesting.csv',index=False)
 common=sorted(set(FU.columns)&set(FD.columns)); F=pd.concat([FU[common],FD[common]]); M=pd.concat([MU,MD],axis=0,join='outer').fillna(0).loc[F.index]; M['regime']=(M.depth>220).astype(float); cov=['depth','ctd_temperature','ctd_salinity','ctd_oxygen','ctd_fluorescence']+[c for c in M.columns if c.startswith('season_') or c.startswith('station_')]; cov=[c for c in cov if c in M and M[c].nunique()>1]; C=M[cov].astype(float).copy()
 for c in C: C[c]=(C[c]-C[c].mean())/C[c].std(ddof=0)
 keep=[]; rank=0
 for c in C.columns:
  nr=np.linalg.matrix_rank(C[keep+[c]].values)
  if nr>rank: keep.append(c); rank=nr
 C=C[keep]; test=X[X.upper_both_features & X.deep_both_features]; rr=[]
 for _,e in test.iterrows():
  p,q=e.parasite_id,e.partner_id
  if p not in F or q not in F: continue
  y=np.log1p(F[p].astype(float)).values; partner=np.log1p(F[q].astype(float)).values; partner=(partner-partner.mean())/partner.std(ddof=0); regime=M.regime.values; Z=pd.DataFrame({'partner':partner,'regime':regime,'interaction':partner*regime},index=F.index); Z=pd.concat([Z,C],axis=1); Z=sm.add_constant(Z,has_constant='add'); fit=sm.OLS(y,Z).fit(cov_type='HC3'); bu=fit.params.partner; bi=fit.params.interaction; rr.append([e.edge_key,e.parasite,e.partner,e.partner_group,e.upper_FW,e.deep_FW,len(y),bu,bu+bi,bi,fit.pvalues.interaction,fit.rsquared_adj])
 R=pd.DataFrame(rr,columns=['edge_key','parasite','partner','partner_group','upper_FW','deep_FW','n','upper_partner_slope','deep_partner_slope','slope_difference_deep_minus_upper','interaction_p','adjusted_R2']); R['interaction_q']=multipletests(R.interaction_p,method='fdr_bh')[1]; R['direction']=np.where(R.slope_difference_deep_minus_upper>0,'stronger below 220 m','stronger above 220 m'); R['outcome']='not significant'; R.loc[(R.interaction_q<.05)&(R.slope_difference_deep_minus_upper>0),'outcome']='stronger below'; R.loc[(R.interaction_q<.05)&(R.slope_difference_deep_minus_upper<0),'outcome']='stronger above'; R=R.sort_values('interaction_q'); R.to_csv(out/'parasite_partner_regime_interaction_tests_with_outcome.csv',index=False)
 tests=[]
 for g in sorted(R.partner_group.unique()):
  ing=R.partner_group.eq(g)
  for direction in ['stronger below','stronger above']:
   hit=R.outcome.eq(direction); tab=[[int((ing&hit).sum()),int((ing&~hit).sum())],[int((~ing&hit).sum()),int((~ing&~hit).sum())]]; odds,p=fisher_exact(tab); tests.append([g,direction,tab[0][0],int(ing.sum()),int(hit.sum()),odds,p])
 T=pd.DataFrame(tests,columns=['partner_group','direction','n_group_direction','n_group_total','n_direction_total','odds_ratio','p']); T['q']=multipletests(T.p,method='fdr_bh')[1]; T.to_csv(out/'partner_guild_direction_enrichment.csv',index=False)
 print(f'Upper parasite-protist edges: {len(U)}; deep: {len(D)}'); print(f'Formal interaction tests: {len(R)}; FDR<0.05: {(R.interaction_q<.05).sum()}')
if __name__=='__main__': main()
