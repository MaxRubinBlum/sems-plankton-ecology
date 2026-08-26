#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np,pandas as pd
BOUNDARY=220.0

def norm(x):
 x=str(x).strip()
 if x.lower() in {'','nan','unclassified','uncultured','unknown','none'}: return ''
 if '__' in x:x=x.split('__',1)[1]
 return x.replace('_',' ').strip()

def tax_label(r):
 for rank in ['Species','Genus','Family','Order','Class','Phylum']:
  v=norm(r.get(rank,''))
  if v:return v,rank
 return 'Unclassified Eukaryota','Kingdom'

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--otu',required=True); ap.add_argument('--taxonomy',required=True); ap.add_argument('--metadata',required=True); ap.add_argument('--outdir',required=True); a=ap.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
 otu=pd.read_csv(a.otu).set_index('OTUID'); tax=pd.read_csv(a.taxonomy).set_index('OTUID').fillna(''); meta=pd.read_csv(a.metadata).set_index('sample-id'); envbase=['depth','ctd_temperature','ctd_salinity','ctd_oxygen','ctd_fluorescence']; ids=[s for s in otu.columns if s in meta.index]; ids=[s for s in ids if meta.loc[s,envbase].notna().all()]
 labels=[]; ranks=[]; groups=[]; parasites=[]
 for _,r in tax.iterrows():
  lab,rk=tax_label(r); lineage=';'.join(norm(r.get(c,'')) for c in ['Phylum','Class','Order','Family','Genus','Species']).lower(); syn=('syndin' in lineage or 'dino-group-i' in lineage or 'dino-group-ii' in lineage); hem=('hematodinium' in lineage); para=syn or hem; group='Metazoa' if 'metazoa' in lineage else ('Syndiniales' if syn else ('Hematodinium' if hem else 'Other eukaryote')); labels.append(lab); ranks.append(rk); groups.append(group); parasites.append(para)
 ann=pd.DataFrame({'taxon_label':labels,'rank':ranks,'ecological_group':groups,'parasite_flag':parasites},index=tax.index); ann['agg_key']=ann.ecological_group+' | '+ann.taxon_label; x=otu[ids].astype(float).copy(); x['agg_key']=ann.reindex(x.index).agg_key.values; agg=x.groupby('agg_key').sum(); prev=(agg>0).mean(axis=1); totals=agg.sum(axis=1); grp=ann.groupby('agg_key').ecological_group.first(); para=ann.groupby('agg_key').parasite_flag.max(); keep=(para&(prev>=.05)&(totals>=500)) | (~para&(prev>=.03)&(totals>=500)); table=agg.loc[keep].T
 mapping=[]; F=pd.DataFrame(index=table.index)
 for i,key in enumerate(table.columns,1):
  fid=f'E18_{i:04d}'; F[fid]=table[key].values; mapping.append([fid,key.split(' | ',1)[1],grp[key],bool(para[key]),float(prev[key]),float(totals[key])])
 mapping=pd.DataFrame(mapping,columns=['feature_id','taxon_label','ecological_group','parasite_flag','prevalence','total_reads']).set_index('feature_id'); M=meta.loc[ids,['depth','ctd_temperature','ctd_salinity','ctd_oxygen','ctd_fluorescence','season','station']].copy(); M['log_depth']=np.log1p(M.depth); M=pd.concat([M.drop(columns=['season','station']),pd.get_dummies(M.season,prefix='season',dtype=int),pd.get_dummies(M.station,prefix='station',dtype=int)],axis=1)
 summaries=[]
 for name,sel in [('upper_le_220m',M.index[M.depth<=BOUNDARY]),('deep_gt_220m',M.index[M.depth>BOUNDARY])]:
  f=F.loc[sel].copy(); prevalence=(f>0).mean(); tt=f.sum(); par=mapping.reindex(f.columns).parasite_flag.astype(bool); minreads=max(100,2*len(f)); keep2=(par&(prevalence>=.05)&(tt>=minreads)) | (~par&(prevalence>=.03)&(tt>=minreads)); f=f.loc[:,keep2]; f=f.loc[:,f.nunique()>1]; md=M.loc[sel].copy(); md=md.loc[:,md.nunique()>1]; f.index.name='sample_id'; md.index.name='sample_id'; f.to_csv(out/f'features_{name}.tsv',sep='\t'); md.to_csv(out/f'metadata_{name}.tsv',sep='\t'); mapping.loc[f.columns].reset_index().to_csv(out/f'mapping_{name}.tsv',sep='\t',index=False); summaries.append([name,len(f),f.shape[1],int(mapping.loc[f.columns].parasite_flag.sum()),minreads])
 pd.DataFrame(summaries,columns=['regime','samples','features','parasite_features','min_total_reads']).to_csv(out/'phase2_input_summary.csv',index=False)
if __name__=='__main__':main()
