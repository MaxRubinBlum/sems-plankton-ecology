#!/usr/bin/env python3
import argparse,pandas as pd

def norm(x):
 x=str(x).strip()
 if x.lower() in {'','nan','unclassified','uncultured','unknown','none'}:return ''
 if '__' in x:x=x.split('__',1)[1]
 return x.replace('_',' ').strip()

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--otu',required=True);ap.add_argument('--taxonomy',required=True);ap.add_argument('--mapping',required=True);ap.add_argument('--candidates',required=True);ap.add_argument('--out-table',required=True);ap.add_argument('--out-ids',required=True);a=ap.parse_args()
 otu=pd.read_csv(a.otu).set_index('OTUID');tax=pd.read_csv(a.taxonomy).set_index('OTUID').fillna('');mp=pd.read_csv(a.mapping,sep='\t').set_index('feature_id');C=pd.read_csv(a.candidates)
 def label(r):
  for rank in ['Species','Genus','Family','Order','Class','Phylum']:
   v=norm(r.get(rank,''))
   if v:return v
  return 'Unclassified Eukaryota'
 def group(r):
  lin=';'.join(norm(r.get(c,'')) for c in ['Phylum','Class','Order','Family','Genus','Species']).lower()
  if 'metazoa' in lin:return 'Metazoa'
  if 'syndin' in lin or 'dino-group-i' in lin or 'dino-group-ii' in lin:return 'Syndiniales'
  if 'hematodinium' in lin:return 'Hematodinium'
  return 'Other eukaryote'
 ann=pd.DataFrame(index=tax.index);ann['taxon_label']=tax.apply(label,axis=1);ann['ecological_group']=tax.apply(group,axis=1);ann['agg_key']=ann.ecological_group+' | '+ann.taxon_label;fidkey={fid:f'{r.ecological_group} | {r.taxon_label}' for fid,r in mp.iterrows()};tot=otu.sum(axis=1);rows=[];ids=[]
 priority=C[C.tier.str.startswith(('Tier 1','Tier 2'))].sort_values('evidence_score',ascending=False).head(20)
 for _,r in priority.iterrows():
  p,q=r.edge_key.split('|');rec=[r.edge_key,r.parasite,r.partner,r.tier,r.literature_class]
  for fid in [p,q]:
   members=ann.index[ann.agg_key==fidkey[fid]].tolist();s=tot.loc[members].sort_values(ascending=False) if members else pd.Series(dtype=float);top=s.head(5);share=(top.iloc[0]/s.sum()) if len(top) and s.sum()>0 else float('nan');rec.extend([fid,len(members),top.index[0] if len(top) else '',share,';'.join(top.index)]);ids.extend(top.index[:3])
  rows.append(rec)
 cols=['edge_key','parasite','partner','tier','literature_class','parasite_feature_id','parasite_constituent_ASVs','parasite_top_ASV','parasite_top_ASV_read_share','parasite_top5_ASVs','partner_feature_id','partner_constituent_ASVs','partner_top_ASV','partner_top_ASV_read_share','partner_top5_ASVs'];pd.DataFrame(rows,columns=cols).to_csv(a.out_table,index=False)
 with open(a.out_ids,'w') as f:
  for x in sorted(set(ids)):f.write(x+'\n')
if __name__=='__main__':main()
