#!/usr/bin/env python3
import argparse,re
from pathlib import Path
import numpy as np,pandas as pd,matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def clean(x):
 x=str(x).replace('Dino-Group-II-Clade-','GII cl. ').replace('Dino-Group-I-Clade-','GI cl. ').replace('Dino-Group-II','GII').replace('Dino-Group-I','GI').replace('-and-','/'); x=re.sub(r'\s+X{1,4}\s+sp\.$','',x); x=re.sub(r'\s+X{1,4}\s+sp\.','',x); return ' '.join(x.replace('-',' ').split())

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--results',required=True); ap.add_argument('--out-prefix',required=True); a=ap.parse_args(); R=pd.read_csv(Path(a.results)/'parasite_partner_regime_interaction_tests_with_outcome.csv')
 assert len(R)==207 and (R.interaction_q<.05).sum()==89 and ((R.interaction_q<.05)&(R.slope_difference_deep_minus_upper>0)).sum()==36 and ((R.interaction_q<.05)&(R.slope_difference_deep_minus_upper<0)).sum()==53
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'axes.linewidth':.8,'pdf.fonttype':42,'ps.fonttype':42}); fig=plt.figure(figsize=(11.4,7.1)); gs=GridSpec(2,2,figure=fig,width_ratios=[1.32,1],height_ratios=[1.08,1],left=.075,right=.985,bottom=.10,top=.965,wspace=.36,hspace=.44)
 ax=fig.add_subplot(gs[:,0])
 for o in ['not significant','stronger above','stronger below']:
  s=R[R.outcome==o]; ax.scatter(s.upper_partner_slope,s.deep_partner_slope,s=16 if o=='not significant' else 27,alpha=.30 if o=='not significant' else .82,label=o)
 allv=np.r_[R.upper_partner_slope,R.deep_partner_slope]; lo,hi=np.nanpercentile(allv,[.5,99.5]); pad=.06*(hi-lo); lo-=pad; hi+=pad; ax.plot([lo,hi],[lo,hi],ls='--',lw=.8); ax.axhline(0,lw=.55,alpha=.45); ax.axvline(0,lw=.55,alpha=.45); ax.set(xlim=(lo,hi),ylim=(lo,hi),xlabel='Partner slope ≤220 m',ylabel='Partner slope >220 m'); ax.legend(frameon=False,fontsize=7.5,loc='upper left'); ax.text(-.10,1.02,'A',transform=ax.transAxes,fontweight='bold',fontsize=11)
 ct=pd.crosstab(R.partner_group,R.outcome)
 for c in ['stronger above','not significant','stronger below']:
  if c not in ct:ct[c]=0
 ct['total']=ct[['stronger above','not significant','stronger below']].sum(axis=1); ct=ct[ct.total>=5].sort_values('total'); ax=fig.add_subplot(gs[0,1]); y=np.arange(len(ct)); left=np.zeros(len(ct))
 for c in ['stronger above','not significant','stronger below']:
  v=ct[c].values; ax.barh(y,v,left=left,height=.64); left+=v
 ax.set_yticks(y,[x.replace(' / ',' /\n') for x in ct.index]); ax.set_xlabel('Tested associations'); ax.text(-.13,1.04,'B',transform=ax.transAxes,fontweight='bold',fontsize=11)
 sig=R[R.interaction_q<.05].copy(); sig['abs_delta']=sig.slope_difference_deep_minus_upper.abs(); top=sig.nlargest(10,'abs_delta').sort_values('slope_difference_deep_minus_upper'); labels=[f'{clean(x)} – {clean(y)}' for x,y in zip(top.parasite,top.partner)]; ax=fig.add_subplot(gs[1,1]); y=np.arange(len(top)); vals=top.slope_difference_deep_minus_upper.values; ax.axvline(0,lw=.7,alpha=.55); ax.hlines(y,0,vals,lw=1.2); ax.scatter(vals,y,s=30,zorder=3); ax.set_yticks(y,labels); ax.set_xlabel('Slope change (>220 m − ≤220 m)'); ax.text(-.13,1.04,'C',transform=ax.transAxes,fontweight='bold',fontsize=11)
 for ax in fig.axes:
  if hasattr(ax,'spines'): ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
 pref=Path(a.out_prefix); fig.savefig(str(pref)+'.png',dpi=600,bbox_inches='tight'); fig.savefig(str(pref)+'.pdf',bbox_inches='tight'); fig.savefig(str(pref)+'.svg',bbox_inches='tight'); plt.close(fig)
if __name__=='__main__': main()
