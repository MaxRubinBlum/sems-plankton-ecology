#!/usr/bin/env python3
"""Reproduce Figure 4: paired 16S–18S community concordance."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.spatial.distance import pdist, squareform
from scipy.linalg import orthogonal_procrustes
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression

ROOT = Path(__file__).resolve().parents[2]
O16 = ROOT / 'data/raw/16S/otu_table.csv'
O18 = ROOT / 'data/raw/18S/otu_table.csv'
META = ROOT / 'data/processed/18S_samples_CTD_chemistry.csv'
OUT = ROOT / 'results/figure4_cross_domain_concordance'
OUT.mkdir(parents=True, exist_ok=True)

AXES = 8
N_PERM = 9999
SEED = 42
ENV = ['depth','ctd_temperature','ctd_salinity','ctd_oxygen','ctd_fluorescence']
REGIMES = ['A_surface','B_nearsurface','C_DCM','D_below_DCM','F_300-600','G_below_600','H_near_bottom']
REGIME_NAMES = ['Surface','Near-surface','DCM','Below DCM','300–600 m','>600 m','Near-bottom']

missing = [p for p in (O16,O18,META) if not p.exists()]
if missing:
    sys.exit('Missing Figure 4 input(s):\n' + '\n'.join(f'  - {p.relative_to(ROOT)}' for p in missing) +
             '\nSee data/raw/README.md and reproducibility/fig04/README.md for required inputs.')

o16 = pd.read_csv(O16).set_index('OTUID')
o18 = pd.read_csv(O18).set_index('OTUID')
meta = pd.read_csv(META).set_index('sample-id')
ids = sorted(set(o16.columns) & set(o18.columns) & set(meta.index))
if len(ids) != 308:
    raise ValueError(f'Expected 308 paired samples for the validated analysis, found {len(ids)}.')

def pcoa(otu, samples, k=AXES):
    x = otu[samples].T.astype(float)
    x = x.div(x.sum(axis=1), axis=0).fillna(0).values
    d = squareform(pdist(x, metric='braycurtis'))
    n = len(samples); h = np.eye(n) - np.ones((n,n))/n
    g = -0.5 * h @ (d**2) @ h
    w,v = np.linalg.eigh(g); ix = np.argsort(w)[::-1]; w=w[ix]; v=v[:,ix]
    keep = w > 1e-10; w=w[keep]; v=v[:,keep]; k=min(k,len(w))
    return v[:,:k] * np.sqrt(w[:k])

def align(a,b):
    k=min(a.shape[1],b.shape[1]); a=a[:,:k]-a[:,:k].mean(0); b=b[:,:k]-b[:,:k].mean(0)
    a/=np.linalg.norm(a); b/=np.linalg.norm(b)
    rot,_=orthogonal_procrustes(b,a); br=b@rot
    return a,br,float(np.sum(a*br))

def protest(a,b,nperm=N_PERM,seed=SEED):
    _,_,obs=align(a,b); rng=np.random.default_rng(seed)
    null=np.empty(nperm)
    for i in range(nperm):
        _,_,null[i]=align(a,b[rng.permutation(len(b))])
    p=(1+np.sum(null>=obs))/(nperm+1)
    return obs,p

def design(m):
    e=m[ENV].copy(); e['log_depth']=np.log1p(e['depth']); e=e.drop(columns='depth')
    e=pd.concat([e,pd.get_dummies(m['season'],drop_first=True,dtype=float),
                 pd.get_dummies(m['station'],drop_first=True,dtype=float)],axis=1)
    return StandardScaler().fit_transform(e.astype(float))

def residualize(y,x):
    return y-LinearRegression().fit(x,y).predict(x)

A=pcoa(o16,ids); B=pcoa(o18,ids); Ag,Bg,r_global=align(A,B); _,p_global=protest(A,B)

good=meta.loc[ids,ENV].notna().all(axis=1)
ids_adj=list(np.array(ids)[good.values]); m=meta.loc[ids_adj]
if len(ids_adj) != 248:
    raise ValueError(f'Expected 248 complete environmental cases, found {len(ids_adj)}.')
A=pcoa(o16,ids_adj); B=pcoa(o18,ids_adj); X=design(m)
Ar=residualize(A,X); Br=residualize(B,X); Aa,Ba,r_adj=align(Ar,Br); _,p_adj=protest(Ar,Br)

rows=[]
for reg in REGIMES:
    ss=[s for s in ids if meta.loc[s,'ds2']==reg]
    a=pcoa(o16,ss,min(AXES,len(ss)-2)); b=pcoa(o18,ss,min(AXES,len(ss)-2))
    rr,pp=protest(a,b)
    sg=[s for s in ss if meta.loc[s,ENV].notna().all()]
    ra=pa=np.nan
    if len(sg)>=20:
        aa=pcoa(o16,sg,min(AXES,len(sg)-2)); bb=pcoa(o18,sg,min(AXES,len(sg)-2)); xx=design(meta.loc[sg])
        ra,pa=protest(residualize(aa,xx),residualize(bb,xx))
    rows.append([reg,len(ss),rr,pp,len(sg),ra,pa])
summary=pd.DataFrame(rows,columns=['regime','n','procrustes_r','p','n_complete_env','adjusted_r','adjusted_p'])
summary.to_csv(OUT/'community_concordance_summary.csv',index=False)
pd.DataFrame([['whole',len(ids),r_global,p_global],['adjusted',len(ids_adj),r_adj,p_adj]],columns=['analysis','n','procrustes_r','p']).to_csv(OUT/'global_concordance_summary.csv',index=False)

cycle=plt.rcParams['axes.prop_cycle'].by_key()['color']; colors={r:cycle[i] for i,r in enumerate(REGIMES)}
fig=plt.figure(figsize=(12.6,8.1)); gs=fig.add_gridspec(2,2,height_ratios=[1.03,.97],hspace=.42,wspace=.25)
ax1=fig.add_subplot(gs[0,0]); ax2=fig.add_subplot(gs[0,1]); ax3=fig.add_subplot(gs[1,:]); rng=np.random.default_rng(SEED)
for i in rng.choice(len(ids),min(85,len(ids)),replace=False):
    ax1.plot([Ag[i,0],Bg[i,0]],[Ag[i,1],Bg[i,1]],color='.75',lw=.35,alpha=.38,zorder=1)
for reg in REGIMES:
    ii=np.array([meta.loc[s,'ds2']==reg for s in ids])
    ax1.scatter(Ag[ii,0],Ag[ii,1],s=19,facecolors='none',edgecolors=colors[reg],linewidths=.65,zorder=3)
    ax1.scatter(Bg[ii,0],Bg[ii,1],s=12,c=[colors[reg]],alpha=.55,linewidths=0,zorder=2)
ax1.text(.02,.97,'a',transform=ax1.transAxes,va='top',fontsize=13,fontweight='bold'); ax1.text(.98,.97,f'$r$ = {r_global:.2f}',transform=ax1.transAxes,ha='right',va='top',fontsize=10)
ax1.set_title('Whole water column',fontsize=10); ax1.set_xlabel('Procrustes dimension 1'); ax1.set_ylabel('Procrustes dimension 2')
for i in rng.choice(len(ids_adj),min(70,len(ids_adj)),replace=False):
    ax2.plot([Aa[i,0],Ba[i,0]],[Aa[i,1],Ba[i,1]],color='.78',lw=.35,alpha=.35,zorder=1)
ax2.scatter(Aa[:,0],Aa[:,1],s=18,facecolors='none',edgecolors='.20',linewidths=.6,zorder=3); ax2.scatter(Ba[:,0],Ba[:,1],s=12,c='.48',alpha=.58,linewidths=0,zorder=2)
ax2.text(.02,.97,'b',transform=ax2.transAxes,va='top',fontsize=13,fontweight='bold'); ax2.text(.98,.97,f'$r$ = {r_adj:.2f}',transform=ax2.transAxes,ha='right',va='top',fontsize=10)
ax2.set_title('Environment-adjusted',fontsize=10); ax2.set_xlabel('Procrustes dimension 1'); ax2.set_ylabel('Procrustes dimension 2')
handles=[Line2D([0],[0],marker='o',linestyle='None',markerfacecolor='none',markeredgecolor='.2',markersize=5,label='16S'),Line2D([0],[0],marker='o',linestyle='None',markerfacecolor='.48',markeredgecolor='.48',markersize=4,label='18S')]
fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.52,.535),frameon=False,ncol=2,fontsize=8)
y=np.arange(len(REGIMES)); raw=summary.procrustes_r.values; adj=summary.adjusted_r.values
for i in range(len(y)):
    if np.isfinite(adj[i]): ax3.plot([raw[i],adj[i]],[y[i],y[i]],color='.72',lw=1.1,zorder=1)
ax3.scatter(raw,y,s=50,facecolors='white',edgecolors='.25',linewidths=1,label='Unadjusted',zorder=3)
mask=np.isfinite(adj); ax3.scatter(adj[mask],y[mask],s=50,c='.35',linewidths=0,label='Adjusted',zorder=4)
labs=[]
for i,nm in enumerate(REGIME_NAMES):
    labs.append(f'{nm}  (n={summary.n.iloc[i]}; adjusted n={summary.n_complete_env.iloc[i]})' if np.isfinite(adj[i]) else f'{nm}  (n={summary.n.iloc[i]})')
ax3.set_yticks(y,labs); ax3.invert_yaxis(); ax3.set_xlim(.45,.95); ax3.set_xlabel('Procrustes correlation ($r$)'); ax3.text(.005,.98,'c',transform=ax3.transAxes,fontsize=13,fontweight='bold',va='top'); ax3.legend(frameon=False,fontsize=8,loc='lower right',ncol=2); ax3.grid(axis='x',lw=.5,alpha=.18)
for ax in (ax1,ax2,ax3): ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
fig.subplots_adjust(left=.13,right=.98,top=.95,bottom=.08)
fig.savefig(OUT/'Figure4_cross_domain_concordance.svg',bbox_inches='tight'); fig.savefig(OUT/'Figure4_cross_domain_concordance.pdf',bbox_inches='tight'); plt.close(fig)
