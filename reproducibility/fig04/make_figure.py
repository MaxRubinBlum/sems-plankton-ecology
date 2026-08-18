#!/usr/bin/env python3
"""Reproduce manuscript Figure 4 from validated result tables.

The manuscript figure deliberately uses original SILVA taxonomy. No GTDB or
manual cross-taxonomy substitutions are made. This script expects the three
CSV inputs documented in inputs.tsv.
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

ROOT = Path(__file__).resolve().parents[2]
WHOLE = ROOT / 'results/eastmed_crossdomain_concordance/16S_18S_whole_community_concordance.csv'
DEPTH = ROOT / 'results/eastmed_crossdomain_concordance/16S_18S_concordance_by_depth_and_regime.csv'
STABLE = ROOT / 'results/eastmed_trait_crossdomain/protist_traits_vs_16S_stable_links.csv'
OUT = ROOT / 'results/fig04'
OUT.mkdir(parents=True, exist_ok=True)

whole = pd.read_csv(WHOLE)
depth = pd.read_csv(DEPTH)
stable = pd.read_csv(STABLE)

cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
depth_order = ['A_surface','B_nearsurface','C_DCM','D_below_DCM','F_300-600','G_below_600','H_near_bottom']
depth_names = {'A_surface':'Surface','B_nearsurface':'Near-surface','C_DCM':'DCM','D_below_DCM':'Below DCM','F_300-600':'300–600 m','G_below_600':'>600 m','H_near_bottom':'Near-bottom'}

dd = depth[depth.subset.isin(depth_order)].copy()
dd['ord'] = dd.subset.map({k:i for i,k in enumerate(depth_order)})
dd = dd.sort_values('ord')
rows = [('Whole community', whole.iloc[0].rho, whole.iloc[0].n_shared, '.25'),
        ('After hydrography', whole.iloc[1].rho, whole.iloc[1].n_shared, '.45')]
for _, r in dd.iterrows():
    rows.append((depth_names[r.subset], r.rho, r.n, cycle[int(r.ord)]))
forest = pd.DataFrame(rows, columns=['label','rho','n','color'])

s = stable[(stable.FDR_meta < 0.05) & (stable.stability_fraction >= 5/7)].copy()
s['absrho'] = s.pooled_within_depth_rho.abs()
rank = (s.groupby('prokaryote').agg(n_traits=('trait','nunique'), max_abs_rho=('absrho','max'), mean_stability=('stability_fraction','mean')).sort_values(['n_traits','max_abs_rho','mean_stability'], ascending=[False,False,False]))
top = rank.head(12).index
traits = ['bona_fide_phototrophy','constitutive_mixotrophy','parasitism','phagotrophy_radiolaria','diplonemid_heterotrophy']
hm = s[s.prokaryote.isin(top)].pivot_table(index='prokaryote', columns='trait', values='pooled_within_depth_rho', aggfunc='first').reindex(columns=traits)
a = hm.fillna(0).to_numpy()
hm = hm.iloc[np.lexsort((-np.max(np.abs(a), axis=1), np.argmax(np.abs(a), axis=1)))]

def display_label(x):
    x = re.sub(r'^(Genus|Family|Order|Class):', '', str(x))
    return x.replace('g__','').replace('f__','').replace('o__','').replace('c__','').replace('_',' ').strip()
hm.index = [display_label(x) for x in hm.index]

fig = plt.figure(figsize=(11.4,6.5))
gs = fig.add_gridspec(1,2,width_ratios=[.92,1.55],wspace=.38)
ax = fig.add_subplot(gs[0,0])
for i,r in forest.iterrows():
    ax.plot([0,r.rho],[i,i],lw=2.2,color='.75')
    ax.scatter(r.rho,i,s=72,color=r.color,edgecolor='white',linewidth=.7,zorder=3)
    ax.text(.985,i,f'n={int(r.n)}',transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=7,color='.4')
ax.axhline(1.5,color='.82',lw=.8)
ax.set_yticks(range(len(forest)),forest.label,fontsize=8.2)
ax.set_xlim(0,1.02); ax.set_xlabel('16S–18S concordance (Spearman ρ)'); ax.invert_yaxis()
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False); ax.grid(axis='x',alpha=.12,lw=.6)

ax = fig.add_subplot(gs[0,1])
norm = TwoSlopeNorm(vmin=-.6,vcenter=0,vmax=.6); cmap=plt.get_cmap('coolwarm')
for i in range(len(hm)):
    for j in range(hm.shape[1]):
        v=hm.iloc[i,j]
        if pd.notna(v): ax.scatter(j,i,s=30+360*(abs(v)/.6)**1.6,c=[cmap(norm(v))],edgecolor='white',linewidth=.7)
ax.set_xlim(-.6,4.6); ax.set_ylim(len(hm)-.6,-.6)
ax.set_xticks(range(5),['Phototrophy','Mixotrophy','Parasitism','Radiolarian\nphagotrophy','Diplonemid\nheterotrophy'],fontsize=8)
ax.set_yticks(range(len(hm)),hm.index,fontsize=7.6); ax.tick_params(length=0)
for x in np.arange(.5,4.6,1): ax.axvline(x,color='.92',lw=.6,zorder=0)
for y in np.arange(.5,len(hm)-.4,1): ax.axhline(y,color='.95',lw=.5,zorder=0)
for sp in ax.spines.values(): sp.set_visible(False)
sm=plt.cm.ScalarMappable(norm=norm,cmap=cmap)
cax=fig.add_axes([.63,.065,.19,.018]); cb=fig.colorbar(sm,cax=cax,orientation='horizontal'); cb.set_label('Within-depth Spearman ρ',fontsize=7.3); cb.ax.tick_params(labelsize=6.5)
lax=fig.add_axes([.83,.035,.13,.075]); lax.axis('off'); vals=[.2,.4,.6]
handles=[plt.scatter([],[],s=30+360*(v/.6)**1.6,color='.55',edgecolor='white',linewidth=.7) for v in vals]
lax.legend(handles,[f'{v:.1f}' for v in vals],title='|ρ|',frameon=False,ncol=3,loc='center',fontsize=6.7,title_fontsize=7,handletextpad=.3,columnspacing=.7,borderpad=0)
fig.subplots_adjust(left=.16,right=.97,top=.97,bottom=.16)
fig.savefig(OUT/'Figure4_cross_domain_concordance_SILVA.svg',bbox_inches='tight')
fig.savefig(OUT/'Figure4_cross_domain_concordance_SILVA.pdf',bbox_inches='tight')
hm.to_csv(OUT/'Figure4_SILVA_selected_trait_links.csv')
forest.drop(columns='color').to_csv(OUT/'Figure4_concordance_summary.csv',index=False)
