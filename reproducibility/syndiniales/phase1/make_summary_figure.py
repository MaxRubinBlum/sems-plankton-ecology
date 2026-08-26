#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np,pandas as pd,matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--results',required=True); ap.add_argument('--out-prefix',required=True); a=ap.parse_args(); r=Path(a.results)
    scan=pd.read_csv(r/'syndiniales_objective_boundary_scan.csv'); pcoa=pd.read_csv(r/'syndiniales_PCoA_scores.csv'); var=pd.read_csv(r/'syndiniales_PCoA_variance.csv'); metrics=pd.read_csv(r/'syndiniales_sample_metrics.csv'); cl=pd.read_csv(r/'syndiniales_clade_vertical_statistics.csv'); seas=pd.read_csv(r/'syndiniales_season_boundary_scans.csv'); s=pd.read_csv(r/'syndiniales_boundary_summary.csv'); bp=float(s.observed_best_threshold_m.iloc[0])
    bins=np.array([0,25,50,80,120,160,200,250,350,500,750,1100,2000]); metrics['bin']=pd.cut(metrics.depth,bins,include_lowest=True); m=metrics.groupby('bin',observed=True).agg(depth=('depth','median'),fraction=('syndiniales_fraction_18S','median'),richness=('syndiniales_richness','median'),shannon=('syndiniales_shannon','median')).reset_index(); top=cl.nlargest(8,'mean_fraction_within_Syndiniales').sort_values('spearman_depth_rho'); top['short']=top.clade.str.replace('Dino-Group-','',regex=False).str.replace('Clade-','cl. ',regex=False).str.replace('-and-','/',regex=False).str.replace('-',' ',regex=False)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'axes.linewidth':.8,'pdf.fonttype':42,'ps.fonttype':42}); fig=plt.figure(figsize=(11,7.6)); gs=GridSpec(2,3,figure=fig,width_ratios=[1.05,1.25,1],left=.075,right=.98,bottom=.09,top=.97,wspace=.38,hspace=.42)
    ax=fig.add_subplot(gs[0,0]); ax.plot(scan.threshold_m,scan.partial_R2,lw=2,label='All samples')
    for sea,x in seas.groupby('season'): ax.plot(x.threshold_m,x.partial_R2,lw=1.25,ls='--',label='Winter' if 'Winter' in sea else 'Summer')
    ax.axvline(bp,lw=.9,ls=':'); ax.set(xlabel='Candidate boundary depth (m)',ylabel='Partial $R^2$'); ax.legend(frameon=False,fontsize=7.5); ax.text(-.17,1.06,'A',transform=ax.transAxes,fontweight='bold',fontsize=11)
    ax=fig.add_subplot(gs[0,1]); sc=ax.scatter(pcoa.PCo1,pcoa.PCo2,c=np.log10(pcoa.depth+1),s=16,alpha=.78,linewidths=0); ax.set(xlabel=f"PCoA1 ({100*var.variance_fraction.iloc[0]:.1f}%)",ylabel=f"PCoA2 ({100*var.variance_fraction.iloc[1]:.1f}%)"); cb=fig.colorbar(sc,ax=ax,pad=.02,fraction=.05); cb.set_label('log$_{10}$(depth + 1)'); ax.text(-.14,1.06,'B',transform=ax.transAxes,fontweight='bold',fontsize=11)
    ax=fig.add_subplot(gs[0,2]); ax.plot(m.fraction,m.depth,marker='o',ms=4,lw=1.4); ax.axhline(bp,lw=.9,ls=':'); ax.invert_yaxis(); ax.set(xlabel='Syndiniales fraction of 18S reads',ylabel='Depth (m)'); ax.text(-.18,1.06,'C',transform=ax.transAxes,fontweight='bold',fontsize=11)
    ax=fig.add_subplot(gs[1,0]); ax.plot(m.richness/m.richness.max(),m.depth,marker='o',ms=3.5,lw=1.35,label='Richness'); ax.plot(m.shannon/m.shannon.max(),m.depth,marker='s',ms=3.2,lw=1.25,label='Shannon'); ax.axhline(bp,lw=.9,ls=':'); ax.invert_yaxis(); ax.set(xlabel='Scaled diversity metric',ylabel='Depth (m)',xlim=(0,1.05)); ax.legend(frameon=False,fontsize=7.5); ax.text(-.17,1.06,'D',transform=ax.transAxes,fontweight='bold',fontsize=11)
    ax=fig.add_subplot(gs[1,1:]); y=np.arange(len(top)); ax.axvline(0,lw=.8); ax.scatter(top.spearman_depth_rho,y,s=38,zorder=3)
    for yy,v in zip(y,top.spearman_depth_rho): ax.plot([0,v],[yy,yy],lw=1.2)
    ax.set_yticks(y,top.short); ax.set(xlabel='Spearman correlation with depth',xlim=(-1,1)); ax.text(-.07,1.06,'E',transform=ax.transAxes,fontweight='bold',fontsize=11)
    for ax in fig.axes:
        if hasattr(ax,'spines'): ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    pref=Path(a.out_prefix); fig.savefig(str(pref)+'.png',dpi=600,bbox_inches='tight'); fig.savefig(str(pref)+'.pdf',bbox_inches='tight'); fig.savefig(str(pref)+'.svg',bbox_inches='tight'); plt.close(fig)
if __name__=='__main__': main()
