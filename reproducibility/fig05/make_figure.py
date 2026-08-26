from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
IN=ROOT/'results/figure5_cross_domain_associations/FlashWeave_14_pairs_grouped_by_18S.csv'
OUT=ROOT/'results/figure5_cross_domain_associations'
d=pd.read_csv(IN)
order=['Syndiniales','Parasitic dinoflagellates','Radiolaria / Acantharea','Diplonemids','MAST','Labyrinthulomycetes','Pelagophytes','Amoebozoa','Unresolved eukaryote']
d['group18_precise']=pd.Categorical(d.group18_precise,categories=order,ordered=True)
d=d.sort_values(['group18_precise','fw_weight'],ascending=[True,False]).reset_index(drop=True)
rows=[]; y=0.; ranges={}
for g in order:
    sub=d[d.group18_precise.astype(str)==g]
    if sub.empty: continue
    y+=.9; start=y
    for _,r in sub.iterrows(): rows.append((y,g,r)); y+=1.15
    ranges[g]=(start,y-1.15); y+=.7
fig,ax=plt.subplots(figsize=(11.2,8.6)); ax.set_xlim(0,1); ax.set_ylim(-.4,y+.5); ax.invert_yaxis(); ax.axis('off')
ax.text(.32,.25,'16S partner',ha='right',va='bottom',fontsize=9.5,fontweight='bold'); ax.text(.68,.25,'18S partner',ha='left',va='bottom',fontsize=9.5,fontweight='bold'); ax.text(.97,.25,'18S affiliation',ha='right',va='bottom',fontsize=9.5,fontweight='bold')
for g,(start,end) in ranges.items(): ax.plot([.03,.95],[start-.55,start-.55],lw=.55,color='.88'); ax.text(.97,(start+end)/2,g,ha='right',va='center',fontsize=8.1,fontweight='bold')
maxw=d.fw_weight.abs().max()
for yy,g,r in rows:
    ax.text(.32,yy,r.tax16,ha='right',va='center',fontsize=8.4); ax.text(.68,yy,r.tax18,ha='left',va='center',fontsize=8.4)
    lw=.9+3.1*abs(r.fw_weight)/maxw; ax.plot([.38,.62],[yy,yy],lw=lw,color='.20',solid_capstyle='round'); ax.text(.50,yy-.19,f'{r.fw_weight:.2f}',ha='center',va='bottom',fontsize=6.6,color='.42')
legend_y=y-.15; ax.text(.03,legend_y,'FlashWeave weight',fontsize=7.3,va='center')
for k,(w,label) in enumerate([(0.15,'0.15'),(.30,'0.30'),(.45,'0.45')]):
    x=.18+k*.14; ax.plot([x,x+.055],[legend_y,legend_y],lw=.9+3.1*w/maxw,color='.20',solid_capstyle='round'); ax.text(x+.065,legend_y,label,fontsize=6.7,va='center')
fig.tight_layout()
for ext in ['svg','pdf']: fig.savefig(OUT/f'Figure5_cross_domain_associations.{ext}',bbox_inches='tight')
plt.close(fig)
