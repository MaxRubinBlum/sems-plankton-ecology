#!/usr/bin/env python3
"""Integrate HaiSec CTD casts with SEMS 16S/18S molecular metadata.

Outputs interpolated temperature, salinity, oxygen and fluorescence at each
molecular sampling depth, plus a per-cast fluorescence-maximum (DCM proxy)
and sample depth relative to that maximum.
"""
from pathlib import Path
import argparse, re
import numpy as np
import pandas as pd

MONTHS={m.lower():i for i,m in enumerate(['january','february','march','april','may','june','july','august','september','october','november','december'],1)}

def norm_station(s):
    s=str(s).strip().upper()
    s=re.sub(r'^[A-Z]+_(?=H\d)', '', s)
    s=s.replace('HS_','H').replace('HS','H')
    return s

def sample_cruise_num(s):
    m=re.match(r'HS(\d{2})_',str(s),re.I)
    return int(m.group(1)) if m else np.nan

def load_ctd(path):
    ctd=pd.read_csv(path,sep='\t',encoding='latin1')
    ctd.columns=[c.strip() for c in ctd.columns]
    ctd['datetime']=pd.to_datetime(ctd['yyyy-mm-ddThh:mm:ss.sss'],errors='coerce')
    ctd['year']=ctd.datetime.dt.year; ctd['month_num']=ctd.datetime.dt.month
    ctd['cruise_num']=ctd.Cruise.str.extract(r'(\d+)$').astype(int)
    ctd['station_norm']=ctd.Station.map(norm_station)
    for c in ['Pressure [db]','Salinity']:
        ctd[c]=pd.to_numeric(ctd[c],errors='coerce')
    for c in ctd.columns:
        if c.startswith('Temperature') or 'Oxygen' in c or 'Fluoresc' in c or c=='Bot. Depth [m]':
            ctd[c]=pd.to_numeric(ctd[c],errors='coerce')
    return ctd

def cast_inventory(ctd):
    return (ctd.groupby(['Cruise','Station','datetime'])
            .agg(year=('year','first'),month=('month_num','first'),
                 lon=('Longitude [degrees_east]','median'),lat=('Latitude [degrees_north]','median'),
                 bottom_depth_m=('Bot. Depth [m]','median'),max_pressure_db=('Pressure [db]','max'),
                 n_points=('Pressure [db]','size')).reset_index())

def descriptors(ctd):
    fcol=next(c for c in ctd.columns if 'Fluoresc' in c); p='Pressure [db]'
    rows=[]
    for (cn,st),z in ctd.groupby(['cruise_num','station_norm']):
        z=z.sort_values(p).dropna(subset=[p]).drop_duplicates(p); maxp=z[p].max()
        dcm=dcmf=np.nan
        if maxp>=30:
            q=z[(z[p]>=2)&(z[p]<=min(200,maxp))].dropna(subset=[fcol]).copy()
            if len(q):
                q['f_smooth']=q[fcol].rolling(5,center=True,min_periods=1).median()
                ix=q.f_smooth.idxmax(); dcm=float(q.loc[ix,p]); dcmf=float(q.loc[ix,'f_smooth'])
        rows.append(dict(cruise_num=cn,station_norm=st,bottom_depth_m=z['Bot. Depth [m]'].median(),
                         max_pressure_db=maxp,dcm_pressure_db=dcm,dcm_fluorescence=dcmf))
    return pd.DataFrame(rows)

def match_and_interpolate(meta,ctd,desc):
    meta=meta.copy(); meta['station_norm']=meta.station.map(norm_station)
    meta['month_num']=meta.month.astype(str).str.lower().map(MONTHS)
    tcol=next(c for c in ctd.columns if c.startswith('Temperature')); ocol=next(c for c in ctd.columns if 'Oxygen' in c); fcol=next(c for c in ctd.columns if 'Fluoresc' in c); p='Pressure [db]'
    rows=[]
    for _,r in meta.iterrows():
        cn=sample_cruise_num(r['sample-id'])
        if pd.notna(cn):
            z=ctd[(ctd.cruise_num==int(cn))&(ctd.station_norm==r.station_norm)]
            method='sample_cruise_id'
        else:
            z=ctd[(ctd.year==r.year)&(ctd.month_num==r.month_num)&(ctd.station_norm==r.station_norm)]
            method='year_month_station'
        rec=r.to_dict(); rec.update(match_method=method,n_candidate_casts=0,cruise='',ctd_interp_ok=False,
            ctd_temperature=np.nan,ctd_salinity=np.nan,ctd_oxygen=np.nan,ctd_fluorescence=np.nan,dcm_pressure_db=np.nan,delta_depth_from_dcm_m=np.nan)
        cruises=z.Cruise.unique(); rec['n_candidate_casts']=len(cruises)
        if len(cruises)==1:
            rec['cruise']=cruises[0]; cn=int(str(cruises[0]).split('_')[-1]); q=z.sort_values(p).dropna(subset=[p]).drop_duplicates(p); dep=float(r.depth)
            if len(q)>=2 and dep>=q[p].min()-2 and dep<=q[p].max()+5:
                def interp(col):
                    w=q[[p,col]].dropna(); return float(np.interp(dep,w[p],w[col])) if len(w)>=2 else np.nan
                ce=desc[(desc.cruise_num==cn)&(desc.station_norm==r.station_norm)]
                dcm=float(ce.dcm_pressure_db.iloc[0]) if len(ce) else np.nan
                rec.update(ctd_temperature=interp(tcol),ctd_salinity=interp('Salinity'),ctd_oxygen=interp(ocol),ctd_fluorescence=interp(fcol),dcm_pressure_db=dcm,delta_depth_from_dcm_m=dep-dcm if np.isfinite(dcm) else np.nan,ctd_interp_ok=True)
        rows.append(rec)
    return pd.DataFrame(rows)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--ctd',required=True); ap.add_argument('--meta16',required=True); ap.add_argument('--meta18',required=True); ap.add_argument('--outdir',default='results/ctd_integration'); a=ap.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    ctd=load_ctd(a.ctd); desc=descriptors(ctd); desc.to_csv(out/'CTD_cast_environmental_descriptors.csv',index=False)
    cast_inventory(ctd).to_csv(out/'CTD_cast_inventory.csv',index=False)
    for label,path in [('16S',a.meta16),('18S',a.meta18)]:
        e=match_and_interpolate(pd.read_csv(path),ctd,desc); e.to_csv(out/f'{label}_samples_with_CTD.csv',index=False)
        print(label, 'samples=',len(e),'interpolated=',int(e.ctd_interp_ok.sum()),'with_DCM=',int(e.delta_depth_from_dcm_m.notna().sum()))
if __name__=='__main__': main()
