#!/usr/bin/env python3
import argparse, numpy as np, pandas as pd

def parasite_code(name):
    s=str(name)
    if 'Dino-Group-I-Clade-2' in s:return 'GI_cl2'
    if 'Dino-Group-II-Clade-6' in s:return 'GII_cl6'
    if 'Dino-Group-II-Clade-7' in s:return 'GII_cl7'
    if 'Euduboscquella' in s:return 'Euduboscquella'
    if 'Hematodinium' in s:return 'Hematodinium'
    if 'Syndinium' in s and 'Syndiniales' not in s:return 'Syndinium'
    return 'other'

def literature_class(r):
    p=parasite_code(r.parasite); partner=str(r.partner).lower(); group=str(r.partner_group)
    if p in {'Hematodinium','Syndinium'}: return 'ecological','known crustacean-parasite lineage; protist link treated as ecological association',-2
    if group=='Radiolaria / Acantharea' and p in {'GI_cl2','GII_cl6','GII_cl7'}: return 'cross-study','same Syndiniales clade–radiolarian host class reported independently at BATS',3
    if p=='Euduboscquella' and group=='Ciliates': return 'documented-host-group','Euduboscquella has documented ciliate hosts',3
    if group=='Dinoflagellates':
        if any(k in partner for k in ['prorocentrum','gymnodinium','gyrodinium','akashiwo','alexandrium']): return 'documented-host-genus','partner genus belongs to directly documented Syndiniales host genera',2
        return 'documented-host-group','dinoflagellates are a well-established Syndiniales host group',1
    if group=='Radiolaria / Acantharea': return 'emerging-host-group','independent sequence/network support; direct infection unconfirmed',1
    if group=='Ciliates': return 'documented-host-group','ciliates include documented Syndiniales hosts',1
    return 'none','no direct host precedent used',0

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--cross-regime',required=True); ap.add_argument('--interaction',required=True); ap.add_argument('--clr-sensitivity',required=True); ap.add_argument('--output',required=True); a=ap.parse_args()
    X=pd.read_csv(a.cross_regime); I=pd.read_csv(a.interaction); C=pd.read_csv(a.clr_sensitivity)
    M=X.merge(I[['edge_key','interaction_q','outcome','slope_difference_deep_minus_upper']],on='edge_key',how='left').merge(C[['edge_key','CLR_nonmetazoan_q','CLR_nonmetazoan_slope_change','same_direction']],on='edge_key',how='left')
    rows=[]
    for _,r in M.iterrows():
        posU=bool(r.upper_supported) and pd.notna(r.upper_resid_rho) and r.upper_resid_rho>=.2 and (r.upper_n_copresent if pd.notna(r.upper_n_copresent) else 0)>=20
        posD=bool(r.deep_supported) and pd.notna(r.deep_resid_rho) and r.deep_resid_rho>=.2 and (r.deep_n_copresent if pd.notna(r.deep_n_copresent) else 0)>=20
        if not(posU or posD): continue
        lc,lb,ls=literature_class(r); score=2*int(posU)+2*int(posD)+int(pd.notna(r.upper_FW) or pd.notna(r.deep_FW))+int(pd.notna(r.interaction_q) and r.interaction_q<.05)+int(pd.notna(r.CLR_nonmetazoan_q) and r.CLR_nonmetazoan_q<.05 and bool(r.same_direction))+ls
        if lc in {'cross-study','documented-host-genus','documented-host-group'} and score>=8:tier='Tier 1: strong candidate host association'
        elif lc in {'cross-study','documented-host-genus','documented-host-group','emerging-host-group'} and score>=6:tier='Tier 2: host-plausible association'
        elif score>=5:tier='Tier 3: robust ecological association'
        else:tier='Exploratory'
        if lc=='ecological':tier='Tier 3: robust ecological association'
        rows.append([r.edge_key,r.parasite,r.partner,r.partner_group,posU,posD,r.upper_resid_rho,r.deep_resid_rho,r.upper_n_copresent,r.deep_n_copresent,r.upper_FW,r.deep_FW,r.interaction_q,r.outcome,r.CLR_nonmetazoan_q,r.same_direction,lc,lb,score,tier])
    cols=['edge_key','parasite','partner','partner_group','positive_upper','positive_deep','upper_resid_rho','deep_resid_rho','upper_n_copresent','deep_n_copresent','upper_FW','deep_FW','regime_interaction_q','regime_outcome','CLR_interaction_q','CLR_same_direction','literature_class','literature_basis','evidence_score','tier']
    pd.DataFrame(rows,columns=cols).to_csv(a.output,index=False)
if __name__=='__main__':main()
