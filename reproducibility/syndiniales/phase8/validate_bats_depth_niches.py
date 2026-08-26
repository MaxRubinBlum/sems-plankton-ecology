#!/usr/bin/env python3
import argparse,pandas as pd
from scipy.stats import binomtest
EXPECTED={'Dino-Group-II-Clade-10-and-11':'photic','Dino-Group-I-Clade-1':'photic','Dino-Group-I-Clade-5':'photic','Dino-Group-II-Clade-6':'aphotic','Dino-Group-II-Clade-7':'aphotic','Dino-Group-I-Clade-2':'aphotic'}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--clade-stats',required=True);ap.add_argument('--output',required=True);a=ap.parse_args();cs=pd.read_csv(a.clade_stats);rows=[]
 for c,e in EXPECTED.items():
  r=cs[cs.clade==c].iloc[0];o='aphotic' if r.spearman_depth_rho>0 else 'photic';rows.append([c,e,r.spearman_depth_rho,o,o==e])
 R=pd.DataFrame(rows,columns=['clade','BATS_reported_niche','EastMed_depth_rho','EastMed_direction','direction_agreement']);R.to_csv(a.output,index=False);print(binomtest(R.direction_agreement.sum(),len(R),p=.5,alternative='greater'))
if __name__=='__main__':main()
