"""Reproduce supplementary environmental PCA from integrated 18S environmental data."""
from pathlib import Path
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

ROOT=Path(__file__).resolve().parents[2]
INPUT=ROOT/'data'/'processed'/'metadata'/'18S_samples_CTD_chemistry.csv'
OUT=ROOT/'results'/'supp_environment'; OUT.mkdir(parents=True,exist_ok=True)
CTD=['ctd_temperature','ctd_salinity','ctd_oxygen','ctd_fluorescence']
CHEM_COMMON=['chem_no3_no2_umol_kg','chem_po4_umol_kg','chem_silicate_umol_kg','chem_pH_total_25C']

def run_pca(df,variables):
    z=df.dropna(subset=variables).copy(); X=StandardScaler().fit_transform(z[variables]); model=PCA().fit(X); scores=model.transform(X)
    z['PC1']=scores[:,0]; z['PC2']=scores[:,1]; loadings=pd.DataFrame(model.components_[:2].T,index=variables,columns=['PC1','PC2']); variance=100*model.explained_variance_ratio_[:2]
    return z,loadings,variance

if __name__=='__main__':
    df=pd.read_csv(INPUT); ctd_scores,ctd_loadings,ctd_var=run_pca(df,CTD); env_scores,env_loadings,env_var=run_pca(df,CTD+CHEM_COMMON)
    ctd_scores.to_csv(OUT/'panel_a_CTD_PCA_scores.csv',index=False); ctd_loadings.to_csv(OUT/'panel_a_CTD_PCA_loadings.csv')
    env_scores.to_csv(OUT/'panel_b_CTD_chemistry_PCA_scores.csv',index=False); env_loadings.to_csv(OUT/'panel_b_CTD_chemistry_PCA_loadings.csv')
    print('Panel a:',len(ctd_scores),ctd_var); print('Panel b:',len(env_scores),env_var)
