# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 18:14:24 2019

@author: Conor Jenkins
"""
import pandas as pd
from scipy.stats import pearsonr 
from scipy.stats import spearmanr

proteomic_df = pd.read_excel('SingleShot_6_thru_12_EggNOG_plus_Viridiplantae_SeQuest.xlsx')

metabolo_df = pd.read_csv('Metabolite_data-from-portal.csv')
metabolo_df = metabolo_df[['Best molecule match 7/18/19','Group Area: 1',
                           'Group Area: 2','Group Area: 3','Group Area: 4',
                           'Group Area: 5','Group Area: 6']]
metabolo_df.rename(columns={'Best molecule match 7/18/19':'Compound',
                            'Group Area: 1':'Sample 1','Group Area: 2':'Sample 2',
                            'Group Area: 3':'Sample 3','Group Area: 4':'Sample 4',
                            'Group Area: 5':'Sample 5','Group Area: 6':'Sample 6'}, 
                            inplace=True)

proteomic_df = proteomic_df[['Accession','Abundances (Grouped): F1',
                           'Abundances (Grouped): F2','Abundances (Grouped): F3',
                           'Abundances (Grouped): F4','Abundances (Grouped): F5',
                           'Abundances (Grouped): F6']]

proteomic_df.rename(columns={'Abundances (Grouped): F1':'Sample 1',
                    'Abundances (Grouped): F2': 'Sample 2', 'Abundances (Grouped): F3':'Sample 3',
                    'Abundances (Grouped): F4':'Sample 4', 'Abundances (Grouped): F5':'Sample 5',
                    'Abundances (Grouped): F6':'Sample 6'}, inplace=True)

proteomic_df.fillna(0)

header = ['Metabolite', 'Protein', 'Pearson_Corr', 'Pearson_pvalue', 'Spearman_Corr', 'Spearman_pvalue']

corr_df = pd.DataFrame(columns=header)

for index, row in metabolo_df.iterrows():
    metabolite = row[0]
    print(metabolite)
    x = row[1:]
    for number, entry in proteomic_df.iterrows():
        protein = entry[0]
        y = entry[1:]
        spear_corr = spearmanr(x,y)
        pear_corr = pearsonr(x,y)
        spearman_correlation = spear_corr[0]
        spearman_pvalue = spear_corr[1]
        pearson_correlation = pear_corr[0]
        pearson_pvalue = pear_corr[1]
        if spearman_pvalue <= 0.05 or pearson_pvalue <= 0.05:
            data_df = pd.DataFrame([[metabolite,protein,pearson_correlation,pearson_pvalue,spearman_correlation,spearman_pvalue]],columns=header)
            corr_df = corr_df.append(data_df)
corr_df.to_csv('results.csv',index=False)
