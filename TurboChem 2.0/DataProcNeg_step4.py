#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 14:11:56 2022

@author: dabrahamsson
"""

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import matplotlib.pyplot as plt


def read_and_prepare():
    df1 = pd.read_csv('1stExpMergedCleanNeg.csv')
    df1 = df1.drop('Unnamed: 0', axis=1)
    df1['chem_id'] = (np.round(df1['Average Mz'], 3).astype(str) + 
                    '@' + np.round(df1['Average Rt(min)'], 3).astype(str))
    df2 = pd.read_csv('neg_combined_csv.csv')
    df = pd.merge(df1, df2, on='chem_id', how='left')
    sol = ['1-octanol', 'butyl acetate', 'chloroform', 'cyclohexane',
           'dichloromethane', 'triolein', 'n-hexane', 'n-octane', 'oleylalcohol',
           'toluene', 'n-undecane']
    return(df, sol)

df, sol = read_and_prepare()
print(df) 

def fillImputed(df):
    loc1 = df.loc[:, 'neg-octlR_log':'neg_undecR_log']
    loc2 = df.loc[:, 'predicted_1-octanol':'predicted_n-undecane']
    loc1.columns = sol
    loc2.columns = sol
    locC = loc1.combine_first(loc2)
    locC.columns = 'imputed_' + locC.columns
    print(locC)
    df = pd.concat([df, locC], axis=1)
    predCol = sorted(df.loc[:, 'imputed_1-octanol':'imputed_n-undecane'].columns)
    df = df.dropna(subset=predCol, axis=0)
    df = df.reset_index(drop=True)
    return(df)
    
df = fillImputed(df)
print(df)

def standardize_1(df):
    scaler = StandardScaler()
    loc1 = df.loc[:, 'octl':'undec'] 
    loc2 = df.loc[:, 'imputed_1-octanol':'imputed_n-undecane']
    def standardize_2(loc):
        loc = loc.T
        scaler.fit(loc) # print(scaler.mean_) to double check the standardization axis
        loc = scaler.transform(loc)
        loc = loc.T
        return(loc)
    loc1 = standardize_2(loc1)
    loc2 = standardize_2(loc2)
    loc1 = pd.DataFrame(loc1, columns = sol)
    loc2 = pd.DataFrame(loc2, columns = sol)
    loc1 = loc1.add_suffix('_db_std')
    loc2 = loc2.add_suffix('_exp_std')
    return(loc1, loc2)

loc1, loc2 = standardize_1(df)

df = pd.concat([df, loc1, loc2], axis=1)
print(loc1, loc2)
print(df)

scaler = StandardScaler()
loc1 = df.loc[:, 'octl':'undec']
oc2 = df.loc[:, '1-octanol_db_std':'n-undecane_db_std']

loc1 = loc1.T
scaler.fit(loc1)
loc1 = scaler.transform(loc1)
loc1 = loc1.T

loci = loc2.T
loci = scaler.inverse_transform(loci)
loci = loci.T
loci = pd.DataFrame(loci, columns = sol)
loci = loci.add_suffix('_exp_rvstd')

df = pd.concat([df, loci], axis=1)

def correlations(df):
    loc1 = df.loc[:, '1-octanol_db_std':'n-undecane_db_std']
    loc2 = df.loc[:, '1-octanol_exp_std':'n-undecane_exp_std']
    loc1.columns = loc1.columns.str.replace('_db_std', '')
    loc2.columns = loc2.columns.str.replace('_exp_std','')
    loc1 = loc1.T
    loc2 = loc2.T
    d1 = loc1.corrwith(loc2)
    d1 = d1**2
    df['Correlation_R2'] = d1 
    return(df, d1)

df, d1 = correlations(df)
print(d1)
print(df)

df.to_csv('1stExpMergedCleanNeg_imputed_rv.csv')






































