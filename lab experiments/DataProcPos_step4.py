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
    df1 = pd.read_csv('1stExpMergedCleanpos1.0_0.1min_0.00145Da.csv')
    df1 = df1.drop('Unnamed: 0', axis=1)
    df1['chem_id'] = (np.round(df1['Average Mz'], 3).astype(str) + 
                    '@' + np.round(df1['Average Rt(min)'], 3).astype(str))
    df2 = pd.read_csv('pos_combined_csv.csv')
    df = pd.merge(df1, df2, on='chem_id', how='left')
    sol = ['1-octanol', 'butyl acetate', 'chloroform', 'cyclohexane',
           'dichloromethane', 'triolein', 'n-hexane', 'n-octane', 'oleylalcohol',
           'toluene', 'n-undecane']
    return(df, sol)

df, sol = read_and_prepare()
print(df) 

def fillImputed(df):
    loc1 = df.loc[:, 'pos-octlR_log':'pos_undecR_log']
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

df.to_csv('1stExpMergedCleanpos1.0_0.1min_0.00145Da_imputed.csv')



def fig1(df):
    over15 = df[df['W_Rstd_flag'] == 0]['Correlation_R2']
    under15 = df[df['W_Rstd_flag'] == 1]['Correlation_R2']
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.boxplot([over15, under15], 
               medianprops = dict(linestyle='-', linewidth=1, color='blue'),
               #patch_artist=True,  # fill with color,
               #boxprops=dict(color='grey', alpha =0.4),
               labels=['over 15%', 'under 15%'])
    
    ax.set_ylabel('Correlation R2')
    ax.set_xlabel('RSTD')
    return(ax)

ax = fig1(df)
ax.figure.savefig('R2_Vs_WRstdflag.png', dpi=300)


def fig2(df):
    z0 = df[df['zeros'] == 0]['Correlation_R2']
    z1 = df[df['zeros'] == 1]['Correlation_R2']
    z2 = df[df['zeros'] == 2]['Correlation_R2']
    z3 = df[df['zeros'] == 3]['Correlation_R2']
    z4 = df[df['zeros'] == 4]['Correlation_R2']
    z5 = df[df['zeros'] == 5]['Correlation_R2']
    z6 = df[df['zeros'] == 6]['Correlation_R2']
    z7 = df[df['zeros'] == 7]['Correlation_R2']
    z8 = df[df['zeros'] == 8]['Correlation_R2']
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.boxplot([z0, z1, z2, z3, z4, z5, z6, z7, z8], 
               medianprops = dict(linestyle='-', linewidth=1, color='blue'),
               #patch_artist=True,  # fill with color,
               #boxprops=dict(color='grey', alpha =0.4),
               flierprops=dict(marker='.', markerfacecolor='blue', markeredgecolor='blue', alpha=0.9), 
               labels=['0', '1', '2', '3', '4', '5', '6', '7', '8'])
    
    ax.set_ylabel('Correlation R2')
    ax.set_xlabel('Zeros')
    return(ax)

ax = fig2(df)
ax.figure.savefig('R2_Vs_zeros.png', dpi=300)












































'''
# Combine first example code
df1 = pd.DataFrame({'A': [None, 0], 'B': [None, 4]})
df2 = pd.DataFrame({'A': [1, 1], 'B': [3, 3]})
df1.combine_first(df2)
'''




'''
fig = sns.boxplot(x='W_Rstd_flag', y='Correlation_R2',
                 data=df, palette='Paired')

fig.figure.savefig('R2_Vs_WRstdflag.png', dpi=300)


fig = sns.boxplot(x='zeros', y='Correlation_R2',
                 data=df, color='darkgrey')
plt.markers("o") 

'''